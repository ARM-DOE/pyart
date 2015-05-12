"""
pyart.map.gates_to_grid
=======================

Generate a Cartesian grid by mapping from radar gates onto the grid.

.. autosummary::
    :toctree: generated/

    map_gates_to_grid
    _detemine_cy_weighting_func
    _parse_grid_origin
    _determine_fields
    _find_offsets
    _find_grid_params
    _parse_roi_func

"""

import numpy as np
from ..config import get_field_name
from ..graph.common import corner_to_point

from ._gate_to_grid_map import GateToGridMapper
from ._gate_to_grid_map import RoIFunction, ConstantRoI, DistBeamRoI, DistRoI


def map_gates_to_grid(
        radars, grid_shape, grid_limits, grid_origin=None,
        grid_origin_alt=None, fields=None, refl_filter_flag=True,
        refl_field=None, max_refl=None, map_roi=True,
        weighting_function='Barnes', toa=17000.0, roi_func='dist_beam',
        constant_roi=500., z_factor=0.05, xy_factor=0.02, min_radius=500.0,
        h_factor=1.0, nb=1.5, bsp=1.0):
    """
    Map gates from one or more radars to a Cartesian grid.

    Generate a Cartesian grid of points for the requested fields from the
    collected points from one or more radars. For each radar gate that is not
    filtered a radius of influence is calculated.  The weighted field values
    for that gate are added to all grid points within that radius.  This
    routine scaled linearly with the number of radar gates and the effective
    grid size.

    Parameters not defined below are identical to those in
    :py:func:`map_to_grid`.

    Parameters
    ----------
    roi_func : str or RoIFunction
        Radius of influence function.  A functions which takes an
        z, y, x grid location, in meters, and returns a radius (in meters)
        within which all collected points will be included in the weighting
        for that grid points. Examples can be found in the
        Typically following strings can use to specify a built in
        radius of influence function:

            * constant: constant radius of influence.
            * dist: radius grows with the distance from each radar.
            * dist_beam: radius grows with the distance from each radar
              and parameter are based of virtual beam sizes.

        A custom RoIFunction can be defined using the RoIFunction class
        and defining a get_roi method which returns the radius.  For efficient
        mapping this class should be implemented in Cython.

    Returns
    -------
    grids : dict
        Dictionary of mapped fields.  The keysof the dictionary are given by
        parameter fields.  Each elements is a `grid_size` float64 array
        containing the interpolated grid for that field.

    See Also
    --------
    grid_from_radars : Map to a grid and return a Grid object
    map_to_grid : Create grid by finding the radius of influence around each
                  grid point.

    """
    if max_refl is None:    # parse max_refl
        max_refl = np.finfo('float32').max
    if grid_origin_alt is None:
        grid_origin_alt = float(radars[0].altitude['data'])

    cy_weighting_function = _detemine_cy_weighting_func(weighting_function)
    grid_origin = _parse_grid_origin(grid_origin, radars)
    fields = _determine_fields(fields, radars, refl_field)
    offsets = _find_offsets(radars, grid_origin, grid_origin_alt)
    grid_starts, grid_steps = _find_grid_params(grid_shape, grid_limits)
    roi_func = _parse_roi_func(roi_func, constant_roi, z_factor, xy_factor,
                               min_radius, h_factor, nb, bsp, offsets)

    # prepare grid storage arrays
    nfields = len(fields)
    grid_sum = np.zeros(grid_shape + (nfields, ), dtype=np.float32)
    grid_wsum = np.zeros(grid_shape + (nfields, ), dtype=np.float32)
    gatemapper = GateToGridMapper(
        grid_shape, grid_starts, grid_steps, grid_sum, grid_wsum)

    # project gates from each radar onto the grid
    for radar, radar_offset in zip(radars, offsets):

        # Copy the field data and masks.
        # TODO method that does not copy field data into new array
        shape = (radar.nrays, radar.ngates, nfields)
        field_data = np.empty(shape, dtype='float32')
        field_mask = np.empty(shape, dtype='uint8')
        for i, field in enumerate(fields):
            fdata = radar.fields[field]['data']
            field_data[:, :, i] = np.ma.getdata(fdata)
            field_mask[:, :, i] = np.ma.getmaskarray(fdata)

        # map the gates onto the grid
        gatemapper.map_gates_to_grid(
            radar.elevation['data'].astype('float32'),
            radar.azimuth['data'].astype('float32'),
            radar.range['data'].astype('float32'),
            field_data, field_mask, radar_offset, toa, roi_func,
            refl_filter_flag, max_refl, cy_weighting_function)

    # create and return the grid dictionary
    mweight = np.ma.masked_equal(grid_wsum, 0)
    msum = np.ma.masked_array(grid_sum, mweight.mask)
    grids = dict(
        [(f, msum[..., i] / mweight[..., i]) for i, f in enumerate(fields)])
    if map_roi:
        roi_array = np.empty(grid_shape, dtype=np.float32)
        gatemapper.find_roi_for_grid(roi_array, roi_func)
        grids['ROI'] = roi_array
    return grids


def _detemine_cy_weighting_func(weighting_function):
    """ Determine cython weight function value. """
    if weighting_function.upper() == 'CRESSMAN':
        cy_weighting_function = 1
    elif weighting_function.upper() == 'BARNES':
        cy_weighting_function = 0
    else:
        raise ValueError('unknown weighting_function')
    return cy_weighting_function


def _parse_grid_origin(grid_origin, radars):
    """ Parse the grid origin parameter, finding origin if not given. """
    if grid_origin is None:
        lat = float(radars[0].latitude['data'])
        lon = float(radars[0].longitude['data'])
        grid_origin = (lat, lon)
    return grid_origin


def _determine_fields(fields, radars, refl_field):
    """ Determine which field should be mapped to the grid. """
    if fields is None:
        fields = set(radars[0].fields.keys())
        for radar in radars[1:]:
            fields = fields.intersection(radar.fields.keys())
        fields = list(fields)

    # find the reflectivity field, check that it is mapped and
    # move it to the front of the fields list
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if refl_field not in fields:
        raise ValueError('reflectivity field not mapped')
    fields.insert(0, fields.pop(fields.index(refl_field)))
    return fields


def _find_offsets(radars, grid_origin, grid_origin_alt):
    """ Find offset between radars and grid origin. """
    # loop over the radars finding offsets from the origin
    offsets = []    # offsets from the grid origin, in meters, for each radar
    for radar in radars:
        radar_lat = float(radar.latitude['data'])
        radar_lon = float(radar.longitude['data'])
        x_disp, y_disp = corner_to_point(grid_origin, (radar_lat, radar_lon))
        z_disp = float(radar.altitude['data']) - grid_origin_alt
        offsets.append((z_disp, y_disp, x_disp))
    return offsets


def _find_grid_params(grid_shape, grid_limits):
    """ Find the starting points and step size of the grid. """
    nz, ny, nx = grid_shape
    zr, yr, xr = grid_limits
    z_start, z_stop = zr
    y_start, y_stop = yr
    x_start, x_stop = xr

    if nz == 1:
        z_step = 0.
    else:
        z_step = (z_stop - z_start) / (nz - 1.)
    if ny == 1:
        y_step = 0.
    else:
        y_step = (y_stop - y_start) / (ny - 1.)
    if nx == 1:
        x_step = 0.
    else:
        x_step = (x_stop - x_start) / (nx - 1.)

    grid_starts = (z_start, y_start, x_start)
    grid_steps = (z_step, y_step, x_step)
    return grid_starts, grid_steps


def _parse_roi_func(roi_func, constant_roi, z_factor, xy_factor, min_radius,
                    h_factor, nb, bsp, offsets):
    """ Return the Radius of influence object. """
    if not isinstance(roi_func, RoIFunction):
        if roi_func == 'constant':
            roi_func = ConstantRoI(constant_roi)
        elif roi_func == 'dist':
            roi_func = DistRoI(z_factor, xy_factor, min_radius, offsets)
        elif roi_func == 'dist_beam':
            roi_func = DistBeamRoI(h_factor, nb, bsp, min_radius, offsets)
        else:
            raise ValueError('unknown roi_func: %s' % roi_func)
    return roi_func
