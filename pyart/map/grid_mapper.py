"""
pyart.map.grid_mapper
=====================

Utilities for mapping radar objects to Cartesian grids.

.. autosummary::
    :toctree: generated/

    grid_from_radars
    map_to_grid
    example_roi_func_constant
    example_roi_func_dist
    _load_nn_field_data
    _gen_roi_func_constant
    _gen_roi_func_dist
    _gen_roi_func_dist_beam

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    NNLocator

"""

import numpy as np
import scipy.spatial

from mpl_toolkits.basemap import pyproj

from ..config import get_fillvalue, get_field_name
from ..io.common import radar_coords_to_cart
from ..io.grid import Grid
from ._load_nn_field_data import _load_nn_field_data
from .ckdtree import cKDTree
from .ball_tree import BallTree


def grid_from_radars(radars, grid_shape, grid_dimensions, **kwargs):
    """
    Map one or more radars to a Cartesian grid, returning a Grid object.

    Additional arguments are passed to :py:func:`map_to_grid`

    Parameters
    ----------
    radars : list
        List of Radar objects which will be mapped to a Cartesian grid.
    grid_shape : tuple
        The (nz, ny, nx) dimension lengths of the grid.
    grid_dimensions : tuple
        The (z, y, x) coordinates of the grid in meters. These can describe
        either a uniform or non-uniform grid.

    Returns
    -------
    grid : Grid
        A :py:class:`pyart.io.Grid` object containing the gridded radar
        data.

    See Also
    --------
    map_to_grid : Map to grid and return a dictionary of radar fields

    """
    # map the radar(s) to a cartesian grid
    grids = map_to_grid(radars, grid_shape, grid_dimensions, **kwargs)

    # create and populate the field dictionary
    fields = {}
    first_radar = radars[0]

    for field in grids.keys():
        if field == 'radius_of_influence':
            fields['radius_of_influence'] = {
                'data': grids['radius_of_influence'],
                'standard_name': 'radius_of_influence',
                'long_name': 'Cressman radius of influence',
                'units': 'meters',
                'least_significant_digit': 1,
                '_FillValue': get_fillvalue()}
        else:
            fields[field] = {'data': grids[field]}
            # copy the metadata from the radar to the grid
            for key in first_radar.fields[field].keys():
                if key == 'data':
                    continue
                fields[field][key] = first_radar.fields[field][key]

    # time dictionaries
    time = {
        'data': np.array([first_radar.time['data'][0]]),
        'units': first_radar.time['units'],
        'calendar': first_radar.time['calendar'],
        'standard_name': first_radar.time['standard_name'],
        'long_name': 'Time in seconds since volume start'}

    time_start = {
        'data': np.array([first_radar.time['data'][0]]),
        'units': first_radar.time['units'],
        'calendar': first_radar.time['calendar'],
        'standard_name': first_radar.time['standard_name'],
        'long_name': 'Time in seconds of volume start'}

    time_end = {
        'data': np.array([first_radar.time['data'][-1]]),
        'units': first_radar.time['units'],
        'calendar': first_radar.time['calendar'],
        'standard_name': first_radar.time['standard_name'],
        'long_name': 'Time in seconds of volume end'}

    # grid coordinate dictionaries
    nz, ny, nx = grid_shape
    z, y, x = grid_dimensions

    xaxis = {'data':  x,
             'long_name': 'X-coordinate in Cartesian system',
             'axis': 'X',
             'units': 'meters'}

    yaxis = {'data': y,
             'long_name': 'Y-coordinate in Cartesian system',
             'axis': 'Y',
             'units': 'meters'}

    zaxis = {'data': z,
             'long_name': 'Z-coordinate in Cartesian system',
             'axis': 'Z',
             'units': 'meters',
             'positive': 'up'}

    # grid origin location dictionaries
    if 'grid_origin' in kwargs:
        lat = np.array([kwargs['grid_origin'][0]])
        lon = np.array([kwargs['grid_origin'][1]])
        alt = first_radar.altitude['data']
    else:
        lat = first_radar.latitude['data']
        lon = first_radar.longitude['data']
        alt = first_radar.altitude['data']

    alt_origin = {'data': alt,
                  'long_name': 'Altitude at grid origin',
                 'units': 'm',
                 'standard_name': 'altitude',
                 }

    lat_origin = {'data': lat,
                  'long_name': 'Latitude at grid origin',
                  'units': 'degrees_north',
                  'standard_name': 'latitude',
                  'valid_min': -90.,
                  'valid_max': 90.
                 }

    lon_origin = {'data': lon,
                  'long_name': 'Longitude at grid origin',
                  'units': 'degrees_east',
                  'standard_name': 'longitude',
                  'valid_min': -180.,
                  'valid_max': 180.
                 }

    # axes dictionary
    axes = {'time': time,
            'time_start': time_start,
            'time_end': time_end,
            'z_disp': zaxis,
            'y_disp': yaxis,
            'x_disp': xaxis,
            'alt': alt_origin,
            'lat': lat_origin,
            'lon': lon_origin}

    # metadata dictionary
    metadata = dict(first_radar.metadata)

    # add radar_{0,1,...}_{lat, lon, alt, instrument_name} key-value pairs
    # to the metadata dictionary.
    for i, radar in enumerate(radars):
        # TODO: add logic here to support moving platform radars
        metadata['radar_{0:d}_lat'.format(i)] = radar.latitude['data'][0]
        metadata['radar_{0:d}_lon'.format(i)] = radar.longitude['data'][0]
        metadata['radar_{0:d}_alt'.format(i)] = radar.altitude['data'][0]
        if 'instrument_name' in radar.metadata:
            i_name = radar.metadata['instrument_name']
        else:
            i_name = ''
        metadata['radar_{0:d}_instrument_name'.format(i)] = i_name

    return Grid(fields, axes, metadata)


class NNLocator:
    """
    Nearest neighbor locator.

    Class for finding the neighbors of a points within a given distance.

    Parameters
    ----------
    data : array_like, (n_sample, n_dimensions)
        Locations of points to be indexed.  Note that if data is a
        C-contiguous array of dtype float64 the data will not be copied.
        Otherwise and internal copy will be made.
    leafsize : int
        The number of points at which the algorithm switches over to
        brute-force.  This can significantly impact the speed of the
        construction and query of the tree.
    algorithm : 'kd_tree' or 'ball_tree'
        Algorithm used to compute the nearest neighbors. 'kd_tree' uses a
        k-d tree, 'ball_tree' a Ball tree.

    """

    def __init__(self, data, leafsize=10, algorithm='kd_tree'):
        """ Initialize. """
        self._algorithm = algorithm

        # build the query tree
        if algorithm == 'kd_tree':
            self.tree = cKDTree(data, leafsize=leafsize)
        elif algorithm == 'ball_tree':
            self.tree = BallTree(data, leaf_size=leafsize)

    def find_neighbors_and_dists(self, q, r):
        """
        Find all neighbors and distances within a given distance.

        Parameters
        ----------
        q : n_dimensional tuple
            Point to query
        r : float
            Distance within which neighbors are returned.

        Returns
        -------
        ind : array of integers
            Indices of the neighbors.
        dist : array of floats
            Distances to the neighbors.

        """
        if self._algorithm == 'kd_tree':
            ind = self.tree.query_ball_point(q, r)
            if len(ind) == 0:
                return ind, 0
            dist = scipy.spatial.minkowski_distance(q, self.tree.data[ind])
            return ind, dist

        elif self._algorithm == 'ball_tree':
            ind, dist = self.tree.query_radius(q, r, return_distance=True)
            return ind[0], dist[0]


def map_to_grid(radars, grid_shape, grid_dimensions, grid_origin=None,
                fields=None, refl_filter_flag=True, refl_field=None,
                max_refl=None, map_roi=True, weighting_function='Cressman',
                toa=17000.0, copy_field_data=True, algorithm='kd_tree',
                leafsize=10.0, roi_func='dist_beam', constant_roi=500.0,
                z_factor=0.05, xy_factor=0.02, min_radius=500.0,
                h_factor=1.0, nb=1.5, bsp=1.0, max_radius=4000.0,
                kappa_star=0.5, data_spacing=1220.0, proj='lcc',
                datum='NAD83', ellps='GRS80'):
    """
    Map one or more radars to a Cartesian grid.

    Generate a Cartesian grid of points for the requested fields from the
    collected points from one or more radars.  The field value for a grid
    point is found by interpolating from the collected points within a given
    radius of influence and weighting these nearby points according to their
    distance from the grid points. Collected points are filtered
    according to a number of criteria so that undesired points are not
    included in the interpolation.

    Parameters
    ----------
    radars : list
        List of Radar objects which will be mapped to a Cartesian grid.
    grid_shape : tuple
        The (nz, ny, nx) dimension lengths of the grid.
    grid_dimensions : tuple
        The (z, y, x) coordinates of the grid in meters. These can describe
        either a uniform or non-uniform grid.
    grid_origin : tuple or None
        Latitude and longitude of the grid origin. A value of None sets the
        origin to the location of the first radar.
    fields : list or None
        List of fields within the radar objects which will be mapped to
        the cartesian grid. None, the default, will map the fields which are
        present in all the radar objects.
    refl_filter_flag : bool
        True to filter the collected points based on the reflectivity field.
        False to perform no filtering.  Gates where the reflectivity field,
        specified by the `refl_field` parameter, is not-finited, masked or
        has a value above the `max_refl` parameter are excluded from the
        grid interpolation.
    refl_field : str
        Name of the field which will be used to filter the collected points.
        A value of None will use the default field name as defined in the
        Py-ART configuration file.
    max_refl : float
        Maximum allowable reflectivity.  Points in the `refl_field` which are
        above is value are not included in the interpolation. None will
        include skip this filtering.
    roi_func : str or function
        Radius of influence function.  A functions which takes an
        z, y, x grid location, in meters, and returns a radius (in meters)
        within which all collected points will be included in the weighting
        for that grid points. Examples can be found in the
        :py:func:`example_roi_func_constant`,
        :py:func:`example_roi_func_dist`, and
        :py:func:`example_roi_func_dist_beam`.
        Alternatively the following strings can use to specify a built in
        radius of influence function:

            * constant: constant radius of influence.
            * dist: radius grows with the distance from each radar.
            * dist_beam: radius grows with the distance from each radar
              and parameter are based of virtual beam sizes.

        The parameters which control these functions are listed in the
        `Other Parameters` section below.
    map_roi : bool
        True to include a radius of influence field in the returned
        dictionary under the 'ROI' key.  This is the value of roi_func at all
        grid points.
    weighting_function : 'Barnes' or 'Cressman'
        Functions used to weight nearby collected points when interpolating a
        grid point.
    toa : float
        Top of atmosphere in meters. Collected points above this height are
        not included in the interpolation.

    Other Parameters
    ----------------
    constant_roi : float
        Radius of influence parameter for the built in 'constant' function.
        This parameter is the constant radius in meter for all grid points.
        This parameter is only used when `roi_func` is `constant`.
    z_factor, xy_factor, min_radius : float
        Radius of influence parameters for the built in 'dist' function.
        The parameter correspond to the radius size increase, in meters,
        per meter increase in the z-dimension from the nearest radar,
        the same foreach meteter in the xy-distance from the nearest radar,
        and the minimum radius of influence in meters. These parameters are
        only used when `roi_func` is 'dist'.
    h_factor, nb, bsp, min_radius : float
        Radius of influence parameters for the built in 'dist_beam' function.
        The parameter correspond to the height scaling, virtual beam width,
        virtual beam spacing, and minimum radius of influence. These
        parameters are only used when `roi_func` is 'dist_mean'.
    max_radius : float
        The maximum radius in meters to search for points. This is only valid
        when weighting_function is 'Barnes'. 
    copy_field_data : bool
        True to copy the data within the radar fields for faster gridding,
        the dtype for all fields in the grid will be float64. False will not
        copy the data which preserves the dtype of the fields in the grid,
        may use less memory but results in significantly slower gridding
        times.  When False gates which are masked in a particular field but
        are not masked in the  `refl_field` field will still be included in
        the interpolation.  This can be prevented by setting this parameter
        to True or by gridding each field individually setting the
        `refl_field` parameter and the `fields` parameter to the field in
        question.  It is recommended to set this parameter to True.
    algorithm : 'kd_tree' or 'ball_tree'
        Algorithms to use for finding the nearest neighbors. 'kd_tree' tends
        to be faster.  This value should only effects the speed of the
        gridding, not the results.
    leafsize : int
        Leaf size passed to the neighbor lookup tree. This can affect the
        speed of the construction and query, as well as the memory required
        to store the tree. The optimal value depends on the nature of the
        problem. This value should only effect the speed of the gridding,
        not the results.
    proj, datum, ellps : str
        Parameters defining the map projection of the analysis domain onto
        Earth's surface. See the pyproj documentation for more information.

    Returns
    -------
    grids : dict
        Dictionary of mapped fields.  The keys of the dictionary are given by
        parameter fields.  Each element is a `grid_size` float64 array
        containing the interpolated grid for that field.

    See Also
    --------
    grid_from_radars : Map to grid and return a Grid object.

    """
    # check the parameters
    if weighting_function.upper() not in ['CRESSMAN', 'BARNES']:
        raise ValueError('Unknown weighting_function: %s' %weighting_function)
    if algorithm.upper() not in ['KD_TREE', 'BALL_TREE']:
        raise ValueError('Unknown algorithm: %s' %algorithm)
    fill_value = get_fillvalue()
    
    # determine what weighting scheme is desired since this affects whether
    # a radius of influence or smoothing parameter are used
    if weighting_function.upper() == 'CRESSMAN':
        is_cressman = True
    else:
        is_cressman = False
    if weighting_function.upper() == 'BARNES':
        is_barnes = True
    else:
        is_barnes = False

    # find the grid origin if not given
    if grid_origin is None:
        lat = float(radars[0].latitude['data'])
        lon = float(radars[0].longitude['data'])
        grid_origin = (lat, lon)

    # fields which should be mapped
    # if None is specified, fields which are available for all radars are
    # mapped
    if fields is None:
        fields = set(radars[0].fields.keys())
        for radar in radars[1:]:
            fields = fields.intersection(radar.fields.keys())
        fields = list(fields)
    nfields = len(fields)

    # determine the number of gates (collected points) for each radar
    nradars = len(radars)
    ngates_per_radar = [r.fields[fields[0]]['data'].size for r in radars]
    total_gates = np.sum(ngates_per_radar)
    gate_offset = np.cumsum([0] + ngates_per_radar)

    # create arrays to hold the gate locations and indicators if the gate
    # should be included in the interpolation
    gate_locations = np.empty((total_gates, 3), dtype=np.float64)
    include_gate = np.ones(total_gates, dtype=np.bool)
    
    # offsets from the grid origin in meters for each radar
    offsets = []

    # create the field lookup tables
    if copy_field_data:
        # copy_field_data == True, lookups are performed on a 2-D copy of
        # all of the field data in all radar objects, which can be a
        # large array. These lookups are fast since the data type is known
        field_data = np.ma.empty((total_gates, nfields), dtype=np.float64)
    else:
        # copy_field_data == False, lookups are performed on a 2-D array
        # object pointing to the radar fields themselves, no copies are made.
        # A table mapping filtered gates to raw gates is created later.
        # Since the data type is not not known, this method is generally
        # slower
        field_data_objs = np.empty((nfields, nradars), dtype='object')
        # We also need to know how many gates from each radar are included
        # in the NNLocator, the filtered_gates_per_radar list records this
        filtered_gates_per_radar = []

    # loop over the radars, calculating the Cartesian coordinates of all gate
    # locations and the offset, and acquire the field data
    for iradar, radar in enumerate(radars):

        # calculate radar offset from the origin
        radar_lat = float(radar.latitude['data'])
        radar_lon = float(radar.longitude['data'])
        proj = pyproj.Proj(proj=proj, datum=datum, ellps=ellps,
                           lat_0=grid_origin[0], lon_0=grid_origin[1],
                           x_0=0.0, y_0=0.0)
        radar_x, radar_y = proj(radar_lon, radar_lat)
        offsets.append((0, radar_y, radar_x)) # XXX should include z

        # calculate the Cartesian locations of all the radar gates
        rg, azg = np.meshgrid(radar.range['data'], radar.azimuth['data'])
        rg, eleg = np.meshgrid(radar.range['data'], radar.elevation['data'])
        xg_loc, yg_loc, zg_loc = radar_coords_to_cart(rg / 1000.0, azg, eleg)
        del rg, azg, eleg

        # add gate locations to gate_locations array
        start, end = gate_offset[iradar], gate_offset[iradar + 1]
        gate_locations[start:end, 0] = zg_loc.ravel()
        gate_locations[start:end, 1] = (yg_loc + radar_y).ravel()
        gate_locations[start:end, 2] = (xg_loc + radar_x).ravel()
        del xg_loc, yg_loc

        # determine which gates should be included in the interpolation
        # this involves removing gates that are above the TOA and where
        # there is missing reflectivity data
        gflags = zg_loc < toa

        if refl_filter_flag:
            if refl_field is None:
                refl_field = get_field_name('corrected_reflectivity')
            refl_data = radar.fields[refl_field]['data']
            
            if np.ma.is_masked(refl_data):
                gflags = np.logical_and(gflags, ~refl_data.mask)
                
            gflags = np.logical_and(gflags, np.isfinite(refl_data))

            if max_refl is not None:
                gflags = np.logical_and(gflags, refl_data < max_refl)
                
            if isinstance(gflags, np.ma.MaskedArray):
                gflags = gflags.data

            del refl_data, zg_loc
            
        include_gate[start:end] = gflags.ravel()

        if not copy_field_data:
            # record the number of gates from the current radar which will be
            # included in the interpolation
            filtered_gates_per_radar.append(gflags.sum())

        del gflags

        # copy and/or store references to field data for lookup
        for ifield, field in enumerate(fields):
            flat_field_data = radar.fields[field]['data'].ravel()
            if copy_field_data:
                field_data[start:end, ifield] = flat_field_data
            else:
                field_data_objs[ifield, iradar] = flat_field_data
        del flat_field_data

    # build field data lookup tables
    if copy_field_data:
        # copy_field_data == True we filter the field data in the
        # same manner as we filtered the gate locations
        filtered_field_data = field_data[include_gate]
    else:
        # copy_field_data == False, build a lookup table which maps from
        # filtered gate number to (radar number, radar gate number)
        # the radar number is given as the quotient of the lookup table
        # value divided by total_gates, the remainder gives the index of
        # the gate in the flat field data array.

        # initialize the lookup table with values from 0 ... total gates
        lookup = np.where(include_gate)[0]

        # number of filtered gates before a given radar
        filtered_gate_offset = np.cumsum([0] + filtered_gates_per_radar)

        # for radars 1 to N-1 add total_gates to the lookup table and
        # subtract the number of gates in all earlier radars
        for i in xrange(1, nradars):
            l_start = filtered_gate_offset[i]
            l_end = filtered_gate_offset[i + 1]
            gates_before = gate_offset[i]
            lookup[l_start:l_end] += (total_gates * i - gates_before)
            
    # unpack the analysis domain parameters
    nz, ny, nx = grid_shape
    z, y, x = grid_dimensions

    # populate the nearest neighbor locator with the filtered gate locations
    nnlocator = NNLocator(gate_locations[include_gate], algorithm=algorithm,
                          leafsize=leafsize)

    if not hasattr(roi_func, '__call__') and is_cressman:
        if roi_func == 'constant':
            roi_func = _gen_roi_func_constant(constant_roi)
        elif roi_func == 'dist':
            roi_func = _gen_roi_func_dist(
                z_factor, xy_factor, min_radius, offsets)
        elif roi_func == 'dist_beam':
            roi_func = _gen_roi_func_dist_beam(
                h_factor, nb, bsp, min_radius, offsets)
        else:
            raise ValueError('Unknown Cressman roi_func: %s' %roi_func)

    # create array to hold interpolated grid data and radius of influence if
    # requested
    grid_data = np.ma.empty((nz, ny, nx, nfields), dtype=np.float64)
    grid_data.set_fill_value(fill_value)

    if map_roi and is_cressman:
        roi = np.empty((nz, ny, nx), dtype=np.float64)

    # interpolate field values onto each analysis grid point
    for iz, iy, ix in np.ndindex(nz, ny, nx):
        
        # determine either the radius of influence (Cressman) or the maximum
        # radius (Barnes) to check for nearest neighbors to the current
        # analysis point
        # find nearest neighbors and their distances
        if is_cressman:
            Rc = roi_func(z[iz], y[iy], x[ix])
            if map_roi:
                roi[iz, iy, ix] = Rc
            
            ind, dist = nnlocator.find_neighbors_and_dists(
                                        (z[iz], y[iy], x[ix]), Rc)
        
        if is_barnes:
            ind, dist = nnlocator.find_neighbors_and_dists(
                                        (z[iz], y[iy], x[ix]), max_radius)
            
        # when there are no neighbors, mark the grid point as bad
        if len(ind) == 0:
            grid_data[iz, iy, ix] = np.ma.masked
            grid_data.data[iz, iy, ix] = fill_value
            continue

        # find the field values for all neighbors
        if copy_field_data:
            # copy_field_data == True, a slice will get the field data
            nn_field_data = filtered_field_data[ind]
        else:
            # copy_field_data == False, use the lookup table to find the
            # radar numbers and gate numbers for the neighbors.  Then
            # use the _load_nn_field_data function to load this data from
            # the field data object array.  This is done in Cython for
            # speed
            r_nums, e_nums = divmod(lookup[ind], total_gates)
            npoints = r_nums.size
            nn_field_data = np.empty((npoints, nfields), np.float64)
            _load_nn_field_data(field_data_objs, nfields, npoints, r_nums,
                                e_nums, nn_field_data)
        
        # perform weighting of neighbors
        dist2 = dist**2
        
        # the Cressman filter (weight) is a function of the radial
        # distance separating an analysis point from a data point
        # (dist) and the radius of influence (Rc) parameter
        if is_cressman:
            Rc2 = Rc**2
            wq = (Rc2 - dist2) / (Rc2 + dist2)
            fq = np.ma.average(nn_field_data, weights=wq, axis=0)
        
        # the Barnes filter (weight) is a function of the radial distance
        # separating an analysis point from a data point (dist) and the
        # smoothing parameter (kappa)
        if is_barnes:
            kappa = kappa_star * (2.0 * data_spacing)**2
            wq = np.exp(-dist2 / kappa)
            fq = np.ma.average(nn_field_data, weights=wq, axis=0)

        grid_data[iz, iy, ix] = fq

    # create and return the grid dictionary
    grids = dict([(f, grid_data[..., i]) for i, f in enumerate(fields)])
    
    if map_roi and is_cressman:
        grids['radius_of_influence'] = roi
    
    return grids


# Radius of influence (RoI) functions

def example_roi_func_constant(zg, yg, xg):
    """
    Example RoI function which returns a constant radius.

    Parameters
    ----------
    zg, yg, xg : float
        Distance from the grid center in meters for the x, y and z axes.

    Returns
    -------
    roi : float
        Radius of influence in meters
    """
    return 500.0 # constant 500 meter radius of influence


def _gen_roi_func_constant(constant_roi):
    """
    Return a RoI function which returns a constant radius.

    See :py:func:`map_to_grid` for a description of the parameters.
    """

    def roi(zg, yg, xg):
        """ constant radius of influence function. """
        return constant_roi

    return roi


def example_roi_func_dist(zg, yg, xg):
    """
    Example RoI function which returns a radius which grows with distance.

    Parameters
    ----------
    zg, yg, xg : float
        Distance from the grid center in meters for the x, y and z axes.

    Returns
    -------
    roi : float

    """
    # RoI function parameters
    z_factor = 0.05         # increase in radius per meter increase in z dim
    xy_factor = 0.02        # increase in radius per meter increase in xy dim
    min_radius = 500.0      # minimum radius
    offsets = ((0, 0, 0), )  # z, y, x offset of grid in meters from radar(s)

    offsets = np.array(offsets)
    zg_off = offsets[:, 0]
    yg_off = offsets[:, 1]
    xg_off = offsets[:, 2]
    r = (z_factor * (zg - zg_off) +
         xy_factor * np.sqrt((xg - xg_off)**2 + (yg - yg_off)**2) +
         min_radius)
    return min(r)


def _gen_roi_func_dist(z_factor, xy_factor, min_radius, offsets):
    """
    Return a RoI function whose radius grows with distance.

    See :py:func:`map_to_grid` for a description of the parameters.
    """
    offsets = np.array(offsets)
    zg_off = offsets[:, 0]
    yg_off = offsets[:, 1]
    xg_off = offsets[:, 2]

    def roi(zg, yg, xg):
        """ dist radius of influence function. """
        r = (z_factor * (zg - zg_off) +
             xy_factor * np.sqrt((xg - xg_off)**2 + (yg - yg_off)**2) +
             min_radius)
        return min(r)

    return roi


def example_roi_func_dist_beam(zg, yg, xg):
    """
    Example RoI function which returns a radius which grows with distance
    and whose parameters are based on virtual beam size.

    Parameters
    ----------
    zg, yg, xg : float
        Distance from the grid center in meters for the x, y and z axes.

    Returns
    -------
    roi : float

    """
    # RoI function parameters
    h_factor = 1.0      # height scaling
    nb = 1.5            # virtual beam width
    bsp = 1.0           # virtual beam spacing
    min_radius = 500.0  # minimum radius in meters
    offsets = ((0, 0, 0), )  # z, y, x offset of grid in meters from radar(s)

    offsets = np.array(offsets)
    zg_off = offsets[:, 0]
    yg_off = offsets[:, 1]
    xg_off = offsets[:, 2]
    r = (h_factor * ((zg - zg_off) / 20.0) +
         np.sqrt((yg - yg_off)**2 + (xg - xg_off)**2) *
         np.tan(nb * bsp * np.pi / 180.0) + min_radius)
    return min(r)


def _gen_roi_func_dist_beam(h_factor, nb, bsp, min_radius, offsets):
    """
    Return a RoI function whose radius which grows with distance
    and whose parameters are based on virtual beam size.

    See :py:func:`map_to_grid` for a description of the parameters.
    """
    offsets = np.array(offsets)
    zg_off = offsets[:, 0]
    yg_off = offsets[:, 1]
    xg_off = offsets[:, 2]

    def roi(zg, yg, xg):
        """ dist_beam radius of influence function. """
        r = (h_factor * ((zg - zg_off) / 20.0) +
             np.sqrt((yg - yg_off)**2 + (xg - xg_off)**2) *
             np.tan(nb * bsp * np.pi / 180.0) + min_radius)
        return min(r)

    return roi
