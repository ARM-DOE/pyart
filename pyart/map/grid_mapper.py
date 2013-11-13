"""
pyart.map.grid_mapper
=====================

Utilities for mapping radar objects to Cartesian grids.

.. autosummary::
    :toctree: generated/

    grid_from_radars
    map_to_grid
    _load_nn_field_data

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    NNLocator

"""

import numpy as np
import scipy.spatial

from ..config import get_fillvalue
from ..graph.common import corner_to_point
from ..io.common import radar_coords_to_cart
from ..io.grid import Grid
from ._load_nn_field_data import _load_nn_field_data
from .ckdtree import cKDTree
from .ball_tree import BallTree


def grid_from_radars(radars, grid_shape, grid_limits, **kwargs):
    """
    Map one or more radars to a Cartesian grid returning a Grid object.

    Additional arguments are passed to :py:func:`map_to_grid`

    Parameters
    ----------
    radars : tuple of Radar objects.
        Radar objects which will be mapped to the Cartesian grid.
    grid_shape : 3-tuple of floats
        Number of points in the grid (x, y, z).
    grid_limits : 3-tuple of 2-tuples
        Minimum and maximum grid location (inclusive) in meters for the
        x, y, z coordinates.

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
    grids = map_to_grid(radars, grid_shape=grid_shape,
                        grid_limits=grid_limits, **kwargs)

    # create and populate the field dictionary
    fields = {}
    first_radar = radars[0]

    for field in grids.keys():
        if field == 'ROI':
            fields['ROI'] = {
                'data': grids['ROI'],
                'standard_name': 'radius_of_influence',
                'long_name': 'Radius of influence for mapping',
                'units': 'm',
                'least_significant_digit': 1,
                'valid_min': 0.,
                'valid_max': 100000.,
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
    nx, ny, nz = grid_shape
    (x0, x1), (y0, y1), (z0, z1) = grid_limits

    xaxis = {'data':  np.linspace(x0, x1, nx),
             'long_name': 'X-coordinate in Cartesian system',
             'axis': 'X',
             'units': 'm'}

    yaxis = {'data': np.linspace(y0, y1, ny),
             'long_name': 'Y-coordinate in Cartesian system',
             'axis': 'Y',
             'units': 'm'}

    zaxis = {'data': np.linspace(z0, z1, nz),
             'long_name': 'Z-coordinate in Cartesian system',
             'axis': 'Z',
             'units': 'm',
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

    altorigin = {'data': alt,
                 'long_name': 'Altitude at grid origin',
                 'units': 'm',
                 'standard_name': 'altitude',
                 }

    latorigin = {'data': lat,
                 'long_name': 'Latitude at grid origin',
                 'units': 'degree_N',
                 'standard_name': 'latitude',
                 'valid_min': -90.,
                 'valid_max': 90.
                 }

    lonorigin = {'data': lon,
                 'long_name': 'Longitude at grid origin',
                 'units': 'degree_E',
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
            'alt': altorigin,
            'lat': latorigin,
            'lon': lonorigin}

    # metadata dictionary
    metadata = dict(first_radar.metadata)

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
        Othersize and internal copy will be made.
    leafsize : int
        The number of points at which the algorithm switches over to
        brute-force.  This can significantly impact the speed of the
        contruction and query of the tree.
    algorithm : 'kd_tree' or 'ball_tree'
        Algorithm used to compute the nearest neigbors.  'kd_tree' uses a
        k-d tree, 'ball_tree' a Ball tree.

    """

    def __init__(self, data, leafsize=10, algorithm='kd_tree'):
        """ initalize. """
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
        ind : array of intergers
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


def map_to_grid(radars, grid_shape=(81, 81, 69),
                grid_limits=((-30000., 20000), (-20000., 20000.), (0, 17000.)),
                grid_origin=None, fields=None,
                refl_field='reflectivity_horizontal', max_refl=190.0,
                qrf_func=None, map_roi=True, weighting_function='Barnes',
                toa=17000.0,
                h_factor=1.0, nb=1.5, bsp=1.0, min_radius=500.0,
                copy_field_data=True, algorithm='kd_tree', leafsize = 10):
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
    radars : tuple of Radar objects.
        Radar objects which will be mapped to the Cartesian grid.
    grid_shape : 3-tuple of floats
        Number of points in the grid (x, y, z).
    grid_limits : 3-tuple of 2-tuples
        Minimum and maximum grid location (inclusive) in meters for the
        x, y, z coordinates.
    grid_origin : (float, float) or None
        Latitude and longitude of grid origin.  None sets the origin
        to the location of the first radar.
    fields : list or None
        List of fields within the radar objects which will be mapped to
        the cartesian grid. None, the default, will map the fields which are
        present in all the radar objects.
    refl_field : str
        Name of the field which will be used to filter the collected points
        based on masking and maximum reflectivity.
    max_refl : float
        Maximum allowable reflectivity.  Points in the refl_field which are
        above is value are not included in the interpolation.
    qrf_func : function or None
        Query radius of influence function.  A functions which takes an
        x, y, z grid location, in meters, and returns a radius (in meters)
        within which all collected points will be included in the weighting
        for that grid points. None will use a function which takes into
        account the h_factor, nb and bsp parameters and increases the radius
        quadratically with elevation.
    map_roi : bool
        True to include a radius of influence field in the returned
        dictionary under the 'ROI' key.  This is the value of qrf_func at all
        grid points.
    weighting_function : 'Barnes' or 'Cressman'
        Functions used to weight nearby collected points when interpolating a
        grid point.
    toa : float
        Top of atmosphere in meters. Collected points above this height are
        not included in the interpolation.

    Other parameters
    ----------------
    h_factor : float
        H factor which influences the increase in the radius of influence as
        elevation increases.  Only used when qrf_func is None.
    nb : float
        Virtual beam width. Only used when qrf_func is None.
    bsp : float
        Virtual beam spacing. Only used when qrf_func is None.
    min_radius : float
        Minimum radius of influence. Only used when qrf_func is None.
    copy_field_data : bool
        True to copy the data within the radar fields for faster gridding,
        False will not copy the data which will use less memory but result in
        significantly fast gridding times.
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


    Returns
    -------
    grids : dict
        Dictionary of mapped fields.  The keysof the dictionary are given by
        parameter fields.  Each elements is a `grid_size` float64 array
        containing the interpolated grid for that field.

    See Also
    --------
    grid_from_radars : Map to grid and return a Grid object.

    """
    # check the parameters
    if weighting_function.upper() not in ['CRESSMAN', 'BARNES']:
        raise ValueError('unknown weighting_function')
    if algorithm not in ['kd_tree', 'ball_tree']:
        raise ValueError('unknow algorithm: %s' % algorithm)
    badval = get_fillvalue()

    # find the grid origin if not given
    if grid_origin is None:
        lat = float(radars[0].latitude['data'])
        lon = float(radars[0].longitude['data'])
        grid_origin = (lat, lon)

    # fields which should be mapped, None for fields which are in all radars
    if fields is None:
        fields = set(radars[0].fields.keys())
        for radar in radars[1:]:
            fields = fields.intersection(radar.fields.keys())
    nfields = len(fields)

    # determine the number of gates (collected points) in each radar
    nradars = len(radars)
    ngates_per_radar = [r.fields[refl_field]['data'].size for r in radars]
    total_gates = np.sum(ngates_per_radar)
    gate_offset = np.cumsum([0] + ngates_per_radar)

    # create arrays to hold the gate locations and indicators if the gate
    # should be included in the interpolation.
    gate_locations = np.empty((total_gates, 3), dtype=np.float64)
    include_gate = np.ones((total_gates), dtype=np.bool)

    # create a field lookup tables
    if copy_field_data:
        # copy_field_data == True, lookups are performed on a 2D copy of
        # all of the field data in all radar objects, this can be a
        # large array.  These lookup are fast as the dtype is know.
        field_data = np.empty((total_gates, nfields), dtype=np.float64)
    else:
        # copy_field_data == False, lookups are performed on a 2D object
        # array pointing to the radar fields themselved, no copies are made.
        # A table mapping filtered gates to raw gates is created later.
        # Since the dtype is not not know this method is slow.
        field_data_objs = np.empty((nfields, nradars), dtype='object')
        # We also need to know how many gates from each radar are included
        # in the NNLocator, the filtered_gates_per_radar list records this
        filtered_gates_per_radar = []

    # loop over the radars finding gate locations and field data
    for iradar, radar in enumerate(radars):

        # calculate radar offset from the origin
        radar_lat = float(radar.latitude['data'])
        radar_lon = float(radar.longitude['data'])
        x_disp, y_disp = corner_to_point(grid_origin, (radar_lat, radar_lon))

        # calculate cartesian locations of gates
        rg, azg = np.meshgrid(radar.range['data'], radar.azimuth['data'])
        rg, eleg = np.meshgrid(radar.range['data'], radar.elevation['data'])
        xg_loc, yg_loc, zg_loc = radar_coords_to_cart(rg / 1000., azg, eleg)
        del rg, azg, eleg

        # add gate locations to gate_locations array
        start, end = gate_offset[iradar], gate_offset[iradar + 1]
        gate_locations[start:end, 0] = (xg_loc + x_disp).flat
        gate_locations[start:end, 1] = (yg_loc + y_disp).flat
        gate_locations[start:end, 2] = zg_loc.flat
        del xg_loc, yg_loc

        # determine which gates should be included in the interpolation
        refl_data = radar.fields[refl_field]['data']
        gflags = zg_loc < toa      # include only those below toa
        gflags = np.logical_and(gflags, np.isfinite(refl_data))

        if max_refl is not None:
            gflags = np.logical_and(gflags, refl_data < max_refl)

        if np.ma.is_masked(refl_data):
            gflags = np.logical_and(gflags, np.logical_not(refl_data.mask))
            gflags = gflags.data
        include_gate[start:end] = gflags.flat

        if not copy_field_data:
            # record the number of gates from the current radar which
            # are included in the interpolation.
            filtered_gates_per_radar.append(gflags.sum())

        del refl_data, gflags, zg_loc

        # copy/store references to field data for lookup
        for ifield, field in enumerate(fields):
            flat_field_data = radar.fields[field]['data'].ravel()
            if copy_field_data:
                field_data[start:end, ifield] = flat_field_data
            else:
                field_data_objs[ifield, iradar] = flat_field_data
        del flat_field_data

    # build field data lookup tables
    if copy_field_data:
        # copy_field_data == True we filtered the field data in the
        # same manner as we will filter the gate locations.
        filtered_field_data = field_data[include_gate]
    else:
        # copy_field_data == True, build a lookup table which maps from
        # filtered gate number to (radar number, radar gate number)
        # the radar number is given as the quotent of the lookup table
        # value divided by total_gates, the remainer gives the index of
        # the gate in the flat field data array.

        # initalize the lookup table with values from 0 ... total gates
        lookup = np.where(include_gate)[0]

        # number of filtered gates before a given radar
        filtered_gate_offset = np.cumsum([0] + filtered_gates_per_radar)

        # for radars 1 to N-1 add total_gates to the lookup table and
        # subtract the number of gates in all ealier radars.
        for i in xrange(1, nradars):
            l_start = filtered_gate_offset[i]
            l_end = filtered_gate_offset[i + 1]
            gates_before = gate_offset[i]
            lookup[l_start:l_end] += (total_gates * i - gates_before)

    # populate the nearest neighbor locator with the filtered gate locations
    nnlocator = NNLocator(gate_locations[include_gate], algorithm=algorithm,
                          leafsize=leafsize)

    # unpack the grid parameters
    nx, ny, nz = grid_shape
    xr, yr, zr = grid_limits
    x_start, x_stop = xr
    y_start, y_stop = yr
    z_start, z_stop = zr
    x_step = (x_stop - x_start) / (nx - 1.)
    y_step = (y_stop - y_start) / (ny - 1.)
    z_step = (z_stop - z_start) / (nz - 1.)

    def standard_qrf(xg, yg, zg):
        return (h_factor * (zg / 20.0) + np.sqrt(yg ** 2 + xg ** 2) *
                np.tan(nb * bsp * np.pi / 180.0) + min_radius).flatten()

    if qrf_func is None:
        qrf_func = standard_qrf

    # create array to hold interpolated grid data and roi if requested
    grid_data = np.ma.empty((nz, ny, nx, nfields), dtype=np.float64)
    grid_data.set_fill_value(badval)

    if map_roi:
        roi = np.empty((nz, ny, nx), dtype=np.float64)

    # interpolate field values for each point in the grid
    for iz, iy, ix in np.ndindex(nz, ny, nx):

        # calculate the grid point
        x = x_start + x_step * ix
        y = y_start + y_step * iy
        z = z_start + z_step * iz
        r = qrf_func(x, y, z)
        if map_roi:
            roi[iz, iy, ix] = r

        #find neighbors and distances
        ind, dist = nnlocator.find_neighbors_and_dists((x, y, z), r)

        if len(ind) == 0:
            # when there are no neighbors, mark the grid point as bad
            grid_data[iz, iy, ix] = np.ma.masked
            grid_data.data[iz, iy, ix] = badval
            continue

        # find the field values for all neighbors
        if copy_field_data:
            # copy_field_data == True, a slice will get the field data.
            nn_field_data = filtered_field_data[ind]
        else:
            # copy_field_data == False, use the lookup table to find the
            # radar numbers and gate numbers for the neighbors.  Then
            # use the _load_nn_field_data function to load this data from
            # the field data object array.  This is done in Cython for speed.
            r_nums, e_nums = divmod(lookup[ind], total_gates)
            npoints = r_nums.size
            nn_field_data = np.empty((npoints, nfields), np.float64)
            _load_nn_field_data(field_data_objs, nfields, npoints, r_nums,
                                e_nums, nn_field_data)

        # preforms weighting of neighbors.
        dist2 = dist * dist
        r2 = r * r

        if weighting_function.upper() == 'CRESSMAN':
            weights = (r2 - dist2) / (r2 + dist2)
            value = np.average(nn_field_data, weights=weights, axis=0)
        elif weighting_function.upper() == 'BARNES':
            w = np.exp(-dist2 / (2.0 * r2)) + 1e-5
            w /= np.sum(w)
            value = np.ma.dot(w, nn_field_data)

        grid_data[iz, iy, ix] = value

    # create and return the grid dictionary
    grids = dict([(f, grid_data[..., i]) for i, f in enumerate(fields)])
    if map_roi:
        grids['ROI'] = roi
    return grids
