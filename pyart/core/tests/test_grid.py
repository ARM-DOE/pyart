""" Unit Tests for Py-ART's core/grid.py module. """

from __future__ import print_function

import functools
import warnings

import numpy as np
from numpy.testing import assert_almost_equal, assert_raises
from numpy.testing.decorators import skipif

try:
    from mpl_toolkits.basemap import pyproj
    _PYPROJ_AVAILABLE = True
except ImportError:
    _PYPROJ_AVAILABLE = False

import pyart
from pyart.lazydict import LazyLoadDict


def test_grid_write_method():
    grid1 = pyart.testing.make_target_grid()

    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_grid.nc'
        grid1.write(tmpfile)
        grid2 = pyart.io.read_grid(tmpfile)

        # check fields
        for field in grid1.fields.keys():
            _check_dicts_similar(grid1.fields[field], grid2.fields[field])

        # check attributes
        _check_dicts_similar(grid1.metadata, grid2.metadata)

        _check_dicts_similar(grid1.time, grid2.time)

        _check_dicts_similar(grid1.origin_latitude, grid2.origin_latitude)
        _check_dicts_similar(grid1.origin_longitude, grid2.origin_longitude)
        _check_dicts_similar(grid1.origin_altitude, grid2.origin_altitude)

        _check_dicts_similar(grid1.x, grid2.x)
        _check_dicts_similar(grid1.y, grid2.y)
        _check_dicts_similar(grid1.z, grid2.z)

        _check_dicts_similar(grid1.point_x, grid2.point_x)
        _check_dicts_similar(grid1.point_y, grid2.point_y)
        _check_dicts_similar(grid1.point_z, grid2.point_z)

        _check_dicts_similar(grid1.projection, grid2.projection)

        _check_dicts_similar(grid1.point_latitude, grid2.point_latitude)
        _check_dicts_similar(grid1.point_longitude, grid2.point_longitude)
        _check_dicts_similar(grid1.point_altitude, grid2.point_altitude)

        assert grid1.nx == grid2.nx
        assert grid1.ny == grid2.ny
        assert grid1.nz == grid2.nz

        # XXX not written at the moment
        # _check_dicts_similar(grid1.radar_latitude, grid2.radar_latitude)
        # _check_dicts_similar(grid1.radar_longitude, grid2.radar_longitude)
        # _check_dicts_similar(grid1.radar_altitude, grid2.radar_altitude)
        # _check_dicts_similar(grid1.radar_time, grid2.radar_time)
        # _check_dicts_similar(grid1.radar_name, grid2.radar_name)
        #  assert grid1.nradar == grid2.nradar


def _check_dicts_similar(dic1, dic2):
    for k, v in dic1.items():
        print("Checking key:", k)
        if k == 'data':
            assert_almost_equal(v, dic2[k])
        else:
            assert dic2[k] == v


def test_grid_class():
    grid = pyart.testing.make_target_grid()

    nz = 2
    ny = 400
    nx = 320
    nradar = 1

    assert grid.nx == nx
    assert grid.ny == ny
    assert grid.nz == nz
    assert grid.nradar == nradar

    assert isinstance(grid.metadata, dict)
    assert isinstance(grid.fields, dict)
    assert isinstance(grid.fields['reflectivity'], dict)
    assert grid.fields['reflectivity']['data'].shape == (nz, ny, nx)

    assert isinstance(grid.time, dict)
    assert grid.time['data'].shape == (1, )

    assert isinstance(grid.origin_longitude, dict)
    assert grid.origin_longitude['data'].shape == (1, )

    assert isinstance(grid.origin_latitude, dict)
    assert grid.origin_latitude['data'].shape == (1, )

    assert isinstance(grid.origin_altitude, dict)
    assert grid.origin_altitude['data'].shape == (1, )

    assert isinstance(grid.x, dict)
    assert grid.x['data'].shape == (nx, )

    assert isinstance(grid.y, dict)
    assert grid.y['data'].shape == (ny, )

    assert isinstance(grid.z, dict)
    assert grid.z['data'].shape == (nz, )

    assert isinstance(grid.point_x, LazyLoadDict)
    assert grid.point_x['data'].shape == (nz, ny, nx)
    assert_almost_equal(
        grid.point_x['data'][0, 0, :], grid.x['data'][:])
    assert_almost_equal(
        grid.point_x['data'][1, 1, :], grid.x['data'][:])

    assert isinstance(grid.point_y, LazyLoadDict)
    assert grid.point_y['data'].shape == (nz, ny, nx)
    assert_almost_equal(
        grid.point_y['data'][0, :, 0], grid.y['data'][:])
    assert_almost_equal(
        grid.point_y['data'][1, :, 1], grid.y['data'][:])

    assert isinstance(grid.point_z, LazyLoadDict)
    assert grid.point_z['data'].shape == (nz, ny, nx)
    assert_almost_equal(
        grid.point_z['data'][:, 0, 0], grid.z['data'][:])
    assert_almost_equal(
        grid.point_z['data'][:, 1, 1], grid.z['data'][:])

    assert isinstance(grid.point_longitude, LazyLoadDict)
    assert grid.point_longitude['data'].shape == (nz, ny, nx)
    assert_almost_equal(grid.point_longitude['data'][0, 0, 159], -98.1, 1)

    grid.init_point_longitude_latitude()
    assert isinstance(grid.point_latitude, LazyLoadDict)
    assert grid.point_latitude['data'].shape == (nz, ny, nx)
    assert_almost_equal(grid.point_latitude['data'][0, 200, 0], 36.7, 1)

    assert isinstance(grid.point_altitude, LazyLoadDict)
    assert grid.point_altitude['data'].shape == (nz, ny, nx)
    assert_almost_equal(grid.point_altitude['data'][0, 0, 0], 300, 0)
    assert_almost_equal(grid.point_altitude['data'][1, 0, 0], 800, 0)

    assert isinstance(grid.radar_latitude, dict)
    assert grid.radar_latitude['data'].shape == (nradar, )

    assert isinstance(grid.radar_longitude, dict)
    assert grid.radar_longitude['data'].shape == (nradar, )

    assert isinstance(grid.radar_altitude, dict)
    assert grid.radar_altitude['data'].shape == (nradar, )

    assert isinstance(grid.radar_time, dict)
    assert grid.radar_time['data'].shape == (nradar, )

    assert isinstance(grid.radar_name, dict)
    assert grid.radar_name['data'].shape == (nradar, )


def test_add_field():
    grid = pyart.testing.make_target_grid()

    assert 'velocity' not in grid.fields
    field_dict = {'data': np.ones((2, 400, 320))}
    grid.add_field('velocity', field_dict)
    assert 'velocity' in grid.fields
    assert grid.fields['velocity']['data'].shape == (2, 400, 320)

    # Existing fields can be replaced if requested
    assert_almost_equal(grid.fields['reflectivity']['data'][0, 0, 0], 0)
    grid.add_field('reflectivity', field_dict, replace_existing=True)
    assert_almost_equal(grid.fields['reflectivity']['data'][0, 0, 0], 1)


def test_add_field_raises():
    grid = pyart.testing.make_target_grid()

    # No 'data' key in field_dict raises a KeyError
    assert_raises(KeyError, grid.add_field, 'field_name', {})

    # Adding a field which already exists raises a ValueError
    field_dict = {'data': np.ones((2, 400, 320))}
    assert_raises(ValueError, grid.add_field, 'reflectivity', field_dict)

    # Incorrect field_shapes raise a ValueError
    field_dict = {'data': np.ones((1, 1, 1))}
    assert_raises(ValueError, grid.add_field, 'field_name', field_dict)


def test_projection_argument():
    grid = pyart.core.Grid(
        time={}, fields={}, metadata={},
        origin_latitude={}, origin_longitude={}, origin_altitude={},
        x={'data': [1]}, y={'data': [1]}, z={'data': [1]},
        radar_latitude={'data': [1]}, projection={})
    assert grid.projection == {}


def test_inconsistent_radar_arguments():

    #  partially instantiate a Grid class with fake data and a radar_latitude
    # argument indicating a single radar
    partial_grid = functools.partial(
        pyart.core.Grid, time={}, fields={}, metadata={},
        origin_latitude={}, origin_longitude={}, origin_altitude={},
        x={'data': [1]}, y={'data': [1]}, z={'data': [1]},
        radar_latitude={'data': [1]})

    # radar_ arguments indicating 2 radars should raise a ValueError
    bad = {'data': [1, 2]}
    assert_raises(ValueError, partial_grid, radar_longitude=bad)
    assert_raises(ValueError, partial_grid, radar_altitude=bad)
    assert_raises(ValueError, partial_grid, radar_time=bad)
    assert_raises(ValueError, partial_grid, radar_name=bad)


def test_init_point_altitude():
    grid = pyart.testing.make_target_grid()
    assert_almost_equal(grid.point_altitude['data'][0, 0, 0], 300, 0)
    assert_almost_equal(grid.point_altitude['data'][1, 0, 0], 800, 0)

    grid.origin_altitude['data'][0] = 555.
    assert_almost_equal(grid.point_altitude['data'][0, 0, 0], 300, 0)
    assert_almost_equal(grid.point_altitude['data'][1, 0, 0], 800, 0)

    grid.init_point_altitude()
    assert_almost_equal(grid.point_altitude['data'][0, 0, 0], 555, 0)
    assert_almost_equal(grid.point_altitude['data'][1, 0, 0], 1055, 0)


def test_init_point_longitude_latitude():
    grid = pyart.testing.make_target_grid()

    assert_almost_equal(grid.point_longitude['data'][0, 0, 159], -98.1, 1)
    assert_almost_equal(grid.point_latitude['data'][0, 200, 0], 36.7, 1)

    grid.origin_longitude['data'][0] = -80.0
    grid.origin_latitude['data'][0] = 40.0
    assert_almost_equal(grid.point_longitude['data'][0, 0, 159], -98.1, 1)
    assert_almost_equal(grid.point_latitude['data'][0, 200, 0], 36.7, 1)

    grid.init_point_longitude_latitude()
    assert_almost_equal(grid.point_longitude['data'][0, 0, 159], -80.0, 1)
    assert_almost_equal(grid.point_latitude['data'][0, 200, 0], 40.0, 1)


@skipif(not _PYPROJ_AVAILABLE)
def test_projection_proj():
    grid = pyart.testing.make_target_grid()
    grid.projection['proj'] = 'aeqd'
    assert isinstance(grid.projection_proj, pyproj.Proj)


def test_projection_proj_raised():
    grid = pyart.testing.make_target_grid()

    def access_projection_proj():
        grid.projection_proj

    assert_raises(ValueError, access_projection_proj)


def test_init_point_x_y_z():
    grid = pyart.testing.make_target_grid()
    assert isinstance(grid.point_x, LazyLoadDict)
    assert grid.point_x['data'].shape == (2, 400, 320)
    assert grid.point_x['data'][0, 0, 0] == grid.x['data'][0]

    grid.x['data'][0] = -500000.
    assert grid.point_x['data'][0, 0, 0] != grid.x['data'][0]

    grid.init_point_x_y_z()
    assert grid.point_x['data'][0, 0, 0] == grid.x['data'][0]


def test_get_point_longitude_latitude():
    grid = pyart.testing.make_target_grid()

    longitude, latitude = grid.get_point_longitude_latitude()
    assert latitude.shape == (400, 320)
    assert longitude.shape == (400, 320)
    assert_almost_equal(latitude[200, 160], 36.75, 2)
    assert_almost_equal(longitude[200, 160], -98.09, 2)

    longitude, latitude = grid.get_point_longitude_latitude(edges=True)
    assert latitude.shape == (401, 321)
    assert longitude.shape == (401, 321)
    assert_almost_equal(latitude[200, 160], 36.74, 2)
    assert_almost_equal(longitude[200, 160], -98.10, 2)


# Remove this test when Grid.axes is Deprecated
def test_grid_axes_attribute():
    # test the deprecated axes Grid attribute
    grid = pyart.testing.make_target_grid()

    nz = 2
    ny = 400
    nx = 320

    axes = grid.axes

    assert isinstance(axes, dict)

    assert isinstance(axes['time'], dict)
    assert axes['time']['data'].shape == (1, )

    assert isinstance(axes['time_start'], dict)
    assert axes['time_start']['data'].shape == (1, )

    assert isinstance(axes['time_end'], dict)
    assert axes['time_end']['data'].shape == (1, )

    assert isinstance(axes['lat'], dict)
    assert axes['lat']['data'].shape == (1, )

    assert isinstance(axes['lon'], dict)
    assert axes['lon']['data'].shape == (1, )

    assert isinstance(axes['alt'], dict)
    assert axes['alt']['data'].shape == (1, )

    assert isinstance(axes['x_disp'], dict)
    assert axes['x_disp']['data'].shape == (nx, )

    assert isinstance(axes['y_disp'], dict)
    assert axes['y_disp']['data'].shape == (ny, )

    assert isinstance(axes['z_disp'], dict)
    assert axes['z_disp']['data'].shape == (nz, )

    assert 'latitude' not in axes
    assert 'longitude' not in axes

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=DeprecationWarning)
        pyart.io.add_2d_latlon_axis(grid)

    assert isinstance(axes['latitude'], dict)
    assert axes['latitude']['data'].shape == (ny, nx)

    assert isinstance(axes['longitude'], dict)
    assert axes['longitude']['data'].shape == (ny, nx)


# Remove this function when Grid.from_legacy_parameters is Deprecated
def make_empty_grid(grid_shape, grid_limits):
    """
    Make an empty grid object without any fields or metadata.

    Parameters
    ----------
    grid_shape : 3-tuple of floats
        Number of points in the grid (z, y, x).
    grid_limits : 3-tuple of 2-tuples
        Minimum and maximum grid location (inclusive) in meters for the
        z, y, x coordinates.

    Returns
    -------
    grid : Grid
        Empty Grid object, centered near the ARM SGP site (Oklahoma).

    """
    time = {
        'data': np.array([0.0]),
        'units': 'seconds since 2000-01-01T00:00:00Z',
        'calendar': 'gregorian',
        'standard_name': 'time',
        'long_name': 'Time in seconds since volume start'}

    time_start = {
        'data': np.array([0.0]),
        'units': 'seconds since 2000-01-01T00:00:00Z',
        'calendar': 'gregorian',
        'standard_name': 'time',
        'long_name': 'Time in seconds since volume start'}

    time_end = {
        'data': np.array([0.0]),
        'units': 'seconds since 2000-01-01T00:00:00Z',
        'calendar': 'gregorian',
        'standard_name': 'time',
        'long_name': 'Time in seconds since volume start'}

    # grid coordinate dictionaries
    nz, ny, nx = grid_shape
    (z0, z1), (y0, y1), (x0, x1) = grid_limits

    xaxis = {'data': np.linspace(x0, x1, nx),
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

    altorigin = {'data': np.array([300.]),
                 'long_name': 'Altitude at grid origin',
                 'units': 'm',
                 'standard_name': 'altitude',
                 }

    latorigin = {'data': np.array([36.74]),
                 'long_name': 'Latitude at grid origin',
                 'units': 'degree_N',
                 'standard_name': 'latitude',
                 'valid_min': -90.,
                 'valid_max': 90.
                 }

    lonorigin = {'data': np.array([-98.1]),
                 'long_name': 'Longitude at grid origin',
                 'units': 'degree_E',
                 'standard_name': 'longitude',
                 'valid_min': -180.,
                 'valid_max': 180.
                 }

    axes = {'time': time,
            'time_start': time_start,
            'time_end': time_end,
            'z_disp': zaxis,
            'y_disp': yaxis,
            'x_disp': xaxis,
            'alt': altorigin,
            'lat': latorigin,
            'lon': lonorigin}

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=DeprecationWarning)
        return pyart.core.Grid.from_legacy_parameters({}, axes, {})


# Remove this test when Grid.from_legacy_parameters is Deprecated
def test_grid_from_legacy_parameters():
    grid_shape = (2, 3, 4)
    grid_limits = ((0, 500), (-400000, 400000), (-300000, 300000))
    grid = make_empty_grid(grid_shape, grid_limits)
    assert grid.nx == 4
    assert grid.ny == 3
    assert grid.nz == 2
