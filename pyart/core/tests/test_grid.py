""" Unit Tests for Py-ART's core/grid.py module. """

from __future__ import print_function

import functools

import numpy as np
from numpy.testing import assert_almost_equal, assert_raises

import pyart
from pyart.lazydict import LazyLoadDict

COMMON_MAP_TO_GRID_ARGS = {
    'grid_shape': (3, 9, 10),
    'grid_limits': ((-400.0, 400.0), (-900.0, 900.0), (-900, 900)),
    'fields': ['reflectivity'],
    'roi_func': lambda z, y, x: 30, }


def test_grid_from_radars():
    radar = pyart.testing.make_target_radar()
    grid = pyart.map.grid_from_radars((radar,), **COMMON_MAP_TO_GRID_ARGS)

    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_grid.nc'
        grid.write(tmpfile)
        grid2 = pyart.io.read_grid(tmpfile)

        # check metadata
        for k, v in grid.metadata.items():
            print("Checking key:", k, "should have value:", v)
            print(grid2.metadata)
            assert grid2.metadata[k] == v

        # check axes
        for axes_key in grid.axes.keys():
            for k, v in grid.axes[axes_key].items():
                print("Checking axes_key:", axes_key, "key:", k)
                if k == 'data':
                    assert np.all(grid.axes[axes_key][k] == v)
                else:
                    assert grid2.axes[axes_key][k] == v

        # check fields
        for field in grid.fields.keys():
            for k, v in grid.fields[field].items():
                print("Checking field:", field, "key:", k)
                if k == 'data':
                    assert np.all(grid.fields[field][k] == v)
                else:
                    assert grid2.fields[field][k] == v


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

    assert isinstance(grid.regular_x, dict)
    assert grid.regular_x['data'].shape == (nx, )

    assert isinstance(grid.regular_y, dict)
    assert grid.regular_y['data'].shape == (ny, )

    assert isinstance(grid.regular_z, dict)
    assert grid.regular_z['data'].shape == (nz, )

    assert isinstance(grid.point_x, LazyLoadDict)
    assert grid.point_x['data'].shape == (nz, ny, nx)
    assert_almost_equal(
        grid.point_x['data'][0, 0, :], grid.regular_x['data'][:])
    assert_almost_equal(
        grid.point_x['data'][1, 1, :], grid.regular_x['data'][:])

    assert isinstance(grid.point_y, LazyLoadDict)
    assert grid.point_y['data'].shape == (nz, ny, nx)
    assert_almost_equal(
        grid.point_y['data'][0, :, 0], grid.regular_y['data'][:])
    assert_almost_equal(
        grid.point_y['data'][1, :, 1], grid.regular_y['data'][:])

    assert isinstance(grid.point_z, LazyLoadDict)
    assert grid.point_z['data'].shape == (nz, ny, nx)
    assert_almost_equal(
        grid.point_z['data'][:, 0, 0], grid.regular_z['data'][:])
    assert_almost_equal(
        grid.point_z['data'][:, 1, 1], grid.regular_z['data'][:])

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


def test_inconsistent_radar_arguments():

    #  partially instantiate a Grid class with fake data and a radar_latitude
    # argument indicating a single radar
    partial_grid = functools.partial(
        pyart.core.Grid, time={}, fields={}, metadata={},
        origin_latitude={}, origin_longitude={}, origin_altitude={},
        regular_x={'data': [1]}, regular_y={'data': [1]},
        regular_z={'data': [1]}, radar_latitude={'data': [1]})

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


def test_init_point_x_y_z():
    grid = pyart.testing.make_target_grid()
    assert isinstance(grid.point_x, LazyLoadDict)
    assert grid.point_x['data'].shape == (2, 400, 320)
    assert grid.point_x['data'][0, 0, 0] == grid.regular_x['data'][0]

    grid.regular_x['data'][0] = -500000.
    assert grid.point_x['data'][0, 0, 0] != grid.regular_x['data'][0]

    grid.init_point_x_y_z()
    assert grid.point_x['data'][0, 0, 0] == grid.regular_x['data'][0]


def test_grid_axes_attribute():
    # test the depreciated axes Grid attribute
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

    pyart.io.add_2d_latlon_axis(grid)

    assert isinstance(axes['latitude'], dict)
    assert axes['latitude']['data'].shape == (ny, nx)

    assert isinstance(axes['longitude'], dict)
    assert axes['longitude']['data'].shape == (ny, nx)
