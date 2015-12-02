""" Unit Tests for Py-ART's core/grid.py module. """

from __future__ import print_function

import numpy as np

import pyart

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
