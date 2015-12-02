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
