""" Unit Tests for Py-ART's core/grid.py module. """

from __future__ import print_function

import tempfile
import os

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

    tmpfile = tempfile.mkstemp(suffix='.nc', dir='.')[1]
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

    os.remove(tmpfile)
