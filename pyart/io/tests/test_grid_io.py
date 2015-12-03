""" Unit Tests for Py-ART's io/grid_io.py module. """

from __future__ import print_function

import netCDF4
from numpy.testing import assert_almost_equal, assert_warns

import pyart


def test_grid_write_read():
    # test the read_grid and write_grid function by performing a
    # write/read roundtrip and comparing the two Grid objects
    grid1 = pyart.testing.make_target_grid()
    grid1.projection['comment'] = 'This is a comment'
    grid1.metadata['comment'] = 'This is another comment'

    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_grid.nc'
        pyart.io.write_grid(tmpfile, grid1)
        grid2 = pyart.io.read_grid(tmpfile)

        # check fields
        for field in grid1.fields.keys():
            _check_dicts_similar(grid1.fields[field], grid2.fields[field])

        # check attributes
        _check_attrs_similar(grid1, grid2, 'metadata')

        _check_attrs_similar(grid1, grid2, 'time')

        _check_attrs_similar(grid1, grid2, 'origin_latitude')
        _check_attrs_similar(grid1, grid2, 'origin_longitude')
        _check_attrs_similar(grid1, grid2, 'origin_altitude')

        _check_attrs_similar(grid1, grid2, 'regular_x')
        _check_attrs_similar(grid1, grid2, 'regular_y')
        _check_attrs_similar(grid1, grid2, 'regular_z')

        _check_attrs_similar(grid1, grid2, 'point_x')
        _check_attrs_similar(grid1, grid2, 'point_y')
        _check_attrs_similar(grid1, grid2, 'point_z')

        _check_attrs_similar(grid1, grid2, 'projection')
        assert grid1.projection['_include_lon_0_lat_0'] is True

        _check_attrs_similar(grid1, grid2, 'point_latitude')
        _check_attrs_similar(grid1, grid2, 'point_longitude')
        _check_attrs_similar(grid1, grid2, 'point_altitude')

        assert grid1.nx == grid2.nx
        assert grid1.ny == grid2.ny
        assert grid1.nz == grid2.nz

        _check_attrs_similar(grid1, grid2, 'radar_latitude')
        _check_attrs_similar(grid1, grid2, 'radar_longitude')
        _check_attrs_similar(grid1, grid2, 'radar_altitude')
        _check_attrs_similar(grid1, grid2, 'radar_time')
        _check_attrs_similar(grid1, grid2, 'radar_name')
        assert grid1.nradar == grid2.nradar


def _check_attrs_similar(grid1, grid2, attr):
    print("Checking attribute:", attr)
    dic1 = getattr(grid1, attr)
    dic2 = getattr(grid2, attr)
    _check_dicts_similar(dic1, dic2)


def _check_dicts_similar(dic1, dic2):
    for k, v in dic1.items():
        print("Checking key:", k)
        if k == 'data':
            if v.dtype.char == 'S':
                assert ''.join(*v) == ''.join(*dic2[k])
            else:
                assert_almost_equal(v, dic2[k])
        else:
            assert dic2[k] == v


def test_grid_write_point_vars():
    grid1 = pyart.testing.make_target_grid()

    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_grid.nc'
        pyart.io.write_grid(tmpfile, grid1, write_point_x_y_z=True,
                            write_point_lon_lat_alt=True)
        dset = netCDF4.Dataset(tmpfile, 'r')

        assert 'point_x' in dset.variables
        assert 'point_y' in dset.variables
        assert 'point_z' in dset.variables
        assert 'point_latitude' in dset.variables
        assert 'point_longitude' in dset.variables
        assert 'point_altitude' in dset.variables
        dset.close()


def test_grid_write_arm_time_vars():
    grid1 = pyart.testing.make_target_grid()

    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_grid.nc'
        pyart.io.write_grid(tmpfile, grid1, arm_time_variables=True)
        dset = netCDF4.Dataset(tmpfile, 'r')
        assert 'base_time' in dset.variables
        assert 'time_offset' in dset.variables
        dset.close()


def test_exclude_fields_argument():
    grid1 = pyart.testing.make_target_grid()
    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_grid.nc'
        pyart.io.write_grid(tmpfile, grid1)
        grid2 = pyart.io.read_grid(tmpfile, exclude_fields=['reflectivity'])
        assert 'reflectivity' in grid1.fields
        assert 'reflectivity' not in grid2.fields


def test_radar_attrs_none():
    grid1 = pyart.testing.make_target_grid()
    grid1.radar_latitude = None
    grid1.radar_longitude = None
    grid1.radar_altitude = None
    grid1.radar_time = None
    grid1.radar_name = None
    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_grid.nc'
        pyart.io.write_grid(tmpfile, grid1)
        grid2 = pyart.io.read_grid(tmpfile)
        assert grid2.radar_latitude is None
        assert grid2.radar_longitude is None
        assert grid2.radar_altitude is None
        assert grid2.radar_time is None
        assert grid2.radar_name is None


def test_bad_shaped_field():
    grid1 = pyart.testing.make_target_grid()
    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_grid.nc'
        pyart.io.write_grid(tmpfile, grid1)

        # add a scalar variable named bad_field
        dset = netCDF4.Dataset(tmpfile, 'a')
        dset.createVariable('bad_field', 'f4', ())
        dset.close()

        # warning should be raised about incorrect shape of bad_field var
        assert_warns(UserWarning, pyart.io.read_grid, tmpfile)
