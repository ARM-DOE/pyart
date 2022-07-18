""" Unit Tests for Py-ART's io/grid_io.py module. """

import netCDF4
import numpy as np
from numpy.testing import assert_almost_equal, assert_warns

import pyart
from pyart.io.common import stringarray_to_chararray


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
        assert 'Conventions' in grid2.metadata
        grid2.metadata.pop('Conventions')
        _check_attrs_similar(grid1, grid2, 'metadata')

        _check_attrs_similar(grid1, grid2, 'time')

        _check_attrs_similar(grid1, grid2, 'origin_latitude')
        _check_attrs_similar(grid1, grid2, 'origin_longitude')
        _check_attrs_similar(grid1, grid2, 'origin_altitude')

        _check_attrs_similar(grid1, grid2, 'x')
        _check_attrs_similar(grid1, grid2, 'y')
        _check_attrs_similar(grid1, grid2, 'z')

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
            if v.dtype.char == 'S' or v.dtype.char == 'U':
                s1 = v.astype('S')
                if s1.ndim == 1:
                    s1 = stringarray_to_chararray(s1)

                s2 = dic2[k].astype('S')
                if s1.ndim == 1:
                    s2 = stringarray_to_chararray(s2)

                assert b''.join(*s1) == b''.join(*s2)
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


def test_write_projection_coordinate_system():
    grid_shape = (2, 3, 4)
    grid_limits = ((0, 500), (-400000, 400000), (-300000, 300000))
    grid = pyart.testing.make_empty_grid(grid_shape, grid_limits)

    with pyart.testing.InTemporaryDirectory():

        # standard projection
        grid.projection['proj'] = 'pyart_aeqd'
        tmpfile = 'tmp_grid1.nc'
        pyart.io.write_grid(tmpfile, grid, write_proj_coord_sys=True)
        dset = netCDF4.Dataset(tmpfile, 'r')
        assert 'ProjectionCoordinateSystem' in dset.variables
        dset.close()

        # explicit proj_coord_sys
        proj_coord_sys = {'foo': 'bar'}
        tmpfile = 'tmp_grid2.nc'
        pyart.io.write_grid(tmpfile, grid, write_proj_coord_sys=True,
                            proj_coord_sys=proj_coord_sys)
        dset = netCDF4.Dataset(tmpfile, 'r')
        assert 'ProjectionCoordinateSystem' in dset.variables
        assert dset.variables['ProjectionCoordinateSystem'].foo == 'bar'
        dset.close()

        # explicit proj_coord_sys
        proj_coord_sys = {'foo': 'bar'}
        tmpfile = 'tmp_grid3.nc'
        pyart.io.write_grid(tmpfile, grid, write_proj_coord_sys=False)
        dset = netCDF4.Dataset(tmpfile, 'r')
        assert 'ProjectionCoordinateSystem' not in dset.variables
        dset.close()


def test_unknown_projection():
    grid_shape = (2, 3, 4)
    grid_limits = ((0, 500), (-400000, 400000), (-300000, 300000))
    grid = pyart.testing.make_empty_grid(grid_shape, grid_limits)
    grid.projection['proj'] = 'null'

    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_grid.nc'
        assert_warns(UserWarning, pyart.io.write_grid, tmpfile, grid,
                     write_proj_coord_sys=True)
        dset = netCDF4.Dataset(tmpfile, 'r')
        assert 'ProjectionCoordinateSystem' not in dset.variables
        dset.close()


def test_make_coordinate_system_dict():
    grid_shape = (2, 3, 4)
    grid_limits = ((0, 500), (-400000, 400000), (-300000, 300000))
    grid = pyart.testing.make_empty_grid(grid_shape, grid_limits)

    # pyart_aeqd
    grid.projection['proj'] = 'pyart_aeqd'
    dic = pyart.io.grid_io._make_coordinatesystem_dict(grid)
    assert dic['grid_mapping_name'] == 'azimuthal_equidistant'
    assert 'semi_major_axis' in dic
    assert 'inverse_flattening' in dic
    assert 'longitude_of_prime_meridian' in dic
    assert 'false_easting' in dic
    assert 'false_northing' in dic
    assert 'longitude_of_projection_origin' in dic
    assert 'latitude_of_projection_origin' in dic

    # aeqd
    grid.projection['proj'] = 'aeqd'
    dic = pyart.io.grid_io._make_coordinatesystem_dict(grid)
    assert dic['grid_mapping_name'] == 'azimuthal_equidistant'
    assert 'semi_major_axis' in dic
    assert 'inverse_flattening' in dic
    assert 'longitude_of_prime_meridian' in dic
    assert 'false_easting' in dic
    assert 'false_northing' in dic
    assert 'longitude_of_projection_origin' in dic
    assert 'latitude_of_projection_origin' in dic

    # tmerc
    grid.projection['proj'] = 'tmerc'
    dic = pyart.io.grid_io._make_coordinatesystem_dict(grid)
    assert dic['grid_mapping_name'] == 'transverse_mercator'
    assert 'longitude_of_central_meridian' in dic
    assert 'latitude_of_projection_origin' in dic
    assert 'scale_factor_at_central_meridian' in dic

    # lcc
    grid.projection['proj'] = 'lcc'
    dic = pyart.io.grid_io._make_coordinatesystem_dict(grid)
    assert dic['grid_mapping_name'] == 'lambert_conformal_conic'
    assert 'standard_parallel' in dic
    assert 'longitude_of_projection_origin' in dic
    assert 'latitude_of_projection_origin' in dic

    # laea
    grid.projection['proj'] = 'laea'
    dic = pyart.io.grid_io._make_coordinatesystem_dict(grid)
    assert dic['grid_mapping_name'] == 'lambert_azimuthal_equal_area'
    assert 'longitude_of_projection_origin' in dic
    assert 'latitude_of_projection_origin' in dic

    # aea
    grid.projection['proj'] = 'aea'
    dic = pyart.io.grid_io._make_coordinatesystem_dict(grid)
    assert dic['grid_mapping_name'] == 'albers_conical_equal_area'
    assert 'standard_parallel' in dic
    assert 'longitude_of_central_meridian' in dic
    assert 'latitude_of_projection_origin' in dic

    # stere
    grid.projection['proj'] = 'stere'
    dic = pyart.io.grid_io._make_coordinatesystem_dict(grid)
    assert dic['grid_mapping_name'] == 'stereographic'
    assert 'longitude_of_projection_origin' in dic
    assert 'latitude_of_projection_origin' in dic
    assert 'scale_factor_at_projection_origin' in dic

    # npstere
    grid.projection['proj'] = 'npstere'
    dic = pyart.io.grid_io._make_coordinatesystem_dict(grid)
    assert dic['grid_mapping_name'] == 'polar_stereographic'
    assert 'longitude_of_projection_origin' in dic
    assert 'latitude_of_projection_origin' in dic
    assert 'standard_parallel' in dic

    # spstere
    grid.projection['proj'] = 'npstere'
    dic = pyart.io.grid_io._make_coordinatesystem_dict(grid)
    assert dic['grid_mapping_name'] == 'polar_stereographic'
    assert 'longitude_of_projection_origin' in dic
    assert 'latitude_of_projection_origin' in dic
    assert 'standard_parallel' in dic

    # ortho
    grid.projection['proj'] = 'ortho'
    dic = pyart.io.grid_io._make_coordinatesystem_dict(grid)
    assert dic['grid_mapping_name'] == 'orthographic'
    assert 'longitude_of_projection_origin' in dic
    assert 'latitude_of_projection_origin' in dic

    # projection that does not map to a CDM transform
    grid.projection['proj'] = 'null'
    dic = pyart.io.grid_io._make_coordinatesystem_dict(grid)
    assert dic is None


def test_write_grid_empty_radar_names():
    # GitHub issue #537
    grid1 = pyart.testing.make_target_grid()
    grid1.nradar = 2
    grid1.radar_name['data'] = np.array(['', ''])
    grid1.radar_latitude = None
    grid1.radar_longitude = None
    grid1.radar_altitude = None
    grid1.radar_time = None

    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_grid.nc'
        pyart.io.write_grid(tmpfile, grid1)
