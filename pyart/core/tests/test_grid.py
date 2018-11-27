""" Unit Tests for Py-ART's core/grid.py module. """

from __future__ import print_function

import functools
import pickle

import numpy as np
from numpy.testing import assert_almost_equal
import pytest

try:
    import pyproj
    _PYPROJ_AVAILABLE = True
except ImportError:
    _PYPROJ_AVAILABLE = False

import pyart
from pyart.lazydict import LazyLoadDict


def test_grid_picklable():
    grid = pyart.testing.make_target_grid()
    picklestring = pickle.dumps(grid)
    grid_new = pickle.loads(picklestring)
    assert 'data' in grid.point_x
    assert 'data' in grid_new.point_x


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
    pytest.raises(KeyError, grid.add_field, 'field_name', {})

    # Adding a field which already exists raises a ValueError
    field_dict = {'data': np.ones((2, 400, 320))}
    pytest.raises(ValueError, grid.add_field, 'reflectivity', field_dict)

    # Incorrect field_shapes raise a ValueError
    field_dict = {'data': np.ones((1, 1, 1))}
    pytest.raises(ValueError, grid.add_field, 'field_name', field_dict)


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
    pytest.raises(ValueError, partial_grid, radar_longitude=bad)
    pytest.raises(ValueError, partial_grid, radar_altitude=bad)
    pytest.raises(ValueError, partial_grid, radar_time=bad)
    pytest.raises(ValueError, partial_grid, radar_name=bad)


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


@pytest.mark.skipif(not _PYPROJ_AVAILABLE,
                    reason="PyProj is not installed.")
def test_projection_proj():
    grid = pyart.testing.make_target_grid()
    grid.projection['proj'] = 'aeqd'
    assert isinstance(grid.projection_proj, pyproj.Proj)


def test_projection_proj_raised():
    grid = pyart.testing.make_target_grid()

    def access_projection_proj():
        grid.projection_proj

    pytest.raises(ValueError, access_projection_proj)


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
