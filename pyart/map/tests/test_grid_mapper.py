""" Unit Tests for Py-ART's map/grid_mapper.py module. """

from __future__ import print_function

import numpy as np
from numpy.testing import assert_almost_equal, assert_raises

import pyart

EXPECTED_CENTER_SLICE = [40, 30, 20, 10, 0, 0, 10, 20, 30, 40]

COMMON_MAP_TO_GRID_ARGS = {
    'grid_shape': (3, 9, 10),
    'grid_limits': ((-400.0, 400.0), (-900.0, 900.0), (-900, 900)),
    'fields': None,
    'roi_func': lambda z, y, x: 30, }


def test_map_to_grid_filter():

    # simulate a radar with bad gates which reports huge reflectivities
    radar = pyart.testing.make_target_radar()
    radar.fields['reflectivity']['data'][0:100, 25] = 99999.0

    # without filtering bad gates leaks through
    gatefilter = pyart.filters.GateFilter(radar)
    grids = pyart.map.map_to_grid(
        (radar,), gatefilters=(gatefilter, ), **COMMON_MAP_TO_GRID_ARGS)
    assert grids['reflectivity'].max() > 41.0

    # with filtering bad gates is supressed
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_above('reflectivity', 41.0)
    grids = pyart.map.map_to_grid(
        (radar,), gatefilters=(gatefilter, ), **COMMON_MAP_TO_GRID_ARGS)
    assert grids['reflectivity'].max() < 41.0


def test_map_to_grid_non_tuple():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_to_grid(radar, **COMMON_MAP_TO_GRID_ARGS)
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(center_slice, EXPECTED_CENTER_SLICE)


def test_map_to_grid_default():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_to_grid((radar,), **COMMON_MAP_TO_GRID_ARGS)
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(center_slice, EXPECTED_CENTER_SLICE)


def test_map_to_grid_cressman():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_to_grid(
        (radar,), (3, 9, 10), ((-400.0, 400.0), (-900.0, 900.0), (-900, 900)),
        roi_func='constant', constant_roi=30., weighting_function='CRESSMAN')
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(center_slice, EXPECTED_CENTER_SLICE)


def test_map_to_grid_constant_roi():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_to_grid(
        (radar,), (3, 9, 10), ((-400.0, 400.0), (-900.0, 900.0), (-900, 900)),
        roi_func='constant', constant_roi=30.)
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(center_slice, EXPECTED_CENTER_SLICE)


def test_map_to_grid_dist_roi():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_to_grid(
        (radar,), (3, 9, 10), ((-400.0, 400.0), (-900.0, 900.0), (-900, 900)),
        roi_func='dist', z_factor=0, xy_factor=0, min_radius=30.)
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(center_slice, EXPECTED_CENTER_SLICE)


def test_map_to_grid_dist_beam_roi():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_to_grid(
        (radar,),
        grid_shape=(3, 9, 10),
        grid_limits=((-400.0, 400.0), (-900.0, 900.0), (-900, 900)),
        fields=['reflectivity'],
        min_radius=30, bsp=0., h_factor=0.)
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(center_slice, EXPECTED_CENTER_SLICE)


def test_map_to_grid_default_two_radars():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_to_grid((radar, radar), **COMMON_MAP_TO_GRID_ARGS)
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(center_slice, EXPECTED_CENTER_SLICE)


def test_map_to_grid_masked_refl_field():
    radar = pyart.testing.make_target_radar()

    # mask the last gate of the first ray
    fdata = radar.fields['reflectivity']['data']
    fdata = np.ma.masked_invalid(fdata)
    fdata.mask[0, -1] = True
    radar.fields['reflectivity']['data'] = fdata

    grids = pyart.map.map_to_grid((radar,), **COMMON_MAP_TO_GRID_ARGS)
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(center_slice, EXPECTED_CENTER_SLICE)


def test_map_to_grid_barnes_filter_masked():

    radar = pyart.testing.make_target_radar()

    # mask and replace with non-sense values a region of radar volume
    fdata = radar.fields['reflectivity']['data']
    fdata = np.ma.array(fdata)
    fdata[100:120, 20:40] = -9999.
    fdata[100:120, 20:40] = np.ma.masked
    radar.fields['reflectivity']['data'] = fdata

    # check that non-sense values are masked
    assert radar.fields['reflectivity']['data'].min() > -1.

    grids = pyart.map.map_to_grid(
        (radar,), weighting_function='Barnes', **COMMON_MAP_TO_GRID_ARGS)

    # check that no masked values were included in the interpolation
    assert grids['reflectivity'].min() > -1.


def test_map_to_grid_no_copy():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_to_grid(
        (radar,), copy_field_data=False, **COMMON_MAP_TO_GRID_ARGS)
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(center_slice, EXPECTED_CENTER_SLICE)


def test_map_to_grid_no_copy_two_radars():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_to_grid(
        (radar, radar), copy_field_data=False, **COMMON_MAP_TO_GRID_ARGS)
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(center_slice, EXPECTED_CENTER_SLICE)


def test_map_to_grid_tiny_grid():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_to_grid(
        (radar,),
        grid_shape=(1, 1, 1),
        grid_limits=((-400.0, 400.0), (-900.0, 900.0), (-900, 900)),
        fields=['reflectivity'])
    assert grids['reflectivity'].shape == (1, 1, 1)
    assert int(grids['reflectivity'][0]) == 40


def test_map_to_grid_errors():
    assert_raises(ValueError, pyart.map.map_to_grid, (None, ), (1, 1, 1),
                  ((-1, 1), (-1, 1), (-1, 1)),
                  weighting_function='foo')
    assert_raises(ValueError, pyart.map.map_to_grid, (None, ), (1, 1, 1),
                  ((-1, 1), (-1, 1), (-1, 1)),
                  algorithm='foo')
    radar = pyart.testing.make_target_radar()
    assert_raises(ValueError, pyart.map.map_to_grid, (radar, ), (1, 1, 1),
                  ((-1, 1), (-1, 1), (-1, 1)),
                  roi_func='foo')


def test_grid_from_radars_errors():
    assert_raises(ValueError, pyart.map.grid_from_radars, (None, ), (1, 1, 1),
                  ((-1, 1), (-1, 1), (-1, 1)), gridding_algo='foo')


def test_grid_from_radars():
    radar = pyart.testing.make_target_radar()
    grid = pyart.map.grid_from_radars((radar,), **COMMON_MAP_TO_GRID_ARGS)

    # check field data
    center_slice = grid.fields['reflectivity']['data'][1, 4, :]
    assert_almost_equal(center_slice, EXPECTED_CENTER_SLICE)

    # check other Grid object attributes
    assert 'ROI' in grid.fields
    assert np.all(grid.fields['ROI']['data'] == 30.)
    assert_almost_equal(grid.x['data'], np.linspace(-900, 900, 10))
    assert_almost_equal(grid.y['data'], np.linspace(-900, 900, 9))
    assert_almost_equal(grid.z['data'], np.linspace(-400, 400, 3))

    # check that grid.radar_ attributes set correctly
    assert isinstance(grid.radar_latitude, dict)
    assert_almost_equal(grid.radar_latitude['data'][0], 36.5)

    assert isinstance(grid.radar_longitude, dict)
    assert_almost_equal(grid.radar_longitude['data'][0], -97.5)

    assert isinstance(grid.radar_altitude, dict)
    assert_almost_equal(grid.radar_altitude['data'][0], 200)

    assert isinstance(grid.radar_time, dict)
    assert_almost_equal(grid.radar_time['data'][0], 0.0)
    assert grid.radar_time['units'] == radar.time['units']

    assert isinstance(grid.radar_name, dict)
    assert grid.radar_name['data'][0] == 'fake_radar'

    assert grid.nradar == 1


def test_unify_times_for_radars():
    radar1 = pyart.testing.make_target_radar()
    radar2 = pyart.testing.make_target_radar()
    radar2.time['units'] = 'seconds since 1989-01-01T00:00:02Z'
    radars = (radar1, radar2)
    times, units = pyart.map.grid_mapper._unify_times_for_radars(radars)
    assert_almost_equal(times[0], 0)
    assert_almost_equal(times[1], 1)
    assert units == 'seconds since 1989-01-01T00:00:01Z'


def test_grid_from_radars_non_tuple():
    radar = pyart.testing.make_target_radar()
    grid = pyart.map.grid_from_radars(radar, **COMMON_MAP_TO_GRID_ARGS)

    # check field data
    center_slice = grid.fields['reflectivity']['data'][1, 4, :]
    assert_almost_equal(center_slice, EXPECTED_CENTER_SLICE)

    # check other Grid object attributes
    assert 'ROI' in grid.fields
    assert np.all(grid.fields['ROI']['data'] == 30.)
    assert_almost_equal(grid.x['data'], np.linspace(-900, 900, 10))
    assert_almost_equal(grid.y['data'], np.linspace(-900, 900, 9))
    assert_almost_equal(grid.z['data'], np.linspace(-400, 400, 3))


def test_grid_from_radars_grid_origin():
    radar = pyart.testing.make_target_radar()
    radar.metadata.pop('instrument_name')
    grid = pyart.map.grid_from_radars((radar,), grid_origin=(36.4, -97.6),
                                      grid_origin_alt=200,
                                      **COMMON_MAP_TO_GRID_ARGS)
    assert_almost_equal(grid.origin_latitude['data'][0], 36.4, 1)
    assert_almost_equal(grid.origin_longitude['data'], -97.6, 1)


def test_example_roi_funcs():
    assert pyart.map.example_roi_func_constant(0, 0, 0) == 500.
    assert pyart.map.example_roi_func_dist(0, 0, 0) == 500.
    assert pyart.map.example_roi_func_dist_beam(0, 0, 0) == 500.
