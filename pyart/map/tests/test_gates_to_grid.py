""" Unit Tests for Py-ART's map/gates_to_grid.py. """

from __future__ import print_function

import numpy as np
from numpy.testing import assert_almost_equal, assert_raises

import pyart

EXPECTED_CENTER_SLICE = [40, 30, 20, 10, 0, 0, 10, 20, 30, 40]

COMMON_MAP_TO_GRID_ARGS = {
    'grid_shape': (3, 9, 10),
    'grid_limits': ((-400.0, 400.0), (-900.0, 900.0), (-900, 900)),
    'fields': None,
    'roi_func': 'constant',
    'constant_roi': 30., }


def test_map_to_grid_filter():

    # simulate a radar with bad gates which reports huge reflectivities
    radar = pyart.testing.make_target_radar()
    radar.fields['reflectivity']['data'][0:100, 25] = 99999.0

    # without filtering bad gates leaks through
    gatefilter = pyart.filters.GateFilter(radar)
    grids = pyart.map.map_gates_to_grid(
        (radar,), gatefilters=(gatefilter, ), **COMMON_MAP_TO_GRID_ARGS)
    assert grids['reflectivity'].max() > 41.0

    # with filtering bad gates is supressed
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_above('reflectivity', 41.0)
    grids = pyart.map.map_gates_to_grid(
        (radar,), gatefilters=(gatefilter, ), **COMMON_MAP_TO_GRID_ARGS)
    assert grids['reflectivity'].max() < 41.0


def test_map_to_grid_non_tuple():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_gates_to_grid(radar,
                                        **COMMON_MAP_TO_GRID_ARGS)
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(np.round(center_slice), EXPECTED_CENTER_SLICE)


def test_map_to_grid_default():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_gates_to_grid((radar,),
                                        **COMMON_MAP_TO_GRID_ARGS)
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(np.round(center_slice), EXPECTED_CENTER_SLICE)


def test_map_to_grid_cressman():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_gates_to_grid(
        (radar,), (3, 9, 10), ((-400.0, 400.0), (-900.0, 900.0), (-900, 900)),
        roi_func='constant', constant_roi=30., weighting_function='CRESSMAN')
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(np.round(center_slice), EXPECTED_CENTER_SLICE)


def test_map_to_grid_constant_roi():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_gates_to_grid(
        (radar,), (3, 9, 10), ((-400.0, 400.0), (-900.0, 900.0), (-900, 900)),
        roi_func='constant', constant_roi=30.)
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(np.round(center_slice), EXPECTED_CENTER_SLICE)


def test_map_to_grid_dist_roi():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_gates_to_grid(
        (radar,), (3, 9, 10), ((-400.0, 400.0), (-900.0, 900.0), (-900, 900)),
        roi_func='dist', z_factor=0, xy_factor=0, min_radius=30.)
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(np.round(center_slice), EXPECTED_CENTER_SLICE)


def test_map_to_grid_dist_beam_roi():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_gates_to_grid(
        (radar,),
        grid_shape=(3, 9, 10),
        grid_limits=((-400.0, 400.0), (-900.0, 900.0), (-900, 900)),
        fields=['reflectivity'],
        min_radius=30, bsp=0., h_factor=0.)
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(np.round(center_slice), EXPECTED_CENTER_SLICE)


def test_map_to_grid_default_two_radars():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_gates_to_grid((radar, radar),
                                        **COMMON_MAP_TO_GRID_ARGS)
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(np.round(center_slice), EXPECTED_CENTER_SLICE)


def test_map_to_grid_masked_refl_field():
    radar = pyart.testing.make_target_radar()

    # mask the last gate of the first ray
    fdata = radar.fields['reflectivity']['data']
    fdata = np.ma.masked_invalid(fdata)
    fdata.mask[0, -1] = True
    radar.fields['reflectivity']['data'] = fdata

    grids = pyart.map.map_gates_to_grid((radar,),
                                        **COMMON_MAP_TO_GRID_ARGS)
    center_slice = grids['reflectivity'][1, 4, :]
    assert_almost_equal(np.round(center_slice), EXPECTED_CENTER_SLICE)


def test_map_to_grid_tiny_grid():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_gates_to_grid(
        (radar,),
        grid_shape=(1, 1, 1),
        grid_limits=((-400.0, 400.0), (-900.0, 900.0), (-900, 900)),
        fields=['reflectivity'])
    assert grids['reflectivity'].shape == (1, 1, 1)
    assert abs(np.round(grids['reflectivity'][0]) - 40.0) < 0.01


def test_grid_from_radars_gates_to_grid():
    radar = pyart.testing.make_target_radar()
    grid = pyart.map.grid_from_radars(
        (radar,), gridding_algo='map_gates_to_grid', **COMMON_MAP_TO_GRID_ARGS)

    # check field data
    center_slice = grid.fields['reflectivity']['data'][1, 4, :]
    assert_almost_equal(np.round(center_slice), EXPECTED_CENTER_SLICE)

    # check other Grid object attributes
    assert 'ROI' in grid.fields
    assert np.all(grid.fields['ROI']['data'] == 30.)
    assert_almost_equal(grid.x['data'], np.linspace(-900, 900, 10))
    assert_almost_equal(grid.y['data'], np.linspace(-900, 900, 9))
    assert_almost_equal(grid.z['data'], np.linspace(-400, 400, 3))


def test_map_to_grid_errors():
    radar = pyart.testing.make_target_radar()

    # invalid weighting_function
    assert_raises(ValueError, pyart.map.map_gates_to_grid, (radar, ),
                  (1, 1, 1), ((-1, 1), (-1, 1), (-1, 1)),
                  weighting_function='foo')

    # invalid roi_func
    assert_raises(ValueError, pyart.map.map_gates_to_grid, (radar, ),
                  (1, 1, 1), ((-1, 1), (-1, 1), (-1, 1)),
                  roi_func='foo')

    # missing reflectivity field
    radar.fields.pop('reflectivity')
    assert_raises(ValueError, pyart.map.map_gates_to_grid, (radar, ),
                  (1, 1, 1), ((-1, 1), (-1, 1), (-1, 1)))


def test_grid_from_radars():
    radar = pyart.testing.make_target_radar()
    grid = pyart.map.grid_from_radars((radar,), **COMMON_MAP_TO_GRID_ARGS)

    # check field data
    center_slice = grid.fields['reflectivity']['data'][1, 4, :]
    assert_almost_equal(np.round(center_slice), EXPECTED_CENTER_SLICE)

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
                                      **COMMON_MAP_TO_GRID_ARGS)
    assert_almost_equal(grid.origin_latitude['data'][0], 36.4, 1)
    assert_almost_equal(grid.origin_longitude['data'][0], -97.6, 1)


def test_example_roi_funcs():
    assert pyart.map.example_roi_func_constant(0, 0, 0) == 500.
    assert pyart.map.example_roi_func_dist(0, 0, 0) == 500.
    assert pyart.map.example_roi_func_dist_beam(0, 0, 0) == 500.
