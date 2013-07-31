""" Unit Tests for Py-ART's map/grid_mapper.py module. """

import numpy as np
from numpy.testing import assert_array_equal

import pyart

EXPECTED_CENTER_SLICE = [40, 30, 20, 10, 0, 0, 10, 20, 30, 40]

COMMON_MAP_TO_GRID_ARGS = {
    'grid_shape': (10, 9, 3),
    'grid_limits': ((-900.0, 900.0), (-900.0, 900.0), (-400, 400)),
    'fields': ['reflectivity_horizontal'],
    'qrf_func': lambda x, y, z: 30, }


def test_map_to_grid_default():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_to_grid((radar,),
                                  **COMMON_MAP_TO_GRID_ARGS)
    center_slice = grids['reflectivity_horizontal'][1, 4, :]
    assert_array_equal(np.round(center_slice), EXPECTED_CENTER_SLICE)


def test_map_to_grid_qrf():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_to_grid(
        (radar,),
        grid_shape=(10, 9, 3),
        grid_limits=((-900.0, 900.0), (-900.0, 900.0), (-400, 400)),
        fields=['reflectivity_horizontal'],
        min_radius=30, bsp=0., h_factor=0.)
    center_slice = grids['reflectivity_horizontal'][1, 4, :]
    assert_array_equal(np.round(center_slice), EXPECTED_CENTER_SLICE)


def test_map_to_grid_masked_refl_field():
    radar = pyart.testing.make_target_radar()

    # mask the last gate of the first ray
    fdata = radar.fields['reflectivity_horizontal']['data']
    fdata = np.ma.masked_invalid(fdata)
    fdata.mask[0, -1] = True
    radar.fields['reflectivity_horizontal']['data'] = fdata

    grids = pyart.map.map_to_grid((radar,),
                                  **COMMON_MAP_TO_GRID_ARGS)
    center_slice = grids['reflectivity_horizontal'][1, 4, :]
    assert_array_equal(np.round(center_slice), EXPECTED_CENTER_SLICE)


def test_map_to_grid_no_copy():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_to_grid((radar,), copy_field_data=False,
                                  **COMMON_MAP_TO_GRID_ARGS)
    center_slice = grids['reflectivity_horizontal'][1, 4, :]
    assert_array_equal(np.round(center_slice), EXPECTED_CENTER_SLICE)


def test_map_to_grid_balltree():
    radar = pyart.testing.make_target_radar()
    grids = pyart.map.map_to_grid((radar,), algorithm='ball_tree',
                                  **COMMON_MAP_TO_GRID_ARGS)
    center_slice = grids['reflectivity_horizontal'][1, 4, :]
    assert_array_equal(np.round(center_slice), EXPECTED_CENTER_SLICE)


def test_grid_from_radars():
    radar = pyart.testing.make_target_radar()
    grid = pyart.map.grid_from_radars((radar,), **COMMON_MAP_TO_GRID_ARGS)

    # check field data
    center_slice = grid.fields['reflectivity_horizontal']['data'][1, 4, :]
    assert_array_equal(np.round(center_slice), EXPECTED_CENTER_SLICE)

    # check other Grid object attributes
    assert 'ROI' in grid.fields
    assert np.all(grid.fields['ROI']['data'] == 30.)
    assert_array_equal(grid.axes['x_disp']['data'],
                       np.linspace(-900, 900, 10))
    assert_array_equal(grid.axes['y_disp']['data'],
                       np.linspace(-900, 900, 9).astype('float64'))
    assert_array_equal(grid.axes['z_disp']['data'],
                       np.linspace(-400, 400, 3).astype('float64'))


def test_grid_from_radars_grid_origin():
    radar = pyart.testing.make_target_radar()
    grid = pyart.map.grid_from_radars((radar,), grid_origin=(36.4, -97.6),
                                      **COMMON_MAP_TO_GRID_ARGS)
    print round(grid.axes['lat']['data'][0], 2)
    print round(grid.axes['lon']['data'][0], 2)
    assert round(grid.axes['lat']['data'][0], 2) == 36.4
    assert round(grid.axes['lon']['data'][0], 2) == -97.6
