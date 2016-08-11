""" Unit Tests for Py-ART's retrieve/advection.py module. """

from numpy.testing import assert_almost_equal
import pyart


def test_grid_displacement_pc():
    grid1 = pyart.testing.make_storm_grid()
    data = grid1.fields['reflectivity']['data'].copy()
    grid1.fields['reflectivity']['data'] = data[:, 5:-5, 3:-3].copy()

    grid2 = pyart.testing.make_storm_grid()
    grid2.fields['reflectivity']['data'] = data[:, 0:-10, 0:-6].copy()

    # test pixels
    displacement = pyart.retrieve.grid_displacement_pc(
        grid1, grid2, 'reflectivity', 0)
    assert displacement == (-5, -3)

    # test distance
    grid1.fields['reflectivity']['valid_min'] = 0
    grid2.fields['reflectivity']['valid_min'] = 0
    displacement = pyart.retrieve.grid_displacement_pc(
        grid1, grid2, 'reflectivity', 0, return_value='distance')
    dx = grid1.x['data'][1] - grid1.x['data'][0]
    dy = grid1.y['data'][1] - grid1.y['data'][0]
    assert_almost_equal(displacement[0], -5*dy, 0)
    assert_almost_equal(displacement[1], -3*dx, 0)

    # test velocity
    grid2.time['data'] += 1
    displacement = pyart.retrieve.grid_displacement_pc(
        grid1, grid2, 'reflectivity', 0, return_value='velocity')
    assert_almost_equal(displacement[0], -5*dy, 0)
    assert_almost_equal(displacement[1], -3*dx, 0)

    # test fallback
    displacement = pyart.retrieve.grid_displacement_pc(
        grid1, grid2, 'reflectivity', 0, return_value='foobar')
    assert displacement == (-5, -3)


def test_grid_shift():

    # create two guassian storms
    grid1 = pyart.testing.make_normal_storm(10.0, [0.0, 0.0])
    grid2 = pyart.testing.make_normal_storm(10.0, [5.0, 5.0])

    # trim one, trim and shift the other
    trimmed_grid2 = pyart.retrieve.grid_shift(grid2, [0.0, 0.0], trim_edges=10)
    shifted_grid1 = pyart.retrieve.grid_shift(grid1, [5.0, 5.0], trim_edges=10)

    # The difference should be nearly zero
    data1 = shifted_grid1.fields['reflectivity']['data'][0]
    data2 = trimmed_grid2.fields['reflectivity']['data'][0]
    assert (data1 - data2).mean() < 1.e-10
