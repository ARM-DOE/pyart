""" Unit tests for the circular_stats.py module. """

import numpy as np
from numpy.testing import assert_almost_equal

from pyart.util import circular_stats


def test_compute_directional_stats():
    field_1d = np.ma.array([1, 1, 3, 4, 5])
    field_1d.mask = [True, False, False, False, False]
    field = np.tile(field_1d, (10, 1))

    mean, nvalid = circular_stats.compute_directional_stats(field, axis=1)
    median, nvalid = circular_stats.compute_directional_stats(
        field, axis=1, avg_type="median"
    )

    assert np.all(mean == 3.25)
    assert np.all(median == 3.5)
    assert np.all(nvalid == 4)

    # Test with larger nb of nvalid_min
    mean, nvalid = circular_stats.compute_directional_stats(field, axis=1, nvalid_min=5)
    assert np.all(mean.mask)

    # Test other direction
    mean, nvalid = circular_stats.compute_directional_stats(field, axis=0)
    median, nvalid = circular_stats.compute_directional_stats(
        field, axis=0, avg_type="median"
    )

    ref = np.ma.array([np.nan, 1, 3, 4, 5], mask=[1, 0, 0, 0, 0])
    assert np.all(mean == ref)
    assert np.all(median == ref)
    assert np.all(nvalid == [0, 10, 10, 10, 10])


def test_mean_of_two_angles():
    mean = circular_stats.mean_of_two_angles(np.pi / 4 + 0.2, 7 * np.pi / 4)
    assert_almost_equal(mean, 0.10, 2)


def test_mean_of_two_angles_deg():
    mean = circular_stats.mean_of_two_angles_deg(359.9, 0.5)
    assert_almost_equal(mean, 0.20, 2)


def test_angles_radians():
    dist = [2 * np.pi - 0.1, 0, 0.1]
    mean = circular_stats.angular_mean(dist)
    assert_almost_equal(mean, 0.0, 2)
    std = circular_stats.angular_std(dist)
    assert_almost_equal(std, 0.08, 2)


def test_angles_degrees():
    dist = [350, 0, 10, 20, 30]
    mean = circular_stats.angular_mean_deg(dist)
    assert_almost_equal(mean, 10.0, 0)
    std = circular_stats.angular_std_deg(dist)
    assert_almost_equal(std, 14.0, 0)


def test_angles_interval():
    dist = [7.5, 8.5, 9.5, 9.5, -9.5, -8.5]
    mean = circular_stats.interval_mean(dist, -10, 10)
    assert_almost_equal(mean, 9.5, 1)
    std = circular_stats.interval_std(dist, -10, 10)
    assert_almost_equal(std, 1.3, 1)
