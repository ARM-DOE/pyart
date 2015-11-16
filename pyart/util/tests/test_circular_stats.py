""" Unit tests for the circular_stats.py module. """

from pyart.util import circular_stats

import numpy as np
from numpy.testing import assert_almost_equal


def test_mean_of_two_angles():
    mean = circular_stats.mean_of_two_angles(np.pi/4 + 0.2, 7 * np.pi / 4)
    assert_almost_equal(mean, 0.10, 2)


def test_mean_of_two_angles_deg():
    mean = circular_stats.mean_of_two_angles_deg(359.9, 0.5)
    assert_almost_equal(mean, 0.20, 2)


def test_angles_radians():
    dist = [2*np.pi-0.1, 0, 0.1]
    mean = circular_stats.angular_mean(dist)
    assert_almost_equal(mean, 0., 2)
    std = circular_stats.angular_std(dist)
    assert_almost_equal(std, 0.08, 2)


def test_angles_degrees():
    dist = [350, 0, 10, 20, 30]
    mean = circular_stats.angular_mean_deg(dist)
    assert_almost_equal(mean, 10., 0)
    std = circular_stats.angular_std_deg(dist)
    assert_almost_equal(std, 14., 0)


def test_angles_interval():
    dist = [7.5, 8.5, 9.5, 9.5, -9.5, -8.5]
    mean = circular_stats.interval_mean(dist, -10, 10)
    assert_almost_equal(mean, 9.5, 1)
    std = circular_stats.interval_std(dist, -10, 10)
    assert_almost_equal(std, 1.3, 1)
