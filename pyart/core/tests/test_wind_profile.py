""" Unit Tests for Py-ART's core/wind_profile.py module. """

from __future__ import print_function

import numpy as np
from numpy.testing import assert_almost_equal, assert_raises

from pyart.core import HorizontalWindProfile


def test_horizontalwindprofile_class():
    height = np.arange(5)
    speed = np.ones(5)
    direction = np.array([0, 90, 180, 270, 45])
    hprofile = HorizontalWindProfile(height, speed, direction)

    assert_almost_equal(hprofile.height, [0, 1, 2, 3, 4])
    assert_almost_equal(hprofile.speed, [1, 1, 1, 1, 1])
    assert_almost_equal(hprofile.direction, [0, 90, 180, 270, 45])
    assert_almost_equal(hprofile.u_wind, [0, -1, 0, 1, -0.707], 3)
    assert_almost_equal(hprofile.v_wind, [-1, 0, 1, 0, -0.707], 3)


def test_horizontalwindprofile_from_u_and_v():
    height = np.arange(5)
    u_wind = [0, -1, 0, 1, -np.sqrt(2)/2.]
    v_wind = [-1, 0, 1, 0, -np.sqrt(2)/2.]
    hprofile = HorizontalWindProfile.from_u_and_v(height, u_wind, v_wind)

    assert_almost_equal(hprofile.height, [0, 1, 2, 3, 4])
    assert_almost_equal(hprofile.speed, [1, 1, 1, 1, 1])
    assert_almost_equal(hprofile.direction, [0, 90, 180, 270, 45])
    assert_almost_equal(hprofile.u_wind, [0, -1, 0, 1, -0.707], 3)
    assert_almost_equal(hprofile.v_wind, [-1, 0, 1, 0, -0.707], 3)


def test_horizontalwindprofile_error():
    height = [1, 2]
    speed = [1, 1]
    direction = [1, 2, 3]
    assert_raises(ValueError, HorizontalWindProfile, height, speed, direction)
