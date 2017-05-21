""" Unit Tests for Py-ART's retrieve/vad.py module. """

from __future__ import print_function

import numpy as np
import pyart

from numpy.testing import assert_almost_equal


def test_velocity_azimuth_display():
    test_radar = pyart.testing.make_target_radar()
    height = np.arange(0, 1000, 200)
    speed = np.ones_like(height) * 5
    direction = np.array([0, 90, 180, 270, 45])
    profile = pyart.core.HorizontalWindProfile(height, speed, direction)
    sim_vel = pyart.util.simulated_vel_from_profile(test_radar, profile)
    test_radar.add_field('velocity', sim_vel, replace_existing=True)

    velocity = 'velocity'
    z_start = 0
    z_end = 10
    z_count = 5

    vad_height = ([0., 2.5, 5., 7.5, 10.])
    vad_speed = ([4.986, 4.940, 4.881, 4.819, 4.758])
    vad_direction = ([359.846, 359.302, 358.586, 357.810, 357.013])
    u_wind = ([0.013, 0.060, 0.120, 0.184, 0.247])
    v_wind = ([-4.986, -4.939, -4.879, -4.815, -4.752])

    vad = pyart.retrieve.velocity_azimuth_display(test_radar,
                                                  velocity,
                                                  z_start, z_end,
                                                  z_count)
    assert_almost_equal(vad.height, vad_height, 3)
    assert_almost_equal(vad.speed, vad_speed, 3)
    assert_almost_equal(vad.direction, vad_direction, 3)
    assert_almost_equal(vad.u_wind, u_wind, 3)
    assert_almost_equal(vad.v_wind, v_wind, 3)
