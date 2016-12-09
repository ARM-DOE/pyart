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
    test_radar.add_field('velocity', sim_vel,
                         replace_existing=True)
    velocity = 'velocity'
    z_start = 0
    z_end = 10
    z_count = 5

    vad_height = ([0., 2.5, 5., 7.5, 10.])
    vad_speed = [4.98665725, 4.94020686, 4.88107152,
                 4.81939374, 4.75851962]
    vad_direction = [359.84659496, 359.30240553, 358.58658589,
                     357.81073051, 357.01353486]
    u_wind = ([0.01335138, 0.06014712, 0.12039762,
               0.18410404, 0.24791911])
    v_wind = ([-4.98663937, -4.9398407, -4.87958641,
               -4.81587601, -4.75205693])

    vad = pyart.retrieve.velocity_azimuth_display(test_radar,
                                                  velocity,
                                                  z_start, z_end,
                                                  z_count)
    assert_almost_equal(vad.height, vad_height, 8)
    assert_almost_equal(vad.speed, vad_speed, 8)
    assert_almost_equal(vad.direction, vad_direction, 8)
    assert_almost_equal(vad.u_wind, u_wind, 8)
    assert_almost_equal(vad.v_wind, v_wind, 8)
