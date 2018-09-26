""" Unit Tests for Py-ART's retrieve/vad.py module. """

import numpy as np
from numpy.testing import assert_allclose

import pyart


def test_velocity_azimuth_display():
    test_radar = pyart.testing.make_target_radar()
    height = np.arange(0.0, 1000.0, 200.0)
    speed = np.ones_like(height) * 5.0
    direction = np.array([0.0, 90.0, 180.0, 270.0, 45.0])
    profile = pyart.core.HorizontalWindProfile(height, speed, direction)
    sim_vel = pyart.util.simulated_vel_from_profile(test_radar, profile)
    test_radar.add_field('velocity', sim_vel,
                         replace_existing=True)

    velocity = 'velocity'
    z_want = np.linspace(0.0, 10.0, 5)

    vad_height = ([0.0, 2.5, 5.0, 7.5, 10.0])
    vad_speed = ([4.9997, 4.9445, 4.8865,
                  4.8182, 4.7497])
    vad_direction = ([89.8712, 90.5113, 91.2084,
                      92.0622, 92.9569])
    u_wind = ([-4.9997, -4.9443,
               -4.8854, -4.8150, -4.7434])
    v_wind = ([-0.0112, 0.0441,
               0.1030, 0.1733, 0.2450])

    vad = pyart.retrieve.velocity_azimuth_display(
        test_radar, velocity, z_want)

    assert_allclose(vad.height, vad_height, rtol=1e-3, atol=1e-1)
    assert_allclose(vad.speed, vad_speed, rtol=1e-3, atol=1e-1)
    assert_allclose(vad.direction, vad_direction,rtol=1e-3, atol=1e-1)
    assert_allclose(vad.u_wind, u_wind, rtol=1e-3, atol=1e-1)
    assert_allclose(vad.v_wind, v_wind, rtol=1e-3, atol=1e-1)
