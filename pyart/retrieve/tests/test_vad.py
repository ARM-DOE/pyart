""" Unit Tests for Py-ART's retrieve/vad.py module. """

import numpy as np
from numpy.testing import assert_almost_equal

import pyart


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
    z_want = np.linspace(0, 10, 5)

    vad_height = ([0., 2.5, 5., 7.5, 10.])
    vad_speed = ([4.999, 4.944, 4.886,
                  4.818, 4.749])
    vad_direction = ([89.871, 90.511, 91.208,
                      92.062, 92.956])
    u_wind = ([-4.999, -4.944, -4.885, -4.815, -4.743])
    v_wind = ([-0.011, 0.044, 0.103, 0.173, 0.245])

    vad = pyart.retrieve.velocity_azimuth_display(
        test_radar, velocity, z_want)

    assert_almost_equal(vad.height, vad_height, 3)
    assert_almost_equal(vad.speed, vad_speed, 3)
    assert_almost_equal(vad.direction, vad_direction, 3)
    assert_almost_equal(vad.u_wind, u_wind, 3)
    assert_almost_equal(vad.v_wind, v_wind, 3)
