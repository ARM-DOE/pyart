"""Unit Tests for Py-ART's retrieve/vad.py module."""

import numpy as np
from numpy.testing import assert_allclose

import pyart


def test_vad_michelson():
    test_radar = pyart.testing.make_target_radar()
    height = np.arange(0.0, 1000.0, 200.0)
    speed = np.ones_like(height) * 5.0
    direction = np.array([0.0, 90.0, 180.0, 270.0, 45.0])
    profile = pyart.core.HorizontalWindProfile(height, speed, direction)
    sim_vel = pyart.util.simulated_vel_from_profile(test_radar, profile)
    test_radar.add_field("velocity", sim_vel, replace_existing=True)

    velocity = "velocity"
    z_want = np.linspace(0.0, 10.0, 5)

    vad_height = [0.0, 2.5, 5.0, 7.5, 10.0]
    vad_speed = [5.0, 4.9483, 4.8888, 4.8260, 4.7667]
    vad_direction = [89.9999, 90.6011, 91.3184, 92.1049, 92.8771]
    u_wind = [-5.0, -4.9481, -4.8875, -4.8227, -4.7607]
    v_wind = [0.0, 0.0519, 0.1125, 0.1773, 0.2393]

    vad = pyart.retrieve.vad_michelson(test_radar, velocity, z_want)

    assert_allclose(vad.height, vad_height, rtol=1e-3, atol=1e-1)
    assert_allclose(vad.speed, vad_speed, rtol=1e-3, atol=1e-1)
    assert_allclose(vad.direction, vad_direction, rtol=1e-3, atol=1e-1)
    assert_allclose(vad.u_wind, u_wind, rtol=1e-3, atol=1e-1)
    assert_allclose(vad.v_wind, v_wind, rtol=1e-3, atol=1e-1)


def test_vad_michelson_masked_gates():
    # Masked gates must be excluded, not treated as zero velocity (#1488)
    test_radar = pyart.testing.make_target_radar()
    height = np.arange(0.0, 1000.0, 200.0)
    speed = np.ones_like(height) * 5.0
    direction = np.ones_like(height) * 90.0
    profile = pyart.core.HorizontalWindProfile(height, speed, direction)
    sim_vel = pyart.util.simulated_vel_from_profile(test_radar, profile)
    rng = np.random.default_rng(1488)
    mask = rng.random(sim_vel["data"].shape) < 0.3
    sim_vel["data"] = np.ma.masked_array(sim_vel["data"], mask)
    test_radar.add_field("velocity", sim_vel, replace_existing=True)

    z_want = np.linspace(0.0, 10.0, 5)
    vad = pyart.retrieve.vad_michelson(test_radar, "velocity", z_want)

    assert_allclose(vad.speed, 5.0, atol=0.05)
    assert_allclose(vad.u_wind, -5.0, atol=0.05)
    assert_allclose(vad.v_wind, 0.0, atol=0.05)


def test_vad_michelson_valid_ray_min():
    # Gates with fewer valid rays than valid_ray_min are left out
    test_radar = pyart.testing.make_target_radar()
    height = np.arange(0.0, 1000.0, 200.0)
    speed = np.ones_like(height) * 5.0
    direction = np.ones_like(height) * 90.0
    profile = pyart.core.HorizontalWindProfile(height, speed, direction)
    sim_vel = pyart.util.simulated_vel_from_profile(test_radar, profile)
    mask = np.ones(sim_vel["data"].shape, dtype=bool)
    mask[::36] = False  # 10 valid rays per gate
    sim_vel["data"] = np.ma.masked_array(sim_vel["data"], mask)
    test_radar.add_field("velocity", sim_vel, replace_existing=True)

    z_want = np.linspace(0.0, 10.0, 5)
    vad = pyart.retrieve.vad_michelson(test_radar, "velocity", z_want)
    assert np.all(np.isnan(vad.speed))

    vad = pyart.retrieve.vad_michelson(test_radar, "velocity", z_want, valid_ray_min=5)
    assert_allclose(vad.speed, 5.0, atol=0.05)


def test_vad_browning():
    test_radar = pyart.testing.make_target_radar()
    height = np.arange(0.0, 1000.0, 200.0)
    speed = np.ones_like(height) * 5.0
    direction = np.array([0.0, 90.0, 180.0, 270.0, 45.0])
    profile = pyart.core.HorizontalWindProfile(height, speed, direction)
    sim_vel = pyart.util.simulated_vel_from_profile(test_radar, profile)
    test_radar.add_field("velocity", sim_vel, replace_existing=True)

    velocity = "velocity"
    z_want = np.linspace(0.0, 10.0, 5)

    vad_height = [0.0, 2.5, 5.0, 7.5, 10.0]
    vad_speed = [4.9728, 4.9465, 4.8802, 4.82951, 4.7578]

    vad_direction = [90.3142, 90.6225, 91.4236, 92.0601, 92.99520]

    u_wind = [-4.9727, -4.9462, -4.8787, -4.8263, -4.7513]
    v_wind = [0.02727, 0.05374, 0.1212, 0.1736, 0.2486]

    vad = pyart.retrieve.vad_browning(test_radar, velocity, z_want)
    assert_allclose(vad.height, vad_height, rtol=1e-3, atol=1e-1)
    assert_allclose(vad.speed, vad_speed, rtol=1e-3, atol=1e-1)
    assert_allclose(vad.direction, vad_direction, rtol=1e-3, atol=1e-1)
    assert_allclose(vad.u_wind, u_wind, rtol=1e-3, atol=1e-1)
    assert_allclose(vad.v_wind, v_wind, rtol=1e-3, atol=1e-1)
