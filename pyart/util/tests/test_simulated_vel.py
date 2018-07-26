""" Unit tests for Py-ART's util/simulated_vel.py module. """

import numpy as np
from numpy.testing import assert_almost_equal
import pytest

import pyart

@pytest.mark.skip(
    reason=("Test does not fail locally, but does when "
            "pytest is ran on all of Py-ART. Can't find issue."))
def test_simulated_vel():

    radar = pyart.testing.make_target_radar()

    # profile of 10 m/s winds out of the west at all heights
    height = np.arange(0, 1000, 30)
    speed = np.ones_like(height) * 10.
    direction = np.ones_like(height) * 270.
    profile = pyart.core.HorizontalWindProfile(height, speed, direction)

    sim_vel = pyart.util.simulated_vel_from_profile(radar, profile)

    assert sim_vel['data'].shape == (360, 50)

    # check simulated velocities along the north, east, south and west rays
    assert_almost_equal(sim_vel['data'][[0, 90, 180, 270], 10],
                        [0, 10, 0, -10], 0)
