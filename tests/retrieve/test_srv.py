""" Unit Tests for Py-ART's retrieve/srv.py module. """

import numpy as np
import pytest

import pyart


def test_storm_relative_velocity():
    # Test all sweeps
    radar = pyart.io.read_nexrad_archive(
        "s3://noaa-nexrad-level2/2022/03/31/KDGX/KDGX20220331_012808_V06"
    )
    sr_data = pyart.retrieve.storm_relative_velocity(radar, direction="NE", speed=20.0)
    test_data = [2.276456117630005, 4.776455879211426, 3.276456117630005]
    np.testing.assert_almost_equal(sr_data[-1][0:3].tolist(), test_data)

    # Test one sweep
    radar_sweep = radar.extract_sweeps([21])
    sr_one_sweep = pyart.retrieve.storm_relative_velocity(
        radar_sweep, direction="NE", speed=20.0
    )
    one_sweep_data = [0.1278250813484192, 4.6278252601623535, 4.6278252601623535]
    np.testing.assert_almost_equal(sr_one_sweep[0][0:3].tolist(), one_sweep_data)

    # Test missing parameters
    pytest.raises(ValueError, pyart.retrieve.storm_relative_velocity, radar)
    pytest.raises(
        ValueError, pyart.retrieve.storm_relative_velocity, radar, u=14.14, v=None
    )
