""" Unit Tests for Py-ART's retrieve/qvp.py module. """

import numpy as np
import pytest
from numpy.testing import assert_almost_equal

import pyart


@pytest.fixture
def test_radar():
    test_radar = pyart.testing.make_empty_ppi_radar(1000, 360, 1)
    test_radar.range["data"] *= 100
    test_radar.elevation["data"] = np.ones(test_radar.elevation["data"].shape) * 10.0
    test_radar.fixed_angle["data"] = np.array([10.0])
    refl = 0.1 * np.arange(test_radar.ngates)
    refl = np.tile(refl, (test_radar.nrays, 1))
    test_radar.add_field("reflectivity", {"data": refl})
    return test_radar


def test_compute_qvp(test_radar):
    qvp = pyart.retrieve.compute_qvp(test_radar, ["reflectivity"], hmax=10000)
    qvp_refl = np.array(
        [
            1.8999999999999875,
            2.200000000000011,
            2.3999999999999853,
            2.7000000000000175,
            3.0,
            3.299999999999978,
            3.599999999999994,
            3.9000000000000132,
            4.200000000000027,
            4.5,
            4.7000000000000295,
            5.0,
            5.299999999999969,
            5.599999999999963,
            5.900000000000039,
            6.2000000000000135,
            6.5,
            6.7,
            7.0,
            7.300000000000015,
        ]
    )

    assert np.allclose(qvp.fields["reflectivity"]["data"][0, 10:30], qvp_refl)
    assert len(qvp.range["data"]) == 200
    assert np.all(qvp.range["data"] == np.arange(25, 10000, 50))
    assert qvp.azimuth["data"][0] == 0


def test_compute_rqvp(test_radar):
    rqvp = pyart.retrieve.compute_rqvp(
        test_radar, ["reflectivity"], hres=50, hmax=5000, rmax=5000, weight_power=-1
    )
    rqvp_refl = np.array(
        [
            1.8999999999999875,
            2.200000000000011,
            2.3999999999999853,
            2.7000000000000175,
            3.0,
            3.299999999999978,
            3.599999999999994,
            3.9000000000000132,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ]
    )
    rqvp_refl_mask = np.array(
        [
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            True,
            True,
            True,
            True,
            True,
            True,
            True,
            True,
            True,
            True,
            True,
            True,
        ]
    )

    assert np.allclose(rqvp.fields["reflectivity"]["data"][0, 10:30], rqvp_refl)
    assert np.allclose(
        rqvp.fields["reflectivity"]["data"][0, 10:30].mask, rqvp_refl_mask
    )

    assert len(rqvp.range["data"]) == 100
    assert np.all(rqvp.range["data"] == np.arange(25, 5000, 50))
    assert rqvp.azimuth["data"][0] == 0


def test_compute_evp(test_radar):
    evp = pyart.retrieve.compute_evp(
        test_radar, ["reflectivity"], -97, 36, delta_rng=40000
    )

    evp_refl = np.array(
        [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            32.7,
            33.7,
            35.10000000000001,
            36.5,
            37.89999999999999,
            39.29999999999999,
            40.70000000000001,
        ]
    )

    evp_refl_mask = np.array(
        [
            [
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
            ]
        ]
    )

    assert np.allclose(evp.fields["reflectivity"]["data"][0, 10:30], evp_refl)
    assert np.allclose(evp.fields["reflectivity"]["data"][0, 10:30].mask, evp_refl_mask)

    assert len(evp.range["data"]) == 40
    assert np.all(evp.range["data"] == np.arange(125, 10000, 250))
    assert evp.azimuth["data"][0] == 0


def test_quasi_vertical_profile():
    test_radar = pyart.testing.make_target_radar()
    height = np.arange(0, 1000, 200)
    speed = np.ones_like(height) * 5
    direction = np.array([0, 90, 180, 270, 45])
    profile = pyart.core.HorizontalWindProfile(height, speed, direction)
    sim_vel = pyart.util.simulated_vel_from_profile(test_radar, profile)
    test_radar.add_field("velocity", sim_vel, replace_existing=True)

    qvp = pyart.retrieve.quasi_vertical_profile(test_radar)

    qvp_height = [
        0.0,
        0.0,
        0.0,
        1.0,
        1.0,
        1.0,
        1.0,
        2.0,
        2.0,
        2.0,
        2.0,
        3.0,
        3.0,
        3.0,
        3.0,
        3.0,
        4.0,
        4.0,
        4.0,
        4.0,
        5.0,
        5.0,
        5.0,
        5.0,
        6.0,
        6.0,
        6.0,
        7.0,
        7.0,
        7.0,
        7.0,
        8.0,
        8.0,
        8.0,
        8.0,
        9.0,
        9.0,
        9.0,
        10.0,
        10.0,
        10.0,
        10.0,
        11.0,
        11.0,
        11.0,
        11.0,
        12.0,
        12.0,
        12.0,
        13.0,
    ]

    qvp_range = [
        0.0,
        20.408,
        40.816,
        61.224,
        81.632,
        102.040,
        122.448,
        142.857,
        163.265,
        183.673,
        204.081,
        224.489,
        244.897,
        265.306,
        285.714,
        306.122,
        326.530,
        346.938,
        367.346,
        387.755,
        408.163,
        428.571,
        448.979,
        469.387,
        489.795,
        510.204,
        530.612,
        551.020,
        571.428,
        591.836,
        612.244,
        632.653,
        653.061,
        673.469,
        693.877,
        714.285,
        734.693,
        755.102,
        775.510,
        795.918,
        816.326,
        836.734,
        857.142,
        877.551,
        897.959,
        918.367,
        938.775,
        959.183,
        979.591,
        1000,
    ]

    qvp_reflectivity = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        10.0,
        20.0,
        20.0,
        20.0,
        20.0,
        20.0,
        20.0,
        20.0,
        20.0,
        20.0,
        20.0,
        30.0,
        30.0,
        30.0,
        30.0,
        30.0,
        30.0,
        30.0,
        30.0,
        30.0,
        30.0,
        40.0,
        40.0,
        40.0,
        40.0,
        40.0,
        40.0,
        40.0,
        40.0,
        40.0,
        40.0,
    ]

    assert_almost_equal(qvp["height"], qvp_height, 3)
    assert_almost_equal(qvp["range"], qvp_range, 3)
    assert_almost_equal(qvp["reflectivity"], qvp_reflectivity, 3)
