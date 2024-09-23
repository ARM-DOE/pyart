""" Unit Tests for Py-ART's retrieve/qvp.py module. """

import datetime

import numpy as np
import pytest
from netCDF4 import num2date
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


def test_project_to_vertical():
    z = np.array([100, 200, 500, 1000, 2000])
    data = np.ma.array(np.linspace(0, 50, len(z)))
    znew = np.arange(0, 2300, 100)
    data_out = pyart.retrieve.qvp.project_to_vertical(data, z, znew)
    data_out2 = pyart.retrieve.qvp.project_to_vertical(
        data, z, znew, interp_kind="nearest"
    )
    ref_out = np.ma.array(
        data=[
            np.nan,
            0.0,
            12.5,
            np.nan,
            np.nan,
            25.0,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            37.5,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            50.0,
            np.nan,
            np.nan,
        ],
        mask=[
            True,
            False,
            False,
            True,
            True,
            False,
            True,
            True,
            True,
            True,
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
            False,
            True,
            True,
        ],
    )
    ref_out2 = np.ma.array(
        data=[
            np.nan,
            0.0,
            12.5,
            12.5,
            25.0,
            25.0,
            25.0,
            25.0,
            37.5,
            37.5,
            37.5,
            37.5,
            37.5,
            37.5,
            37.5,
            37.5,
            50.0,
            50.0,
            50.0,
            50.0,
            50.0,
            np.nan,
            np.nan,
        ],
        mask=[
            True,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
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
        ],
    )

    assert np.all(data_out == ref_out)
    assert np.all(data_out2 == ref_out2)


def test_find_rng_index():
    vec = np.arange(0, 1000, 10)
    idx = pyart.retrieve.qvp.find_rng_index(vec, 140)
    assert idx == 14


def test_get_target_elevations():
    test_radar = pyart.testing.make_target_radar()
    test_radar.elevation["data"][0] = 2
    target_elevations, el_tol = pyart.retrieve.qvp.get_target_elevations(test_radar)
    assert target_elevations[-1] == 2
    assert np.all(target_elevations[0:-1] == 0.75)
    assert el_tol == 0


def test_find_nearest_gate(test_radar):
    ind_ray, ind_rng, azi, rng = pyart.retrieve.qvp.find_nearest_gate(
        test_radar, 36.4, -97.4
    )

    assert ind_ray == 141.0
    assert ind_rng == 145.0
    assert azi == 141.0
    assert (
        abs(rng - 14514.514) < 1e-3
    )  # Allow for a small tolerance in floating-point comparison


def test_find_neighbour_gates(test_radar):
    inds_ray, inds_rng = pyart.retrieve.qvp.find_neighbour_gates(
        test_radar, 141, 14514, delta_azi=3, delta_rng=200
    )

    assert np.all(inds_ray == [139, 140, 141, 142, 143])
    assert np.all(inds_rng == [143, 144, 145, 146])


def test_find_ang_index():
    ang_vec = np.arange(0, 360)
    idx = pyart.retrieve.qvp.find_ang_index(ang_vec, 36)
    assert idx == 36


def test__create_qvp_object(test_radar):
    qvp = pyart.retrieve.qvp._create_qvp_object(
        test_radar, ["reflectivity"], hmax=5000, hres=500
    )
    assert np.all(qvp.range["data"] == np.arange(250, 5000, 500))
    assert "reflectivity" in qvp.fields
    assert qvp.sweep_mode["data"] == ["qvp"]
    assert qvp.fixed_angle["data"] == 10


def test__create_along_coord_object(test_radar):
    acoord = pyart.retrieve.qvp._create_along_coord_object(
        test_radar, ["reflectivity"], np.arange(0, 10000, 100), 10, "ALONG_AZI"
    )

    assert np.all(acoord.range["data"] == np.arange(0, 10000, 100))
    assert acoord.sweep_mode["data"] == ["ALONG_AZI"]


def test__update_qvp_metadata(test_radar):
    qvp = pyart.retrieve.qvp._create_qvp_object(
        test_radar, ["reflectivity"], hmax=5000, hres=500
    )
    new_time = datetime.datetime(2024, 6, 10)
    qvp = pyart.retrieve.qvp._update_qvp_metadata(qvp, new_time, 10, 45)
    start_time = num2date(0, qvp.time["units"], qvp.time["calendar"])
    time_offset = (new_time - start_time).total_seconds()

    assert np.all(qvp.gate_longitude["data"] == 10)
    assert np.all(qvp.gate_latitude["data"] == 45)
    assert qvp.time["data"] == time_offset


def test__update_along_coord_metadata(test_radar):
    acoord = pyart.retrieve.qvp._create_along_coord_object(
        test_radar, ["reflectivity"], np.arange(0, 10000, 100), 10, "ALONG_AZI"
    )
    new_time = datetime.datetime(2024, 6, 10)
    acoord = pyart.retrieve.qvp._update_along_coord_metadata(acoord, new_time, 5, 200)

    start_time = num2date(0, acoord.time["units"], acoord.time["calendar"])
    time_offset = (new_time - start_time).total_seconds()

    assert acoord.time["data"] == time_offset
    assert acoord.elevation["data"] == [5]
    assert acoord.azimuth["data"] == [200]
