""" Unit Tests for Py-ART's util/radar_utils.py module. """

import numpy as np

import pyart


def test_is_vpt():
    radar = pyart.testing.make_empty_ppi_radar(10, 36, 3)
    assert not pyart.util.is_vpt(radar)
    pyart.util.to_vpt(radar)
    assert pyart.util.is_vpt(radar)


def test_join_radar():
    radar1 = pyart.testing.make_empty_ppi_radar(10, 36, 3)
    radar1.instrument_parameters = {"nyquist_velocity": {"data": np.array([6] * 3)}}
    field = {"data": np.ones((36 * 3, 10))}
    radar1.add_field("f1", field)
    radar1.add_field("f2", field)
    radar2 = pyart.testing.make_empty_ppi_radar(10, 36, 4)
    radar2.instrument_parameters = {"nyquist_velocity": {"data": np.array([8] * 4)}}
    field = {"data": np.ones((36 * 4, 10))}
    radar2.add_field("f1", field)

    radar3 = pyart.util.join_radar(radar1, radar2)
    assert "f1" in radar3.fields
    assert len(radar3.instrument_parameters["nyquist_velocity"]["data"]) == 7


def test_to_vpt():
    # single scan
    radar = pyart.testing.make_empty_ppi_radar(10, 36, 3)
    radar.instrument_parameters = {"prt_mode": {"data": np.array(["fixed"] * 3)}}
    pyart.util.to_vpt(radar)
    assert pyart.util.is_vpt(radar)
    assert radar.nsweeps == 1
    assert radar.azimuth["data"][10] == 0.0
    assert radar.elevation["data"][0] == 90.0
    assert len(radar.instrument_parameters["prt_mode"]["data"]) == 1

    # multiple scans
    radar = pyart.testing.make_empty_ppi_radar(10, 36, 3)
    radar.instrument_parameters = {"prt_mode": {"data": np.array(["fixed"] * 3)}}
    pyart.util.to_vpt(radar, False)
    assert pyart.util.is_vpt(radar)
    assert radar.nsweeps == 108
    assert radar.azimuth["data"][10] == 10.0
    assert radar.elevation["data"][0] == 90.0
    assert len(radar.instrument_parameters["prt_mode"]["data"]) == 108


def test_determine_sweeps():
    # ppi
    radar = pyart.testing.make_empty_ppi_radar(10, 36, 3)
    radar.elevation["data"] = radar.elevation["data"] * np.ceil(
        np.arange(1, 36 * 3 + 1) / 36
    )
    pyart.util.determine_sweeps(radar)
    assert np.all(radar.sweep_end_ray_index["data"] == [35, 71, 107])
    assert np.all(radar.sweep_start_ray_index["data"] == [0, 36, 72])
    assert len(radar.sweep_number["data"]) == 3
    assert radar.nsweeps == 3

    # rhi
    radar = pyart.testing.make_empty_rhi_radar(10, 25, 5)
    radar.azimuth["data"] = (
        radar.azimuth["data"] * np.ceil(np.arange(1, 25 * 5 + 1) / 25) * 25
    )
    pyart.util.determine_sweeps(radar)
    assert np.all(radar.sweep_end_ray_index["data"] == [24, 49, 74, 99, 124])
    assert np.all(radar.sweep_start_ray_index["data"] == [0, 25, 50, 75, 100])
    assert len(radar.sweep_number["data"]) == 5
    assert radar.nsweeps == 5


def test_subset_radar():
    radar = pyart.testing.make_empty_ppi_radar(10, 36, 3)
    field = {"data": np.ones((36 * 3, 10))}
    radar.add_field("f1", field)
    radar.add_field("f2", field)
    azi_min = 10
    azi_max = 100
    rng_min = 200
    rng_max = 800
    ele_min = 0.75
    ele_max = 0.75
    radarcut = pyart.util.radar_utils.subset_radar(
        radar,
        ["f1"],
        rng_min=200,
        rng_max=800,
        ele_min=0.75,
        ele_max=0.75,
        azi_min=10,
        azi_max=100,
    )
    # assert correct domain and correct fields
    assert radarcut.azimuth["data"].min() >= azi_min
    assert radarcut.azimuth["data"].max() <= azi_max
    assert radarcut.range["data"].min() >= rng_min
    assert radarcut.range["data"].max() <= rng_max
    assert radarcut.elevation["data"].min() >= ele_min
    assert radarcut.elevation["data"].max() <= ele_max
    assert list(radarcut.fields) == ["f1"]


def test_ma_broadcast_to():
    buf = np.ma.zeros(5)
    buf.mask = [1, 1, 0, 0, 0]
    buf_broad = pyart.util.radar_utils.ma_broadcast_to(buf, (10, 5))
    assert buf_broad.shape == (10, 5)
    assert buf_broad.mask.shape == (10, 5)
    expected_mask = np.tile(np.array([True, True, False, False, False]), [10, 1])
    assert np.all(expected_mask == buf_broad.mask)


# read in example file
radar = pyart.io.read_nexrad_archive(pyart.testing.NEXRAD_ARCHIVE_MSG31_FILE)


def test_image_mute_radar():
    # image mute example file
    mute_radar = pyart.util.image_mute_radar(
        radar,
        field="reflectivity",
        mute_field="cross_correlation_ratio",
        mute_threshold=0.97,
    )

    # check that fields are added to the radar object
    assert "nonmuted_reflectivity" in mute_radar.fields.keys()
    assert "muted_reflectivity" in mute_radar.fields.keys()

    # check that the number of points in the muted and non muted reflectivity
    # fields have the same number of points as the rhoHV field
    n_rhohv = np.sum(~mute_radar.fields["cross_correlation_ratio"]["data"].mask)
    n_mutez = np.sum(~mute_radar.fields["muted_reflectivity"]["data"].mask)
    n_nonmutez = np.sum(~mute_radar.fields["nonmuted_reflectivity"]["data"].mask)

    assert n_mutez + n_nonmutez == n_rhohv
