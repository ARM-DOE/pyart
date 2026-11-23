"""Unit Tests for Py-ART's util/radar_utils.py module."""

import numpy as np
import pytest

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


@pytest.fixture(scope="module")
def msg31_radar():
    return pyart.io.read_nexrad_archive(pyart.testing.NEXRAD_ARCHIVE_MSG31_FILE)


def _original_sweep_indices(radar, selected_radar):
    """Return source sweep indices retained in the extracted Radar."""
    return selected_radar.sweep_number["data"]


def test_select_unique_elevation_sweeps_returns_radar(msg31_radar):
    selected = pyart.util.select_unique_elevation_sweeps(msg31_radar)

    assert selected is not msg31_radar
    assert selected.nsweeps == 14
    np.testing.assert_array_equal(
        _original_sweep_indices(msg31_radar, selected),
        [0, 2, *range(4, 16)],
    )


def test_select_unique_elevation_sweeps_keeps_first_representative(msg31_radar):
    selected = pyart.util.select_unique_elevation_sweeps(
        msg31_radar,
        strategy="first",
    )

    np.testing.assert_array_equal(
        _original_sweep_indices(msg31_radar, selected),
        [0, 2, *range(4, 16)],
    )


def test_select_unique_elevation_sweeps_keeps_last_representative(msg31_radar):
    selected = pyart.util.select_unique_elevation_sweeps(
        msg31_radar,
        strategy="last",
    )

    np.testing.assert_array_equal(
        _original_sweep_indices(msg31_radar, selected),
        [1, 3, *range(4, 16)],
    )


@pytest.mark.parametrize(
    "waveform_type, expected_sweeps",
    [
        ("contiguous_surveillance", [0, 2]),
        ("contiguous_doppler_with_ambiguity_resolution", [1, 3]),
        ("batch", [4, 5, 6, 7, 8]),
        ("contiguous_doppler_without_ambiguity_resolution", list(range(9, 16))),
    ],
)
def test_select_unique_elevation_sweeps_filters_waveform(
    msg31_radar, waveform_type, expected_sweeps
):
    selected = pyart.util.select_unique_elevation_sweeps(
        msg31_radar,
        waveform_type=waveform_type,
    )

    np.testing.assert_array_equal(
        _original_sweep_indices(msg31_radar, selected),
        expected_sweeps,
    )


def test_select_unique_elevation_sweeps_filters_both_cd_variants(msg31_radar):
    selected = pyart.util.select_unique_elevation_sweeps(
        msg31_radar,
        waveform_type=[
            "contiguous_doppler_with_ambiguity_resolution",
            "contiguous_doppler_without_ambiguity_resolution",
        ],
    )

    np.testing.assert_array_equal(
        _original_sweep_indices(msg31_radar, selected),
        [1, 3, *range(9, 16)],
    )


def test_select_unique_elevation_sweeps_preserves_sweep_metadata(msg31_radar):
    selected = pyart.util.select_unique_elevation_sweeps(
        msg31_radar,
        waveform_type=[
            "contiguous_doppler_without_ambiguity_resolution",
            "contiguous_doppler_with_ambiguity_resolution",
            "batch",
        ],
    )

    expected_indices = np.array([1, 3, *range(4, 16)])
    np.testing.assert_array_equal(
        selected.sweep_number["data"],
        expected_indices,
    )

    for key, metadata in msg31_radar.sweep_metadata.items():
        np.testing.assert_array_equal(
            selected.sweep_metadata[key]["data"],
            metadata["data"][expected_indices],
        )


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"atol": -0.1}, "atol must be non-negative"),
        ({"strategy": "middle"}, "strategy must be"),
    ],
)
def test_select_unique_elevation_sweeps_rejects_invalid_options(
    msg31_radar, kwargs, message
):
    with pytest.raises(ValueError, match=message):
        pyart.util.select_unique_elevation_sweeps(msg31_radar, **kwargs)
