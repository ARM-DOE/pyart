""" Unit tests for bias and noise module. """

import dask
import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_almost_equal, assert_array_equal
from open_radar_data import DATASETS

import pyart

radar = pyart.io.read(pyart.testing.NEXRAD_ARCHIVE_MSG31_FILE)
kazr_file = DATASETS.fetch("sgpkazrgeC1.a1.20190529.000002.cdf")
radar_kazr = pyart.aux_io.read_kazr(kazr_file)


def test_correct_bias():
    corr_field_expected = [-32.0, -32.0, -32.0]
    corr_test = pyart.correct.correct_bias(radar, field_name=None)

    assert_allclose(corr_test["data"][0][0:3], corr_field_expected, atol=1e-14)
    assert_allclose(corr_test["data"][-1][0:3], corr_field_expected, atol=1e-14)
    assert corr_test["long_name"] == "Corrected reflectivity"

    radar.add_field("corrected_reflectivity", corr_test, replace_existing=True)

    corr_test_2 = pyart.correct.correct_bias(radar, field_name="corrected_reflectivity")
    assert corr_test_2["long_name"] == "Corrected reflectivity"

    bias = 0
    field_name = "foo"
    pytest.raises(KeyError, pyart.correct.correct_bias, radar, bias, field_name)


def test_calc_zdr_offset():
    xsapr_test_file = DATASETS.fetch("sgpxsaprcfrvptI4.a1.20200205.100827.nc")
    ds = pyart.io.read(xsapr_test_file)
    gatefilter = pyart.filters.GateFilter(ds)
    gatefilter.exclude_below("cross_correlation_ratio_hv", 0.995)
    gatefilter.exclude_above("cross_correlation_ratio_hv", 1)
    gatefilter.exclude_below("reflectivity", 10)
    gatefilter.exclude_above("reflectivity", 30)

    results = pyart.correct.calc_zdr_offset(
        ds,
        zdr_var="differential_reflectivity",
        gatefilter=gatefilter,
        height_range=(1000, 3000),
    )
    assert_almost_equal(results["bias"], 2.69, decimal=2)
    assert_almost_equal(results["profile_reflectivity"][15], 14.37, decimal=2)


def test_calc_noise_floor():
    expected = [-46.25460013, -46.48371626, -46.3314618, -46.82639895, -46.76403711]

    result = pyart.correct.calc_noise_floor(radar_kazr, "reflectivity_copol", "range")

    assert_almost_equal(result[0:5], expected, decimal=3)

    bad_radar = "foo"
    pytest.raises(
        ValueError,
        pyart.correct.calc_noise_floor,
        bad_radar,
        "reflectivity_copol",
        "range",
    )


def test_range_correction():
    expected_first = [-40.92386, -39.07553, -29.10794, -26.25786, -26.27254]
    expected_last = [-43.42107, -23.74393, -22.77576, -16.700004, -19.68343]

    result = pyart.correct.range_correction(radar_kazr, "reflectivity_copol", "range")
    assert_almost_equal(result[0][0:5], expected_first, decimal=3)
    assert_almost_equal(result[-1][0:5], expected_last, decimal=3)

    bad_radar = "foo"
    pytest.raises(
        ValueError,
        pyart.correct.range_correction,
        bad_radar,
        "reflectivity_copol",
        "range",
    )

    radar_kazr.range.pop("units")
    pytest.warns(
        UserWarning,
        pyart.correct.range_correction,
        radar_kazr,
        "reflectivity_copol",
        "range",
    )


def test_cloud_threshold():
    expected = [-46.25460, -46.48371, -46.33146, -46.82639, -46.76403]

    result = pyart.correct.range_correction(radar_kazr, "reflectivity_copol", "range")

    task = [dask.delayed(pyart.correct.cloud_threshold)(row) for row in result]
    noise = dask.compute(*task)
    noise = np.array(noise, dtype=float)

    assert_almost_equal(noise[0:5], expected, decimal=3)


def test_calc_cloud_mask():
    expected_cloud_1_first = [1, 1, 1, 1, 1]
    expected_cloud_1_last = [0, 1, 1, 1, 1]
    expected_cloud_2_first = [0, 0, 0, 0, 0]
    expected_cloud_2_last = [0, 0, 0, 1, 1]

    radar_mask = pyart.correct.calc_cloud_mask(
        radar_kazr, "reflectivity_copol", "range"
    )

    assert "cloud_mask_1" in radar_mask.fields.keys()
    assert "cloud_mask_2" in radar_mask.fields.keys()
    assert_array_equal(
        expected_cloud_1_first, radar_mask.fields["cloud_mask_1"]["data"][0][0:5]
    )
    assert_array_equal(
        expected_cloud_1_last, radar_mask.fields["cloud_mask_1"]["data"][-1][0:5]
    )
    assert_array_equal(
        expected_cloud_2_first, radar_mask.fields["cloud_mask_2"]["data"][0][0:5]
    )
    assert_array_equal(
        expected_cloud_2_last, radar_mask.fields["cloud_mask_2"]["data"][-1][0:5]
    )

    assert (
        radar_mask.fields["cloud_mask_1"]["long_name"]
        == "Cloud mask 1 (linear profile)"
    )
    assert radar_mask.fields["cloud_mask_2"]["long_name"] == "Cloud mask 2 (2D box)"

    assert radar_mask.fields["cloud_mask_2"]["units"] == "1"

    assert radar_mask.fields["cloud_mask_2"]["flag_meanings"] == ["no_cloud", "cloud"]

    bad_radar = "foo"
    pytest.raises(
        ValueError,
        pyart.correct.calc_cloud_mask,
        bad_radar,
        "reflectivity_copol",
        "range",
    )

    bad_field = 42
    pytest.raises(ValueError, pyart.correct.calc_cloud_mask, radar, bad_field, "range")
