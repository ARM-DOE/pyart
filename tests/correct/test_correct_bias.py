""" Unit tests for bias and noise module. """

import pytest
import pyart
from open_radar_data import DATASETS
import numpy as np

from numpy.testing import assert_allclose
radar = pyart.io.read(pyart.testing.NEXRAD_ARCHIVE_MSG31_FILE)


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
    xsapr_test_file = DATASETS.fetch('sgpxsaprcfrvptI4.a1.20200205.100827.nc')
    ds = pyart.io.read(xsapr_test_file)
    gatefilter = pyart.filters.GateFilter(ds)
    gatefilter.exclude_below('cross_correlation_ratio_hv', 0.995)
    gatefilter.exclude_above('cross_correlation_ratio_hv', 1)
    gatefilter.exclude_below('reflectivity', 10)
    gatefilter.exclude_above('reflectivity', 30)

    results = pyart.correct.calc_zdr_offset(ds, zdr_var='differential_reflectivity', gatefilter=gatefilter,
                                        height_range=(1000, 3000))
    np.testing.assert_almost_equal(results['bias'], 2.69, decimal=2)
    print(results['profile_reflectivity'][1:16])
    np.testing.assert_almost_equal(results['profile_reflectivity'][15], 14.37, decimal=2)
