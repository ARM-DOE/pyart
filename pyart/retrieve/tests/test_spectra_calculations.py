""" Unit tests for spectra_calculations.py """

from numpy.testing import assert_allclose
import pytest

from pyart.retrieve import spectra_moments
from pyart.testing import make_target_spectra_radar

try:
    import xarray as xr
    _XARRAY_AVAILABLE = True
except ImportError:
    _XARRAY_AVAILABLE = False

@pytest.mark.skipif(not _XARRAY_AVAILABLE,
                    reason="Xarray is not installed.")
def test_spectra_moments():
    radar = make_target_spectra_radar()
    fields = spectra_moments(radar)
    assert len(fields.keys()) == 5
    assert_allclose(
        fields['reflectivity']['data'][0, 0:5],
        [69.04937465, 69.04937465, 69.04937465, 69.04937465, 69.04937465],
        atol=1e-14)
    assert_allclose(
        fields['velocity']['data'][0, 0:5],
        [-7.47376351, -7.47376351, -7.47376351, -7.47376351, -7.47376351],
        atol=1e-14)
    assert_allclose(
        fields['spectrum_width']['data'][0, 0:5],
        [1.83767467, 1.83767467, 1.83767467, 1.83767467, 1.83767467],
        atol=1e-14)
    assert_allclose(
        fields['skewness']['data'][0, 0:5],
        [1.59127727, 1.59127727, 1.59127727, 1.59127727, 1.59127727],
        atol=1e-14)
    assert_allclose(
        fields['kurtosis']['data'][0, 0:5],
        [7.26160787, 7.26160787, 7.26160787, 7.26160787, 7.26160787],
        atol=1e-14)
