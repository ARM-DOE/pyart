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
        [71.34390136, 71.34390136, 71.34390136, 71.34390136, 71.34390136],
        atol=1e-14)
    assert_allclose(
        fields['velocity']['data'][0, 0:5],
        [2.73382043e-16, 2.73382043e-16, 2.73382043e-16, 2.73382043e-16,
         2.73382043e-16],
        atol=1e-14)
    assert_allclose(
        fields['spectrum_width']['data'][0, 0:5],
        [2.85047687, 2.85047687, 2.85047687, 2.85047687, 2.85047687],
        atol=1e-14)
    assert_allclose(
        fields['skewness']['data'][0, 0:5],
        [9.44294386e-17, 9.44294386e-17, 9.44294386e-17, 9.44294386e-17,
         9.44294386e-17],
        atol=1e-14)
    assert_allclose(
        fields['kurtosis']['data'][0, 0:5],
        [2.79911197, 2.79911197, 2.79911197, 2.79911197, 2.79911197],
        atol=1e-14)
