""" Unit tests for pyart's retrieve/kdp_proc.py module. """

import numpy as np

from pyart.retrieve import kdp_proc
from pyart.filters import GateFilter
from pyart.testing import sample_objects
from pyart.config import get_field_name


def test_kdp_maesaka_linear_psidp(slope=0.002, maxiter=100):
    radar = _make_linear_psidp_radar(slope=slope)
    kdp_dict, phidpf_dict, phidpr_dict = kdp_proc.kdp_maesaka(
        radar, maxiter=maxiter, check_outliers=False)

    assert np.allclose(np.diff(kdp_dict['data'][0]), 0.0, atol=0.1)
    assert np.allclose(kdp_dict['data'], 1000.0 * slope / 2.0, atol=0.1)

    return


def test_kdp_maesaka_all_excluded(first_guess=0.01, maxiter=100):
    radar = _make_linear_psidp_radar()
    gatefilter = GateFilter(radar)
    gatefilter.exclude_all()
    kdp_dict, phidpf_dict, phidpr_dict = kdp_proc.kdp_maesaka(
        radar, gatefilter=gatefilter, first_guess=first_guess, maxiter=maxiter,
        check_outliers=False)

    assert np.allclose(kdp_dict['data'][0], 0.0, atol=first_guess)

    return


def _make_linear_psidp_radar(slope=0.002):
    """
    Create single-ray radar with linear differential phase profile with
    specified slope.

    Parameters
    ----------
    slope : float, optional
        Slope of differential phase profile in deg/m. Radar range gates cover
        0-1000 m, inclusive, with 10 m gate spacings.

    Returns
    -------
    radar : Radar
        Radar with linear differential phase profile in deg.

    """
    radar = sample_objects.make_empty_ppi_radar(101, 1, 1)
    psidp_dict = {
        'data': np.atleast_2d(np.linspace(0.0, slope * 1000.0, radar.ngates))
        }
    radar.add_field(get_field_name('differential_phase'), psidp_dict)

    return radar
