""" Unit tests for pyart.retrieve.kdp_proc module. """

import numpy as np

from pyart.retrieve import kdp_proc
from pyart.filters import GateFilter
from pyart.testing import sample_objects
from pyart.config import get_field_name

# Should include a test to make sure KDP -> slope of a linear differential
# phase profile.


def test_kdp_maesaka_linear_psidp(maxiter=100):
    """ Test KDP -> constant for linear PSIDP profile. """
    radar = _make_linear_psidp_radar()
    kdp_dict, phidpf_dict, phidpr_dict = kdp_proc.kdp_maesaka(
        radar, maxiter=maxiter, check_outliers=False)

    assert np.allclose(np.diff(kdp_dict['data'][0]), 0.0, atol=0.1)

    return


def test_kdp_maesaka_all_excluded(first_guess=0.001, maxiter=100):
    """ Test for first guess field if all gates are excluded. """
    radar = _make_linear_psidp_radar()
    gatefilter = GateFilter(radar)
    gatefilter.exclude_all()
    kdp_dict, phidpf_dict, phidpr_dict = kdp_proc.kdp_maesaka(
        radar, gatefilter=gatefilter, first_guess=first_guess, maxiter=maxiter,
        check_outliers=False)

    assert np.allclose(kdp_dict['data'][0], 0.0, atol=first_guess)

    return

def _make_linear_psidp_radar(slope=None):
    """ Create single-ray radar with linear differential phase profile. """
    radar = sample_objects.make_empty_ppi_radar(101, 1, 1)
    psidp_dict = {
        'data': np.atleast_2d(np.linspace(0.0, 100.0, radar.ngates))
        }
    radar.add_field(get_field_name('differential_phase'), psidp_dict)

    return radar


