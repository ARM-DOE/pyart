"""
A script for comparing the e Vulpiani, Maesaka and Schneedbeli methods
for KDP estimation with sample data.
"""

import numpy as np
import matplotlib.pyplot as plt
import pytest

from pyart.retrieve.kdp_proc import kdp_maesaka, kdp_vulpiani, kdp_schneebeli
from pyart.testing import sample_objects
from pyart.config import get_field_name


def _make_real_psidp_radar():
    """
    Create single-ray radar with linear differential phase profile with
    specified slope.
    ---
    Returns
    -------
    radar : Radar
            PyART radar instance with differential phase profile in deg.
    """
    psidp = np.array([[-2.33313751e+00, 1.80617523e+00, 7.17742920e-01,
                     1.82811661e+01, 1.89352417e+01, 1.67904205e+01]])
    psidp = np.ma.array(psidp)
    radar = sample_objects.make_empty_ppi_radar(len(psidp[0]), 1, 1)
    psidp_dict = {
                  'data': psidp,
                  }
    radar.add_field(get_field_name('differential_phase'), psidp_dict)
    # Define real ranges
    radar.range['data'] = 75*np.arange(0, len(psidp[0]))
    return radar


@pytest.mark.mpl_image_compare(tolerance=30)
def test_compare_kdp_estimation_methods():
    # Get profile of noisy psidp
    prof_psidp = _make_real_psidp_radar()
    # Maesaka method
    kdp_mae, phidpf_mae, phidp_mae = kdp_maesaka(prof_psidp,
                                                 maxiter=1000,
                                                 check_outliers=False,
                                                 parallel=False)

    # Vulpiani method (note windsize is just a guess here..)
    kdp_vulp, phidp_vulp = kdp_vulpiani(prof_psidp,
                                        windsize=2,
                                        n_iter=20,
                                        parallel=False,
                                        band='X')
    # Kalman filter method
    kdp_schnee, kdp_std_schnee, phidp_schnee = kdp_schneebeli(prof_psidp,
                                                              parallel=False,
                                                              band='X')
    # Create figure
    fig = plt.figure(figsize=(10, 10))
    plt.subplot(2, 1, 1)
    plt.grid(True)
    plt.title('Kdp estimation')
    ranges = prof_psidp.range['data']
    plt.plot(ranges, kdp_mae['data'][0])
    plt.plot(ranges, kdp_vulp['data'][0])
    plt.plot(ranges, kdp_schnee['data'][0])
    plt.xlabel('Range [m]')
    plt.ylabel('Kdp [deg/km]')
    plt.legend(['Maesaka', 'Vulpiani', 'Schneebeli'], loc=0)
    plt.subplot(2, 1, 2)
    plt.grid(True)
    plt.title('Reconstructed Phidp')
    ranges = prof_psidp.range['data']
    phidp_mae = 0.5 * (phidp_mae['data'][0] + phidp_mae['data'][0])
    plt.plot(ranges, phidp_mae)
    plt.plot(ranges, phidp_vulp['data'][0])
    plt.plot(ranges, phidp_schnee['data'][0])
    plt.plot(ranges, prof_psidp.fields['differential_phase']['data'][0])
    plt.xlabel('Range [m]')
    plt.ylabel('Diff. phase [deg]')
    plt.legend(['Maesaka', 'Vulpiani', 'Schneebeli', 'Real Psidp'], loc=0)
    try:
        return fig
    finally:
        plt.close(fig)
