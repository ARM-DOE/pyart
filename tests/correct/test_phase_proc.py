""" Unit Tests for Py-ART's correct/proc_phase.py module. """

# Can also be run as a script to create a ray_plot.png file
# python test_phase_proc.py
# This file should look like the reference_ray_plot.png file

# Use:
# python test_phase_proc.py -r
# to recreate the reference_rays.npz and reference_ray_plot.png files

import os
import warnings

import numpy as np
import pytest

try:
    import cvxopt
    cvxopt_available = True
except ImportError:
    cvxopt_available = False

try:
    import glpk
    glpk_available = True
except ImportError:
    glpk_available = False


try:
    import cylp.cy
    cylp_available = True
except ImportError:
    cylp_available = False

import pyart


PATH = os.path.dirname(__file__)
REFERENCE_RAYS_FILE = os.path.join(PATH, 'reference_rays.npz')


@pytest.mark.skipif(not glpk_available,
                    reason="GLPK is not installed.")
def test_phase_proc_lp_glpk():
    radar, phidp, kdp = perform_phase_processing()
    ref = np.load(REFERENCE_RAYS_FILE)
    assert _ratio(ref['reference_phidp'], phidp['data']) <= 0.01
    assert _ratio(ref['reference_kdp'], kdp['data']) <= 0.01
    assert _ratio(ref['reference_unfolded_phidp'],
                  radar.fields['unfolded_differential_phase']['data']) <= 0.01


@pytest.mark.skipif(not cvxopt_available,
                    reason="CVXOPT is not installed.")
def test_phase_proc_lp_cvxopt():
    from cvxopt import solvers
    solvers.options['LPX_K_MSGLEV'] = 0     # supress screen output
    radar, phidp, kdp = perform_phase_processing('cvxopt')
    ref = np.load(REFERENCE_RAYS_FILE)
    assert _ratio(ref['reference_phidp'], phidp['data']) <= 0.01
    assert _ratio(ref['reference_kdp'], kdp['data']) <= 0.01
    assert _ratio(ref['reference_unfolded_phidp'],
                  radar.fields['unfolded_differential_phase']['data']) <= 0.01


@pytest.mark.skipif(not cylp_available,
                    reason="CyLP is not installed.")
def test_phase_proc_lp_cylp():
    with warnings.catch_warnings():
        # ignore FutureWarnings as CyLP emits a number of these
        warnings.simplefilter("ignore", category=FutureWarning)
        radar, phidp, kdp = perform_phase_processing('cylp')
    ref = np.load(REFERENCE_RAYS_FILE)
    assert _ratio(ref['reference_phidp'], phidp['data']) <= 0.01
    assert _ratio(ref['reference_kdp'], kdp['data']) <= 0.01
    assert _ratio(ref['reference_unfolded_phidp'],
                  radar.fields['unfolded_differential_phase']['data']) <= 0.01

@pytest.mark.skipif(not cylp_available,
                    reason="CyLP is not installed.")
def test_phase_proc_lp_gf_cylp():
    with warnings.catch_warnings():
        # ignore FutureWarnings as CyLP emits a number of these
        warnings.simplefilter("ignore", category=FutureWarning)
        radar, phidp, kdp = perform_phase_processing_gf('cylp')
    ref = np.load(REFERENCE_RAYS_FILE)
    assert _ratio(ref['reference_phidp'], phidp['data']) <= 0.05
    assert _ratio(ref['reference_kdp'], kdp['data']) <= 0.05
    assert _ratio(ref['reference_unfolded_phidp'],
                  radar.fields['unfolded_differential_phase']['data']) <= 0.05


@pytest.mark.skipif(not cylp_available,
                    reason="CyLP is not installed.")
def test_phase_proc_lp_cylp_mp():
    with warnings.catch_warnings():
        # ignore FutureWarnings as CyLP emits a number of these
        warnings.simplefilter("ignore", category=FutureWarning)
        radar, phidp, kdp = perform_phase_processing('cylp_mp')
    ref = np.load(REFERENCE_RAYS_FILE)
    assert _ratio(ref['reference_phidp'], phidp['data']) <= 0.01
    assert _ratio(ref['reference_kdp'], kdp['data']) <= 0.01
    assert _ratio(ref['reference_unfolded_phidp'],
                  radar.fields['unfolded_differential_phase']['data']) <= 0.01


def _ratio(a1, a2):
    """ Ratio the sum of the abs difference vs sum abs of two vectors. """
    abs_residues = np.abs(a1 - a2).sum()
    avg_abs_sum = 0.5 * np.abs(a1).sum() + 0.5 * np.abs(a2).sum()
    return abs_residues / avg_abs_sum


def perform_phase_processing(LP_solver='pyglpk'):
    """ Perform LP phase processing on a single ray radar. """
    radar = pyart.testing.make_single_ray_radar()
    phidp, kdp = pyart.correct.phase_proc_lp(radar, 0.0, LP_solver=LP_solver)
    return radar, phidp, kdp

def perform_phase_processing_gf(LP_solver='pyglpk'):
    """ Perform LP phase processing on a single ray radar. """
    radar = pyart.testing.make_single_ray_radar()
    my_gatefilter = pyart.filters.GateFilter(radar)
    my_gatefilter.exclude_below('normalized_coherent_power', 0.5)
    my_gatefilter.exclude_below('cross_correlation_ratio', 0.8)
    phidp, kdp = pyart.correct.phase_proc_lp_gf(radar, gatefilter=my_gatefilter,
                                                LP_solver=LP_solver, doc=15,
                                                ncpts=20, system_phase=-140.1)
    return radar, phidp, kdp


def save_reference_rays(radar, phidp, kdp):
    """ Save the phase processed rays to REFERENCE_RAY_FILE. """
    unfolded = radar.fields['unfolded_differential_phase']
    np.savez(
        REFERENCE_RAYS_FILE,
        reference_phidp=phidp['data'],
        reference_kdp=kdp['data'],
        reference_unfolded_phidp=unfolded['data'])


def make_plot(range_km, unfolded_phidp, refl, phidp, kdp, filename):
    """
    Make and save a plot of PhiDP, reflectivity and KDP.

    Parameters
    ----------
    range_km : array
        Gate range in kilometers.
    unfolded_phidp, refl, phidp, kdp : dict
        Dictionary containing the unfolded phiDP, reflectivity,
        LP filtered phiDP and KDP in the data key
    filename :
        Filename to save the file to.

    """

    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=[10, 5])
    ax = fig.add_subplot(111)

    # filtered phidp and unfolded phidp
    p1, = ax.plot(range_km, phidp['data'][0], 'b-')
    p2, = ax.plot(range_km, unfolded_phidp['data'][0], 'g-')

    # set labels
    ax.set_ylim(0, 250)
    ax.set_ylabel('Differential phase shift (degrees)')
    ax.set_xlabel('Range (km)')

    # plot KDP and reflectivity on second axis
    ax2 = ax.twinx()
    p3, = ax2.plot(range_km, kdp['data'][0], 'r-')
    p4, = ax2.plot(range_km, refl['data'][0]/10.)

    # decorate and save
    ax2.yaxis.grid(color='gray', linestyle='dashed')
    ax.legend([p1, p2, p3, p4],
              ["Filtered phiDP", "Unfolded phiDP", 'KDP', 'Z/10.0'],
              loc='upper left')
    fig.savefig(filename)


if __name__ == "__main__":
    import sys

    radar, phidp, kdp = perform_phase_processing()
    filename = 'ray_plot.png'
    if sys.argv[-1] == '-r':    # regenerate reference files
        save_reference_rays(radar, phidp, kdp)
        filename = 'reference_ray_plot.png'

    make_plot(radar.range['data'] / 1000.0,
              radar.fields['unfolded_differential_phase'],
              radar.fields['reflectivity'],
              phidp, kdp, filename)
