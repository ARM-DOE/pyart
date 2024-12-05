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
    import cvxopt  # noqa

    cvxopt_available = True
except ImportError:
    cvxopt_available = False

try:
    import glpk  # noqa

    glpk_available = True
except ImportError:
    glpk_available = False


try:
    import cylp.cy  # noqa

    cylp_available = True
except ImportError:
    cylp_available = False

import pyart

PATH = os.path.dirname(__file__)
REFERENCE_RAYS_FILE = os.path.join(PATH, "reference_rays.npz")


@pytest.mark.skipif(not glpk_available, reason="GLPK is not installed.")
def test_phase_proc_lp_glpk():
    radar, phidp, kdp = perform_phase_processing()
    ref = np.load(REFERENCE_RAYS_FILE)
    assert _ratio(ref["reference_phidp"], phidp["data"]) <= 0.01
    assert _ratio(ref["reference_kdp"], kdp["data"]) <= 0.01
    assert (
        _ratio(
            ref["reference_unfolded_phidp"],
            radar.fields["unfolded_differential_phase"]["data"],
        )
        <= 0.01
    )


@pytest.mark.skipif(not cvxopt_available, reason="CVXOPT is not installed.")
def test_phase_proc_lp_cvxopt():
    from cvxopt import solvers

    solvers.options["LPX_K_MSGLEV"] = 0  # supress screen output
    radar, phidp, kdp = perform_phase_processing("cvxopt")
    ref = np.load(REFERENCE_RAYS_FILE)
    assert _ratio(ref["reference_phidp"], phidp["data"]) <= 0.01
    assert _ratio(ref["reference_kdp"], kdp["data"]) <= 0.01
    assert (
        _ratio(
            ref["reference_unfolded_phidp"],
            radar.fields["unfolded_differential_phase"]["data"],
        )
        <= 0.01
    )


@pytest.mark.skipif(not cylp_available, reason="CyLP is not installed.")
def test_phase_proc_lp_cylp():
    with warnings.catch_warnings():
        # ignore FutureWarnings as CyLP emits a number of these
        warnings.simplefilter("ignore", category=FutureWarning)
        radar, phidp, kdp = perform_phase_processing("cylp")
    ref = np.load(REFERENCE_RAYS_FILE)
    assert _ratio(ref["reference_phidp"], phidp["data"]) <= 0.01
    assert _ratio(ref["reference_kdp"], kdp["data"]) <= 0.01
    assert (
        _ratio(
            ref["reference_unfolded_phidp"],
            radar.fields["unfolded_differential_phase"]["data"],
        )
        <= 0.01
    )


@pytest.mark.skipif(not cylp_available, reason="CyLP is not installed.")
def test_phase_proc_lp_gf_cylp():
    with warnings.catch_warnings():
        # ignore FutureWarnings as CyLP emits a number of these
        warnings.simplefilter("ignore", category=FutureWarning)
        radar, phidp, kdp = test_perform_phase_processing_gf("cylp")
    ref = np.load(REFERENCE_RAYS_FILE)
    assert _ratio(ref["reference_phidp"], phidp["data"]) <= 0.05
    assert _ratio(ref["reference_kdp"], kdp["data"]) <= 0.05
    assert (
        _ratio(
            ref["reference_unfolded_phidp"],
            radar.fields["unfolded_differential_phase"]["data"],
        )
        <= 0.05
    )


@pytest.mark.skipif(not cylp_available, reason="CyLP is not installed.")
def test_phase_proc_lp_cylp_mp():
    with warnings.catch_warnings():
        # ignore FutureWarnings as CyLP emits a number of these
        warnings.simplefilter("ignore", category=FutureWarning)
        radar, phidp, kdp = perform_phase_processing("cylp_mp")
    ref = np.load(REFERENCE_RAYS_FILE)
    assert _ratio(ref["reference_phidp"], phidp["data"]) <= 0.01
    assert _ratio(ref["reference_kdp"], kdp["data"]) <= 0.01
    assert (
        _ratio(
            ref["reference_unfolded_phidp"],
            radar.fields["unfolded_differential_phase"]["data"],
        )
        <= 0.01
    )


def _ratio(a1, a2):
    """Ratio the sum of the abs difference vs sum abs of two vectors."""
    abs_residues = np.abs(a1 - a2).sum()
    avg_abs_sum = 0.5 * np.abs(a1).sum() + 0.5 * np.abs(a2).sum()
    return abs_residues / avg_abs_sum


def perform_phase_processing(LP_solver="cvxopt"):
    """Perform LP phase processing on a single ray radar."""
    radar = pyart.testing.make_single_ray_radar()
    phidp, kdp = pyart.correct.phase_proc_lp(radar, 0.0, LP_solver=LP_solver)
    return radar, phidp, kdp


@pytest.mark.skipif(not cvxopt_available, reason="CVXOPT is not installed.")
def test_perform_phase_processing_gf(LP_solver="cvxopt"):
    """Perform LP phase processing on a single ray radar."""
    radar = pyart.testing.make_single_ray_radar()
    my_gatefilter = pyart.filters.GateFilter(radar)
    my_gatefilter.exclude_below("normalized_coherent_power", 0.5)
    my_gatefilter.exclude_below("cross_correlation_ratio", 0.8)
    phidp, kdp = pyart.correct.phase_proc_lp_gf(
        radar,
        gatefilter=my_gatefilter,
        LP_solver=LP_solver,
        doc=15,
        ncpts=20,
        system_phase=-140.1,
    )
    assert phidp["data"].max() < 360.0
    assert kdp["data"].max() < 6.0
