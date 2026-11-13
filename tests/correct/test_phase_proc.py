"""Unit Tests for Py-ART's correct/proc_phase.py module."""

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

try:
    import highspy  # noqa

    highspy_available = True
except ImportError:
    highspy_available = False

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


@pytest.mark.skipif(not highspy_available, reason="HiGHS is not installed.")
def test_phase_proc_lp_highspy():
    radar, phidp, kdp = perform_phase_processing("highspy")
    ref = np.load(REFERENCE_RAYS_FILE)
    assert _ratio(ref["reference_phidp"], phidp["data"]) <= 0.011
    assert _ratio(ref["reference_kdp"], kdp["data"]) <= 0.01
    assert (
        _ratio(
            ref["reference_unfolded_phidp"],
            radar.fields["unfolded_differential_phase"]["data"],
        )
        <= 0.01
    )


@pytest.mark.skipif(not highspy_available, reason="HiGHS is not installed.")
def test_phase_proc_lp_highspy_max_kdp():
    _, _, baseline_kdp = perform_phase_processing("highspy")
    _, _, kdp = perform_phase_processing("highspy", max_kdp=5.0)
    assert kdp["valid_max"] == 5.0
    assert np.array_equal(
        np.ma.getmaskarray(kdp["data"]), np.ma.getmaskarray(baseline_kdp["data"])
    )
    assert kdp["data"].compressed().max() <= 5.0 + 1e-10


def _ratio(a1, a2):
    """Ratio the sum of the abs difference vs sum abs of two vectors."""
    abs_residues = np.abs(a1 - a2).sum()
    avg_abs_sum = 0.5 * np.abs(a1).sum() + 0.5 * np.abs(a2).sum()
    return abs_residues / avg_abs_sum


def perform_phase_processing(LP_solver="cvxopt", **kwargs):
    """Perform LP phase processing on a single ray radar."""
    radar = pyart.testing.make_single_ray_radar()
    phidp, kdp = pyart.correct.phase_proc_lp(radar, 0.0, LP_solver=LP_solver, **kwargs)
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


@pytest.mark.skipif(not highspy_available, reason="HiGHS is not installed.")
def test_phase_proc_lp_gf_highspy():
    radar = pyart.testing.make_single_ray_radar()
    my_gatefilter = pyart.filters.GateFilter(radar)
    my_gatefilter.exclude_below("normalized_coherent_power", 0.5)
    my_gatefilter.exclude_below("cross_correlation_ratio", 0.8)
    phidp, kdp = pyart.correct.phase_proc_lp_gf(
        radar,
        gatefilter=my_gatefilter,
        LP_solver="highspy",
        doc=15,
        ncpts=20,
        system_phase=-140.1,
    )
    assert phidp["data"].max() < 360.0
    assert kdp["data"].max() < 6.0


@pytest.mark.skipif(not highspy_available, reason="HiGHS is not installed.")
def test_phase_proc_lp_gf_highspy_max_kdp():
    radar = pyart.testing.make_single_ray_radar()
    my_gatefilter = pyart.filters.GateFilter(radar)
    my_gatefilter.exclude_below("normalized_coherent_power", 0.5)
    my_gatefilter.exclude_below("cross_correlation_ratio", 0.8)
    _, baseline_kdp = pyart.correct.phase_proc_lp_gf(
        radar,
        gatefilter=my_gatefilter,
        LP_solver="highspy",
        doc=15,
        ncpts=20,
        system_phase=-140.1,
    )
    _, kdp = pyart.correct.phase_proc_lp_gf(
        radar,
        gatefilter=my_gatefilter,
        LP_solver="highspy",
        doc=15,
        ncpts=20,
        system_phase=-140.1,
        max_kdp=5.0,
    )
    assert kdp["valid_max"] == 5.0
    assert np.array_equal(
        np.ma.getmaskarray(kdp["data"]),
        np.ma.getmaskarray(baseline_kdp["data"]),
    )
    assert kdp["data"].compressed().max() <= 5.0 + 1e-10


def _make_two_sweep_radar(ngates=100, rays_per_sweep=10):
    """Radar with PhiDP = 50 deg in sweep 0 and 80 deg in sweep 1."""
    radar = pyart.testing.make_empty_ppi_radar(ngates, rays_per_sweep, 2)
    shape = (2 * rays_per_sweep, ngates)

    phidp = np.full(shape, 50.0)
    phidp[rays_per_sweep:, :] = 80.0

    radar.add_field("differential_phase", {"data": phidp})
    radar.add_field("normalized_coherent_power", {"data": np.ones(shape)})
    radar.add_field("cross_correlation_ratio", {"data": np.ones(shape)})
    return radar


def test_det_sys_phase_default_is_first_sweep():
    radar = _make_two_sweep_radar()
    assert pyart.correct.phase_proc.det_sys_phase(radar) == pytest.approx(50.0)


def test_det_sys_phase_sweep_selection():
    radar = _make_two_sweep_radar()
    sys_phase = pyart.correct.phase_proc.det_sys_phase
    assert sys_phase(radar, sweep=0) == pytest.approx(50.0)
    assert sys_phase(radar, sweep=1) == pytest.approx(80.0)
    # negative index selects the last sweep (lowest sweep in top-down volumes)
    assert sys_phase(radar, sweep=-1) == pytest.approx(80.0)


def test_det_sys_phase_sweep_ignores_other_sweeps():
    """Invalid data in sweep 0 must not affect the estimate from sweep 1."""
    radar = _make_two_sweep_radar()
    radar.fields["normalized_coherent_power"]["data"][:10, :] = 0.0
    sys_phase = pyart.correct.phase_proc.det_sys_phase
    assert sys_phase(radar, sweep=0) is None
    assert sys_phase(radar, sweep=1) == pytest.approx(80.0)


def test_det_sys_phase_gf_sweep_selection():
    radar = _make_two_sweep_radar()
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_below("normalized_coherent_power", 0.5)
    gatefilter.exclude_below("cross_correlation_ratio", 0.8)
    sys_phase_gf = pyart.correct.phase_proc.det_sys_phase_gf
    assert sys_phase_gf(radar, gatefilter) == pytest.approx(50.0)
    assert sys_phase_gf(radar, gatefilter, sweep=0) == pytest.approx(50.0)
    assert sys_phase_gf(radar, gatefilter, sweep=1) == pytest.approx(80.0)
    assert sys_phase_gf(radar, gatefilter, sweep=-1) == pytest.approx(80.0)


def test_det_sys_phase_gf_sweep_ignores_other_sweeps():
    radar = _make_two_sweep_radar()
    radar.fields["normalized_coherent_power"]["data"][:10, :] = 0.0
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_below("normalized_coherent_power", 0.5)
    sys_phase_gf = pyart.correct.phase_proc.det_sys_phase_gf
    assert sys_phase_gf(radar, gatefilter, sweep=0) is None
    assert sys_phase_gf(radar, gatefilter, sweep=1) == pytest.approx(80.0)


def test_det_sys_phase_gf_first_gate_none_and_float():
    """first_gate=None and a float (the old default) must both work."""
    radar = _make_two_sweep_radar()
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.include_all()
    sys_phase_gf = pyart.correct.phase_proc.det_sys_phase_gf
    assert sys_phase_gf(radar, gatefilter, first_gate=None) == pytest.approx(50.0)
    assert sys_phase_gf(radar, gatefilter, first_gate=30.0) == pytest.approx(50.0)


def _make_rhi_radar(elevations, ngates=60, gate_spacing=250.0, altitude=0.0):
    """Synthetic RHI radar with one ray per requested elevation."""
    elevations = np.asarray(elevations, dtype=float)
    nrays = len(elevations)
    radar = pyart.testing.make_empty_rhi_radar(ngates, nrays, 1)
    radar.range["data"] = np.arange(ngates) * gate_spacing
    radar.elevation["data"] = elevations
    radar.altitude["data"] = np.array([altitude])
    # fixed_angle is the azimuth for an RHI; the factory default (0.75)
    # already matches radar.azimuth, so it is left alone.

    shape = (nrays, ngates)
    radar.add_field(
        "differential_phase",
        {"data": np.tile(np.linspace(0.0, 30.0, ngates), (nrays, 1))},
    )
    radar.add_field("reflectivity", {"data": np.full(shape, 30.0)})
    radar.add_field("normalized_coherent_power", {"data": np.ones(shape)})
    radar.add_field("cross_correlation_ratio", {"data": np.ones(shape)})
    return radar


def test_fzl_index_vertical_beam():
    """At 90 deg elevation gate height equals range."""
    ranges = np.arange(400) * 250.0
    # gates at 0, 250, ..., 3750 are below 4000 m; 4000 m is not
    assert pyart.correct.phase_proc.fzl_index(4000.0, ranges, 90.0, 0.0) == 15


def test_fzl_index_below_radar_returns_zero():
    """Freezing level under the radar: no gates qualify, return 0 not 6."""
    ranges = np.arange(40) * 250.0
    assert pyart.correct.phase_proc.fzl_index(4000.0, ranges, 5.0, 5000.0) == 0


def test_fzl_index_low_elevation_covers_all_gates():
    ranges = np.arange(40) * 250.0
    assert pyart.correct.phase_proc.fzl_index(4000.0, ranges, 0.0, 0.0) == 39


def test_det_process_range_rays_rhi_varies_with_elevation():
    radar = _make_rhi_radar([0.0, 30.0, 90.0], ngates=400)
    gate_ends, start, end = pyart.correct.phase_proc.det_process_range_rays(
        radar, 0, 4000.0, doc=0
    )
    assert (start, end) == (0, 3)
    assert gate_ends[0] == 399
    assert gate_ends[2] == 15
    assert np.all(np.diff(gate_ends) < 0)


def test_det_process_range_rays_respects_doc():
    radar = _make_rhi_radar([0.0, 90.0], ngates=400)
    gate_ends, _, _ = pyart.correct.phase_proc.det_process_range_rays(
        radar, 0, 4000.0, doc=10
    )
    assert gate_ends[0] == 390
    assert gate_ends[1] == 15


def test_det_process_range_rhi_uses_minimum_over_rays():
    radar = _make_rhi_radar([0.0, 30.0, 90.0], ngates=400)
    gate_end, start, end = pyart.correct.phase_proc.det_process_range(
        radar, 0, 4000.0, doc=0
    )
    assert gate_end == 15
    assert (start, end) == (0, 3)


def test_det_process_range_rhi_ignores_fixed_angle():
    """RHI gate_end must come from per-ray elevations, not fixed_angle."""
    radar = _make_rhi_radar([0.0, 0.0, 0.0], ngates=400)
    radar.fixed_angle["data"] = np.array([90.0])  # misleading on purpose
    gate_end, _, _ = pyart.correct.phase_proc.det_process_range(radar, 0, 4000.0, doc=0)
    assert gate_end == 399


def test_det_process_range_ppi_uses_fixed_angle():
    radar = pyart.testing.make_empty_ppi_radar(400, 5, 1)
    radar.range["data"] = np.arange(400) * 250.0
    radar.altitude["data"] = np.array([0.0])
    radar.fixed_angle["data"] = np.array([90.0])
    gate_end, start, end = pyart.correct.phase_proc.det_process_range(
        radar, 0, 4000.0, doc=0
    )
    assert gate_end == 15
    assert (start, end) == (0, 5)


def test_unfold_ray_all_excluded_returns_none():
    phidp = np.linspace(0.0, 50.0, 100)
    notmeteo = np.ones(100, dtype=bool)
    out = pyart.correct.phase_proc._unfold_ray(phidp, notmeteo, 2, None, 0.0)
    assert out is None


def test_unfold_ray_valid_ray():
    phidp = np.linspace(0.0, 50.0, 100)
    notmeteo = np.zeros(100, dtype=bool)
    out = pyart.correct.phase_proc._unfold_ray(phidp, notmeteo, 2, None, 0.0)
    assert out is not None
    assert out.shape == (100,)
    assert np.all(np.isfinite(out))
    assert out[0] == pytest.approx(0.0)


def test_unfold_ray_does_not_modify_notmeteo():
    phidp = np.linspace(0.0, 50.0, 100)
    notmeteo = np.zeros(100, dtype=bool)
    notmeteo[:5] = True
    notmeteo[-1] = True
    before = notmeteo.copy()
    pyart.correct.phase_proc._unfold_ray(phidp, notmeteo, 2, None, 0.0)
    assert np.array_equal(notmeteo, before)


def test_get_phidp_unf_gf_does_not_modify_gatefilter():
    radar = _make_two_sweep_radar()
    radar.fields["normalized_coherent_power"]["data"][:, 0] = 0.0
    radar.fields["normalized_coherent_power"]["data"][:, -1] = 0.0
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_below("normalized_coherent_power", 0.5)
    before = gatefilter.gate_excluded.copy()
    assert before[:, 0].all() and before[:, -1].all()  # filter is non-trivial

    pyart.correct.phase_proc.get_phidp_unf_gf(radar, gatefilter, sys_phase=0.0)
    assert np.array_equal(before, gatefilter.gate_excluded)


def test_get_phidp_unf_gf_fully_filtered_gives_zeros():
    radar = _make_two_sweep_radar()
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_all()
    out = pyart.correct.phase_proc.get_phidp_unf_gf(radar, gatefilter, sys_phase=0.0)
    assert np.all(np.isfinite(out))
    assert np.all(out == 0.0)


def test_get_phidp_unf_gf_ncpts_none_matches_two():
    """ncpts=None must fall back to 2 rather than zeroing every ray."""
    radar = _make_two_sweep_radar()
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.include_all()
    unf = pyart.correct.phase_proc.get_phidp_unf_gf
    out_none = unf(radar, gatefilter, ncpts=None, sys_phase=0.0)
    out_two = unf(radar, gatefilter, ncpts=2, sys_phase=0.0)
    assert np.allclose(out_none, out_two)
    assert np.any(out_none != 0.0)


@pytest.mark.skipif(not cvxopt_available, reason="CVXOPT is not installed.")
def test_phase_proc_lp_gf_rhi_masks_gates_above_fzl():
    from cvxopt import solvers

    solvers.options["show_progress"] = False
    # 250 m gates, fzl = 2000 m: ray 0 (0 deg) keeps almost all gates,
    # ray 1 (90 deg) is cut after gate 7
    radar = _make_rhi_radar([0.0, 90.0], ngates=60)
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.include_all()

    phidp, kdp = pyart.correct.phase_proc_lp_gf(
        radar,
        gatefilter=gatefilter,
        LP_solver="cvxopt",
        fzl=2000.0,
        doc=0,
        ncpts=2,
        system_phase=0.0,
    )

    assert np.ma.isMaskedArray(kdp["data"])
    assert np.ma.isMaskedArray(phidp["data"])
    kdp_mask = np.ma.getmaskarray(kdp["data"])
    # 90 deg ray: everything past the cutoff is masked, the start is not
    assert kdp_mask[1, 8:].all()
    assert not kdp_mask[1, :7].any()
    # 0 deg ray: processed over (almost) the whole ray
    assert not kdp_mask[0, :50].any()
    assert np.isfinite(kdp["data"].compressed()).all()


@pytest.mark.skipif(not cvxopt_available, reason="CVXOPT is not installed.")
def test_phase_proc_lp_gf_rhi_below_fzl_is_masked_not_error():
    """Radar above the freezing level: nothing processed, no exception."""
    from cvxopt import solvers

    solvers.options["show_progress"] = False
    radar = _make_rhi_radar([10.0, 20.0], ngates=60, altitude=6000.0)
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.include_all()

    phidp, kdp = pyart.correct.phase_proc_lp_gf(
        radar,
        gatefilter=gatefilter,
        LP_solver="cvxopt",
        fzl=2000.0,
        doc=0,
        ncpts=2,
        system_phase=0.0,
    )
    assert np.ma.getmaskarray(kdp["data"]).all()


@pytest.mark.skipif(not cvxopt_available, reason="CVXOPT is not installed.")
def test_phase_proc_lp_gf_does_not_modify_gatefilter():
    radar = pyart.testing.make_single_ray_radar()
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_below("normalized_coherent_power", 0.5)
    gatefilter.exclude_below("cross_correlation_ratio", 0.8)
    before = gatefilter.gate_excluded.copy()
    pyart.correct.phase_proc_lp_gf(
        radar,
        gatefilter=gatefilter,
        LP_solver="cvxopt",
        doc=15,
        ncpts=20,
        system_phase=-140.1,
    )
    assert np.array_equal(before, gatefilter.gate_excluded)
