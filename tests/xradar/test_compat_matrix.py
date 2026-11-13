"""
Radar-vs-Xradar compatibility test matrix.

This module is the acceptance gate for the xradar-compatibility program: it
asserts that a fixed set of Py-ART algorithms produce numerically equal
results whether handed a classic ``pyart.core.Radar`` or the equivalent
``pyart.xradar.Xradar`` wrapper (and, for gridded data, ``pyart.core.Grid``
vs ``pyart.xradar.Xgrid``). A regression in the accessor should be caught
here rather than in a downstream bug report.

Fixture design
---------------
The brief calls for opening one real ``open_radar_data`` file two ways (as
in ``tests/xradar/test_accessor.py``) and comparing directly. In practice
this does not hold up for exact, per-gate numerical comparison: for every
real cfradial1 file inspected (single- and multi-sweep), ``pyart.io.read``
and ``xradar.io.open_cfradial1_datatree`` disagree on ray counts per sweep
and/or on where a sweep starts (e.g. which ray straddling the 0/360 azimuth
wrap belongs to which sweep). This is a property of the two independent
*readers* disagreeing on ray bookkeeping for messy, real-world antenna
dwell data -- not a defect in the ``Xradar``/``Xgrid`` accessor classes
under test here, and not something a field-name normalization can paper
over.

To get an apples-to-apples radar pair, ``radar_pair`` instead builds one
synthetic multi-sweep ``pyart.core.Radar`` (clean, non-duplicated azimuths)
with ``reflectivity``, ``velocity`` and ``differential_reflectivity``
fields, writes it once to a temporary CfRadial file with
``pyart.io.write_cfradial``, and opens *that* file both directly
(``pyart.io.read``) and through ``xradar.io.open_cfradial1_datatree`` ->
``Xradar``. This is the same "round trip through a temp cfradial file"
technique the brief itself prescribes for the RHI case, applied to the main
pair as well, and it eliminates ray-bookkeeping disagreements entirely
(verified: azimuths, slices and field values line up ray-for-ray). Field
names do not differ between the two read paths for this fixture (both see
the names given to the synthetic Radar), so no renaming is required --
noted here since the brief anticipates a rename step.

A second field, ``reflectivity_masked``, carries real missing gates (a
``numpy.ma.MaskedArray`` with a masked edge) and exists solely to exercise
``pyart.correct.despeckle_field``, which needs genuinely masked input on
the ``Radar`` side to avoid a numpy ``nomask`` collapse (a pre-existing
Py-ART despeckle quirk, confirmed to reproduce on a vanilla ``Radar`` with
an *unmasked* field too, independent of xradar). That case surfaces a real
``Xradar`` accessor defect: see ``test_despeckle_field`` below.

``grid_pair`` builds a small, fully-covered grid (no out-of-range gaps) so
that Grid/Xgrid comparisons are not tripped up by a separate, unrelated
pre-existing gap: ``pyart.io.write_grid`` does not write a netCDF
``_FillValue`` attribute for masked grid cells, so any masked value
surviving a Grid -> file -> Xgrid round trip reads back as the raw
``9.9692e+36`` CDF fill sentinel rather than NaN. That is a
``pyart.io.write_grid`` limitation predating xradar integration, out of
scope for this task; the fixture simply avoids masked cells so it does not
interfere with the accessor comparisons under test.
"""

import numpy as np
import pytest
import xarray as xr
import xradar as xd
from numpy.testing import assert_allclose, assert_array_equal

import pyart

NSWEEPS = 3
RAYS_PER_SWEEP = 60
NGATES = 20
NRAYS = NSWEEPS * RAYS_PER_SWEEP


def _build_synthetic_radar():
    """Build a small, clean multi-sweep PPI radar with several fields."""
    radar = pyart.testing.make_empty_ppi_radar(NGATES, RAYS_PER_SWEEP, NSWEEPS)
    radar.range["data"] = np.linspace(1000, 20000, NGATES).astype("float32")

    rng = np.random.default_rng(20260722)
    ramp = np.linspace(0, 6 * np.pi, NRAYS)[:, None]
    reflectivity = (20 + 10 * np.sin(ramp) + rng.normal(0, 1, (NRAYS, NGATES))).astype(
        "float32"
    )
    velocity = rng.normal(0, 5, (NRAYS, NGATES)).astype("float32")
    zdr = rng.normal(0.5, 0.2, (NRAYS, NGATES)).astype("float32")

    # A masked variant of reflectivity (masked far gates) used only for the
    # despeckle_field case, which needs real masked input on the Radar side.
    reflectivity_masked = np.ma.masked_array(reflectivity.copy(), mask=False)
    reflectivity_masked.mask[:, -3:] = True

    radar.fields = {
        "reflectivity": pyart.config.get_metadata("reflectivity"),
        "reflectivity_masked": pyart.config.get_metadata("reflectivity"),
        "velocity": pyart.config.get_metadata("velocity"),
        "differential_reflectivity": pyart.config.get_metadata(
            "differential_reflectivity"
        ),
    }
    radar.fields["reflectivity"]["data"] = reflectivity
    radar.fields["reflectivity_masked"]["data"] = reflectivity_masked
    radar.fields["velocity"]["data"] = velocity
    radar.fields["differential_reflectivity"]["data"] = zdr
    # Real Py-ART readers (nexrad, sigmet, uf, mdv, ...) set '_FillValue' on
    # every field at ingest time; this synthetic fixture does the same so
    # that write_cfradial() records it as a variable attribute and xarray's
    # CF decoding (used when reading the file back as an Xradar object) can
    # recover the mask, matching real-world round trips.
    for field in radar.fields.values():
        field["_FillValue"] = pyart.config.get_fillvalue()
    return radar


@pytest.fixture(scope="module")
def radar_pair(tmp_path_factory):
    """(Radar, Xradar) pair built from one synthetic source file.

    See the module docstring for why a synthetic + round-trip source is
    used instead of a real ``open_radar_data`` file.
    """
    radar = _build_synthetic_radar()
    tmp_dir = tmp_path_factory.mktemp("compat_matrix_ppi")
    filename = str(tmp_dir / "synthetic_ppi.nc")
    pyart.io.write_cfradial(filename, radar)

    radar_from_file = pyart.io.read(filename)
    dtree = xd.io.open_cfradial1_datatree(filename, optional=False)
    xradar_obj = pyart.xradar.Xradar(dtree)
    return radar_from_file, xradar_obj


@pytest.fixture(scope="module")
def grid_pair(radar_pair):
    """(Grid, Grid) pair, one from Radar, one from Xradar, over one field.

    A large ``roi_func`` radius is used so every grid point falls inside
    some gate's region of influence -- see module docstring for why gaps
    are avoided here.
    """
    radar, xradar_obj = radar_pair
    common_kwargs = dict(
        grid_shape=(3, 10, 10),
        grid_limits=((0, 3000), (-3000, 3000), (-3000, 3000)),
        fields=["reflectivity"],
        gridding_algo="map_to_grid",
        roi_func=lambda z, y, x: 1.0e6,
    )
    grid_from_radar = pyart.map.grid_from_radars((radar,), **common_kwargs)
    grid_from_xradar = pyart.map.grid_from_radars((xradar_obj,), **common_kwargs)
    return grid_from_radar, grid_from_xradar


@pytest.fixture(scope="module")
def xgrid_pair(grid_pair, tmp_path_factory):
    """(Grid, Xgrid) pair via the write_grid/read round trip used elsewhere
    in tests/xradar/test_accessor.py (test_grid_write_read, test_xgrid_*)."""
    grid, _ = grid_pair
    tmp_dir = tmp_path_factory.mktemp("compat_matrix_grid")
    filename = str(tmp_dir / "synthetic_grid.nc")
    pyart.io.write_grid(filename, grid)
    grid_ds = xr.open_dataset(filename, decode_times=False)
    xgrid = pyart.xradar.Xgrid(grid_ds)
    yield grid, xgrid
    xgrid.ds.close()


@pytest.fixture(scope="module")
def rhi_pair(tmp_path_factory):
    """(Radar, Xradar) pair for an RHI scan, via a tempfile cfradial round
    trip (there is no direct in-memory Radar -> DataTree conversion)."""
    radar = pyart.testing.make_empty_rhi_radar(15, 40, 2)
    nrays = 40 * 2
    ramp = np.linspace(0, 4 * np.pi, nrays)[:, None]
    data = (20 + 5 * np.sin(ramp) * np.ones((nrays, 15))).astype("float32")
    # Mask the far gates on every ray (same technique as radar_pair's
    # reflectivity_masked) so despeckle_field has real missing gates to
    # exclude on both the Radar and Xradar side -- an all-valid field would
    # leave the gatefilter's excluded-gate mask trivially all False, which
    # np.ma represents as a bare scalar (no shape) rather than a full array,
    # a separate pre-existing GateFilter.exclude_gates() shape-check quirk
    # unrelated to the Xradar accessor defects under test here.
    data = np.ma.masked_array(data, mask=False)
    data.mask[:, -3:] = True
    radar.fields = {"reflectivity": pyart.config.get_metadata("reflectivity")}
    radar.fields["reflectivity"]["data"] = data
    radar.fields["reflectivity"]["_FillValue"] = pyart.config.get_fillvalue()

    tmp_dir = tmp_path_factory.mktemp("compat_matrix_rhi")
    filename = str(tmp_dir / "synthetic_rhi.nc")
    pyart.io.write_cfradial(filename, radar)

    radar_from_file = pyart.io.read(filename)
    dtree = xd.io.open_cfradial1_datatree(filename, optional=False)
    xradar_obj = pyart.xradar.Xradar(dtree, scan_type="rhi")
    return radar_from_file, xradar_obj


# ---------------------------------------------------------------------------
# 1. calculate_velocity_texture
# ---------------------------------------------------------------------------
def test_calculate_velocity_texture(radar_pair):
    radar, xradar_obj = radar_pair
    kwargs = dict(vel_field="velocity", nyq=10.0, check_nyq_uniform=False)
    texture_radar = pyart.retrieve.calculate_velocity_texture(radar, **kwargs)
    texture_xradar = pyart.retrieve.calculate_velocity_texture(xradar_obj, **kwargs)
    assert_allclose(texture_radar["data"], texture_xradar["data"])


# ---------------------------------------------------------------------------
# 2. kdp_maesaka -- substituted: the fixture has no PHIDP field (see module
#    docstring), so this uses calculate_snr_from_reflectivity per the
#    brief's own suggested substitution.
# ---------------------------------------------------------------------------
def test_calculate_snr_from_reflectivity(radar_pair):
    radar, xradar_obj = radar_pair
    snr_radar = pyart.retrieve.calculate_snr_from_reflectivity(
        radar, refl_field="reflectivity"
    )
    snr_xradar = pyart.retrieve.calculate_snr_from_reflectivity(
        xradar_obj, refl_field="reflectivity"
    )
    assert_allclose(snr_radar["data"], snr_xradar["data"])


# ---------------------------------------------------------------------------
# 3. despeckle_field
# ---------------------------------------------------------------------------
def test_despeckle_field(radar_pair):
    radar, xradar_obj = radar_pair
    gatefilter_radar = pyart.correct.despeckle_field(radar, "reflectivity_masked")
    gatefilter_xradar = pyart.correct.despeckle_field(xradar_obj, "reflectivity_masked")
    assert_array_equal(gatefilter_radar.gate_excluded, gatefilter_xradar.gate_excluded)


# ---------------------------------------------------------------------------
# 4. calc_zdr_offset
# ---------------------------------------------------------------------------
def test_calc_zdr_offset(radar_pair):
    radar, xradar_obj = radar_pair
    result_radar = pyart.correct.calc_zdr_offset(
        radar,
        gatefilter=pyart.filters.GateFilter(radar),
        zdr_var="differential_reflectivity",
    )
    result_xradar = pyart.correct.calc_zdr_offset(
        xradar_obj,
        gatefilter=pyart.filters.GateFilter(xradar_obj),
        zdr_var="differential_reflectivity",
    )
    assert_allclose(result_radar["bias"], result_xradar["bias"])
    assert_allclose(result_radar["profile_zdr"], result_xradar["profile_zdr"])


# ---------------------------------------------------------------------------
# 5. GateFilter.exclude_below + Xradar.add_filter smoke path
# ---------------------------------------------------------------------------
def test_gatefilter_exclude_below(radar_pair):
    radar, xradar_obj = radar_pair
    gatefilter_radar = pyart.filters.GateFilter(radar)
    gatefilter_radar.exclude_below("reflectivity", 15)
    gatefilter_xradar = pyart.filters.GateFilter(xradar_obj)
    gatefilter_xradar.exclude_below("reflectivity", 15)
    assert_array_equal(gatefilter_radar.gate_excluded, gatefilter_xradar.gate_excluded)


def test_xradar_add_filter_smoke(radar_pair):
    _, xradar_obj = radar_pair
    gatefilter_xradar = pyart.filters.GateFilter(xradar_obj)
    gatefilter_xradar.exclude_below("reflectivity", 15)
    xradar_obj.add_filter(gatefilter_xradar, replace_existing=False)
    assert "filtered_reflectivity" in xradar_obj.fields
    filtered = xradar_obj.fields["filtered_reflectivity"]["data"]
    assert np.ma.is_masked(filtered) or np.ma.getmaskarray(filtered).any()


# ---------------------------------------------------------------------------
# 6. get_nyquist_vel -- skipped: the fixture (like the two real files
#    surveyed while designing it) carries no nyquist_velocity instrument
#    parameter on either side after a cfradial round trip, so there is
#    nothing to compare (both raise/lack the value symmetrically).
# ---------------------------------------------------------------------------
def test_get_nyquist_vel_skipped_when_absent(radar_pair):
    radar, xradar_obj = radar_pair
    if radar.instrument_parameters and "nyquist_velocity" in (
        radar.instrument_parameters
    ):
        pytest.skip("fixture unexpectedly has nyquist_velocity; case not exercised")
    with pytest.raises(LookupError):
        radar.get_nyquist_vel(0)
    with pytest.raises(LookupError):
        xradar_obj.get_nyquist_vel(0)


# ---------------------------------------------------------------------------
# 7. Iterators: iter_field, iter_azimuth, iter_slice
# ---------------------------------------------------------------------------
def test_iterators_equal(radar_pair):
    radar, xradar_obj = radar_pair

    for sweep, (field_radar, field_xradar) in enumerate(
        zip(radar.iter_field("reflectivity"), xradar_obj.iter_field("reflectivity"))
    ):
        assert_allclose(
            np.ma.filled(field_radar, np.nan), field_xradar, err_msg=f"sweep {sweep}"
        )

    for sweep, (az_radar, az_xradar) in enumerate(
        zip(radar.iter_azimuth(), xradar_obj.iter_azimuth())
    ):
        assert_allclose(az_radar, az_xradar, err_msg=f"sweep {sweep}")

    for sweep, (slice_radar, slice_xradar) in enumerate(
        zip(radar.iter_slice(), xradar_obj.iter_slice())
    ):
        assert slice_radar == slice_xradar, f"sweep {sweep}"


# ---------------------------------------------------------------------------
# 8. get_slice / get_start / get_end / get_start_end
# ---------------------------------------------------------------------------
def test_slice_start_end_equal(radar_pair):
    radar, xradar_obj = radar_pair
    for sweep in range(radar.nsweeps):
        assert radar.get_slice(sweep) == xradar_obj.get_slice(sweep)
        assert radar.get_start(sweep) == xradar_obj.get_start(sweep)
        assert radar.get_end(sweep) == xradar_obj.get_end(sweep)
        assert radar.get_start_end(sweep) == xradar_obj.get_start_end(sweep)


# ---------------------------------------------------------------------------
# 9. RHI handling
# ---------------------------------------------------------------------------
def test_rhi_scan_type(rhi_pair):
    radar, xradar_obj = rhi_pair
    assert radar.scan_type == "rhi"
    assert xradar_obj.scan_type == "rhi"


def test_rhi_despeckle_field(rhi_pair):
    radar, xradar_obj = rhi_pair
    gatefilter_radar = pyart.correct.despeckle_field(radar, "reflectivity")
    gatefilter_xradar = pyart.correct.despeckle_field(xradar_obj, "reflectivity")
    assert gatefilter_radar.gate_excluded.shape == gatefilter_xradar.gate_excluded.shape
    assert_array_equal(gatefilter_radar.gate_excluded, gatefilter_xradar.gate_excluded)


def test_rhi_add_field(rhi_pair):
    """add_field must not crash on RHI Xradar objects (each RHI sweep has a
    single azimuth, so hardcoding "azimuth" for drop_duplicates/dims raises
    a ValueError on conflicting sizes). Values added through the Xradar
    accessor should match the equivalent Radar.add_field path."""
    radar, xradar_obj = rhi_pair
    new_data = radar.fields["reflectivity"]["data"] + 1.0
    dic = pyart.config.get_metadata("reflectivity")
    dic["data"] = new_data.copy()
    radar.add_field("reflectivity_plus_one", dic, replace_existing=True)

    xradar_dic = pyart.config.get_metadata("reflectivity")
    xradar_dic["data"] = new_data.copy()
    xradar_obj.add_field("reflectivity_plus_one", xradar_dic, replace_existing=True)

    assert "reflectivity_plus_one" in xradar_obj.fields
    assert_allclose(
        radar.fields["reflectivity_plus_one"]["data"],
        xradar_obj.fields["reflectivity_plus_one"]["data"],
    )
    for sweep in range(radar.nsweeps):
        assert_allclose(
            radar.get_field(sweep, "reflectivity_plus_one"),
            xradar_obj.get_field(sweep, "reflectivity_plus_one"),
            err_msg=f"sweep {sweep}",
        )


# ---------------------------------------------------------------------------
# 10. Gridding: grid_from_radars(Radar) vs grid_from_radars(Xradar)
# ---------------------------------------------------------------------------
def test_grid_from_radars_equal(grid_pair):
    grid_from_radar, grid_from_xradar = grid_pair
    assert_allclose(
        grid_from_radar.fields["reflectivity"]["data"],
        grid_from_xradar.fields["reflectivity"]["data"],
    )


# ---------------------------------------------------------------------------
# 11. Xgrid.to_xarray() vs Grid.to_xarray()
# ---------------------------------------------------------------------------
def test_xgrid_to_xarray_matches_grid(xgrid_pair):
    grid, xgrid = xgrid_pair
    grid_ds = grid.to_xarray()
    xgrid_ds = xgrid.to_xarray()
    assert_allclose(grid_ds["reflectivity"].values, xgrid_ds["reflectivity"].values)


# ---------------------------------------------------------------------------
# 12. Xgrid.add_field then read back
# ---------------------------------------------------------------------------
def test_xgrid_add_field_roundtrip(xgrid_pair):
    _, xgrid = xgrid_pair
    new_field = dict(xgrid.fields["reflectivity"])
    new_field["data"] = new_field["data"] + 1.0
    xgrid.add_field("reflectivity_plus1", new_field, replace_existing=True)
    assert "reflectivity_plus1" in xgrid.fields
    assert_allclose(
        xgrid.fields["reflectivity_plus1"]["data"],
        xgrid.fields["reflectivity"]["data"] + 1.0,
    )


# ---------------------------------------------------------------------------
# 13. get_point_longitude_latitude
# ---------------------------------------------------------------------------
def test_xgrid_get_point_longitude_latitude_matches_grid(xgrid_pair):
    grid, xgrid = xgrid_pair
    lon_grid, lat_grid = grid.get_point_longitude_latitude()
    lon_xgrid, lat_xgrid = xgrid.get_point_longitude_latitude()
    assert_allclose(lon_grid, lon_xgrid)
    assert_allclose(lat_grid, lat_xgrid)
