import io

import numpy as np
import pyproj
import pytest
import xarray as xr
import xradar as xd
from numpy.testing import assert_allclose, assert_almost_equal
from open_radar_data import DATASETS

import pyart

filename = DATASETS.fetch("cfrad.20080604_002217_000_SPOL_v36_SUR.nc")
filename_multiple_sweeps = DATASETS.fetch("swx_20120520_0641.nc")
cfradial2_file = DATASETS.fetch("cfrad2.20080604_002217_000_SPOL_v36_SUR.nc")

COMMON_MAP_TO_GRID_ARGS = {
    "grid_shape": (3, 9, 10),
    "grid_limits": ((-400.0, 400.0), (-900.0, 900.0), (-900, 900)),
    "fields": None,
    "gridding_algo": "map_to_grid",
    "roi_func": lambda z, y, x: 30,
}


def test_get_field(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)
    reflectivity = radar.get_field(0, "DBZ")
    assert reflectivity.shape == (480, 996)


def test_get_azimuth(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)
    azimuths = radar.get_azimuth(0)
    assert azimuths.shape == (480,)


def test_instrument_parameters(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)
    assert radar.instrument_parameters["instrument_name"] == "SPOLRVP8"
    assert_allclose(
        radar.instrument_parameters["sweep_fixed_angle"]["data"][0], np.array(0.4999)
    )
    assert radar.instrument_parameters["sweep_group_name"]["data"][0] == "sweep_0"


def test_coords(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)
    assert_allclose(radar.xradar.coords["latitude"].values, np.array(22.52669907))
    assert_allclose(radar.xradar.coords["longitude"].values, np.array(120.4335022))


def test_get_gate_x_y_z(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)
    x, y, z = radar.get_gate_x_y_z(0)
    assert x.shape == (480, 996)
    assert y.shape == (480, 996)
    assert z.shape == (480, 996)


def test_get_gate_lat_lon(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)
    lat, lon, alt = radar.get_gate_lat_lon_alt(0)
    # Check lat, lon, and alt values
    assert lat.shape == (480, 996)
    assert_allclose(lat.min(), 21.183521)
    assert lon.shape == (480, 996)
    assert_allclose(lon.min(), 118.97935)
    assert alt.shape == (480, 996)
    assert_allclose(alt.min(), 45.0000017)


def test_add_field(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)
    new_field = radar.fields["DBZ"]
    radar.add_field("reflectivity", new_field)
    assert "reflectivity" in radar.fields
    assert radar["sweep_0"]["reflectivity"].shape == radar["sweep_0"]["DBZ"].shape
    # Test add field like as well
    radar.add_field_like("reflectivity", "test", new_field["data"])
    assert "test" in radar.fields
    assert "data" in radar.fields["test"]
    assert radar.fields["test"]["units"] == "dBZ"


def test_grid(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)
    grid = pyart.map.grid_from_radars(
        (radar,),
        grid_shape=(1, 11, 11),
        grid_limits=((2000, 2000), (-100_000.0, 100_000.0), (-100_000.0, 100_000.0)),
        fields=["DBZ"],
    )
    assert_allclose(grid.x["data"], np.arange(-100_000, 120_000, 20_000))
    assert_allclose(grid.fields["DBZ"]["data"][0, -1, 0], np.array(-0.511), rtol=1e-03)


def test_extract_sweeps():
    dtree = xd.io.open_cfradial1_datatree(
        filename_multiple_sweeps,
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)
    sweeps = [0, 5, 6, 8]

    test_radar = radar.extract_sweeps(sweeps)

    assert test_radar.nsweeps == 4
    assert "sweep_0" in test_radar.xradar.children.keys()
    assert "sweep_7" not in test_radar.xradar.children.keys()

    assert np.round(test_radar.latitude["data"][0], 4) == 36.4908

    pytest.raises(ValueError, radar.extract_sweeps, [0, 50])
    pytest.raises(ValueError, radar.extract_sweeps, [-1, 1])


def _check_attrs_similar(grid1, grid2, attr):
    print("Checking attribute:", attr)
    dic1 = getattr(grid1, attr)
    dic2 = getattr(grid2, attr)
    _check_dicts_similar(dic1, dic2)


def _check_dicts_similar(dic1, dic2):
    for k, v in dic1.items():
        print("Checking key:", k)
        if k == "data":
            assert_almost_equal(v, dic2[k])
        else:
            assert dic2[k] == v


def test_grid_write_read():
    # test the read_grid and write_grid function by performing a
    # write/read roundtrip and comparing the two Grid objects
    grid1 = pyart.testing.make_target_grid()
    grid1.projection["comment"] = "This is a comment"
    grid1.metadata["comment"] = "This is another comment"

    with pyart.testing.InTemporaryDirectory():
        tmpfile = "tmp_grid2.nc"
        pyart.io.write_grid(tmpfile, grid1)
        grid2 = xr.open_dataset(tmpfile, decode_times=False)
        grid2 = pyart.xradar.Xgrid(grid2)

        # check fields
        for field in grid1.fields.keys():
            _check_dicts_similar(grid1.fields[field], grid2.fields[field])

        # check attributes
        assert "Conventions" in grid2.metadata
        grid2.metadata.pop("Conventions")
        _check_attrs_similar(grid1, grid2, "metadata")

        _check_attrs_similar(grid1, grid2, "time")

        _check_attrs_similar(grid1, grid2, "origin_latitude")
        _check_attrs_similar(grid1, grid2, "origin_longitude")
        _check_attrs_similar(grid1, grid2, "origin_altitude")

        _check_attrs_similar(grid1, grid2, "x")
        _check_attrs_similar(grid1, grid2, "y")
        _check_attrs_similar(grid1, grid2, "z")

        _check_attrs_similar(grid1, grid2, "point_x")
        _check_attrs_similar(grid1, grid2, "point_y")
        _check_attrs_similar(grid1, grid2, "point_z")

        _check_attrs_similar(grid1, grid2, "projection")
        assert grid1.projection["_include_lon_0_lat_0"] is True

        _check_attrs_similar(grid1, grid2, "point_latitude")
        _check_attrs_similar(grid1, grid2, "point_longitude")
        _check_attrs_similar(grid1, grid2, "point_altitude")

        assert grid1.nx == grid2.nx
        assert grid1.ny == grid2.ny
        assert grid1.nz == grid2.nz

        _check_attrs_similar(grid1, grid2, "radar_latitude")
        _check_attrs_similar(grid1, grid2, "radar_longitude")
        _check_attrs_similar(grid1, grid2, "radar_altitude")
        _check_attrs_similar(grid1, grid2, "radar_time")
        assert grid1.radar_name["data"] == grid2.radar_name["data"]
        assert grid1.nradar == grid2.nradar
        grid2.ds.close()


def test_to_xradar_object():
    """Test the to_radar_object"""
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = dtree.pyart.to_radar()
    reflectivity = radar.get_field(0, "DBZ")
    assert reflectivity.shape == (480, 996)
    assert radar.nsweeps == 9
    assert radar.ngates == 996


def test_cfradial2_data():
    """
    Test swapping the dimensions when the first is time
    """
    dtree = xr.open_datatree(
        cfradial2_file,
    )

    # Before the transformation, time should be in the dimensions
    assert "time" in dtree["sweep_0"].dims

    # After the transformation, azimuth should be moved from coord to dim
    radar = dtree.pyart.to_radar()
    assert "azimuth" in radar["sweep_0"].dims

    # Ensure georeferencing works as expected
    x, y, z = radar.get_gate_x_y_z(0)
    assert x.shape == (480, 996)


def test_grid_xradar():
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = dtree.pyart.to_radar()
    grid = pyart.map.grid_from_radars([radar], **COMMON_MAP_TO_GRID_ARGS)
    assert (10, 9, 3) == (grid.nx, grid.ny, grid.nz)
    assert_allclose(np.min(grid.fields["DBZ"]["data"]), -1.3201784)


def test_to_pyart_radar_passes_radar_through():
    radar = pyart.testing.make_target_radar()
    result = pyart.xradar.to_pyart_radar(radar)
    assert result is radar


def test_to_pyart_radar_passes_xradar_through():
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    xradar_obj = pyart.xradar.Xradar(dtree)
    result = pyart.xradar.to_pyart_radar(xradar_obj)
    assert result is xradar_obj


def test_to_pyart_radar_wraps_datatree():
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    result = pyart.xradar.to_pyart_radar(dtree)
    assert isinstance(result, pyart.xradar.Xradar)


def test_to_pyart_radar_raises_typeerror_on_unsupported_type():
    with pytest.raises(TypeError):
        pyart.xradar.to_pyart_radar({"not": "a radar"})


def test_get_elevation(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)
    elevations = radar.get_elevation(0)
    s = radar.get_slice(0)
    assert elevations.shape == (480,)
    assert_allclose(elevations, radar.elevation["data"][s])
    # copy=True should return an equal but distinct array
    elevations_copy = radar.get_elevation(0, copy=True)
    assert elevations_copy is not radar.elevation["data"][s]
    assert_allclose(elevations_copy, elevations)


def test_get_gate_area(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)
    area = radar.get_gate_area(0)
    nrays_sweep = radar.get_slice(0).stop - radar.get_slice(0).start
    assert area.shape == (nrays_sweep - 1, radar.ngates - 1)
    assert np.all(area > 0)

    # locked known value for the first gate area of sweep 0
    assert_allclose(area[0, 0], 441.78647, rtol=1e-5)


def test_info(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)

    out = io.StringIO()
    radar.info(level="compact", out=out)
    compact_output = out.getvalue()
    assert "DBZ:" in compact_output
    assert "nsweeps: 9" in compact_output
    assert "nrays:" in compact_output

    out = io.StringIO()
    radar.info(level="standard", out=out)
    standard_output = out.getvalue()
    assert "DBZ:" in standard_output
    assert "nsweeps: 9" in standard_output
    assert "ngates: 996" in standard_output
    assert "sweep_mode:" in standard_output

    with pytest.raises(ValueError):
        radar.info(level="invalid")


def test_rays_per_sweep(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)
    rays_per_sweep = radar.rays_per_sweep
    assert "data" in rays_per_sweep
    expected = (
        radar.sweep_end_ray_index["data"] - radar.sweep_start_ray_index["data"] + 1
    )
    assert_allclose(rays_per_sweep["data"], expected)
    assert rays_per_sweep["data"].shape == (radar.nsweeps,)

    # init_rays_per_sweep should reset the attribute
    radar.init_rays_per_sweep()
    assert_allclose(radar.rays_per_sweep["data"], expected)


def test_init_gate_altitude_alias(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)
    radar.init_gate_altitude()
    assert_allclose(
        radar.gate_altitude["data"], radar.altitude["data"] + radar.gate_z["data"]
    )


def test_xgrid_projection_proj():
    grid = pyart.testing.make_target_grid()

    with pyart.testing.InTemporaryDirectory():
        tmpfile = "tmp_grid_proj.nc"
        pyart.io.write_grid(tmpfile, grid)
        grid_ds = xr.open_dataset(tmpfile, decode_times=False)
        xgrid = pyart.xradar.Xgrid(grid_ds)

        # default pyart_aeqd projection cannot be turned into a Proj
        with pytest.raises(ValueError):
            xgrid.projection_proj

        xgrid.projection["proj"] = "aeqd"
        proj = xgrid.projection_proj
        assert isinstance(proj, pyproj.Proj)
        grid_ds.close()


def test_xgrid_write_roundtrip():
    grid = pyart.testing.make_target_grid()

    with pyart.testing.InTemporaryDirectory():
        tmpfile = "tmp_grid_source.nc"
        pyart.io.write_grid(tmpfile, grid)
        grid_ds = xr.open_dataset(tmpfile, decode_times=False)
        xgrid = pyart.xradar.Xgrid(grid_ds)

        out_file = "tmp_grid_written.nc"
        xgrid.write(out_file)

        grid_read = pyart.io.read_grid(out_file)
        assert_allclose(
            grid_read.fields["reflectivity"]["data"],
            xgrid.fields["reflectivity"]["data"],
        )
        grid_ds.close()


def test_optional_platform_scan_attrs_default_to_none(filename=filename):
    # None of these variables are present in this fixture file, so the
    # attributes must exist (no AttributeError) and be None, matching
    # Radar's behavior for optional attrs absent from the source file.
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)

    optional_attrs = [
        "antenna_transition",
        "rays_are_indexed",
        "ray_angle_res",
        "target_scan_rate",
        "scan_rate",
        "altitude_agl",
        "rotation",
        "tilt",
        "roll",
        "drift",
        "heading",
        "pitch",
        "georefs_applied",
        "radar_calibration",
    ]
    for attr in optional_attrs:
        assert hasattr(radar, attr)
        assert getattr(radar, attr) is None


def test_optional_platform_scan_attrs_populated_when_present(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )

    # Inject a scan-rate variable (per-ray) into each sweep dataset with
    # known values and attributes.
    scan_rate_attrs = {"units": "degrees/second", "long_name": "antenna scan rate"}
    heading_attrs = {"units": "degrees", "long_name": "platform heading"}
    for sweep in dtree.match("sweep_*"):
        sweep_ds = dtree[sweep].to_dataset()
        n = sweep_ds.sizes["azimuth"]
        sweep_ds["scan_rate"] = (("azimuth",), np.full(n, 5.0, dtype="float64"))
        sweep_ds["scan_rate"].attrs = scan_rate_attrs
        sweep_ds["heading"] = (("azimuth",), np.full(n, 45.0, dtype="float64"))
        sweep_ds["heading"].attrs = heading_attrs
        dtree[sweep].ds = sweep_ds

    radar = pyart.xradar.Xradar(dtree)

    assert radar.scan_rate is not None
    assert radar.scan_rate["data"].shape == (radar.nrays,)
    assert_allclose(radar.scan_rate["data"], 5.0)
    assert radar.scan_rate["units"] == "degrees/second"
    assert radar.scan_rate["long_name"] == "antenna scan rate"

    assert radar.heading is not None
    assert radar.heading["data"].shape == (radar.nrays,)
    assert_allclose(radar.heading["data"], 45.0)
    assert radar.heading["units"] == "degrees"
    assert radar.heading["long_name"] == "platform heading"


def test_field_named_z_not_shadowed_by_georeference_coordinate(filename=filename):
    # Regression test for GH #1765: a moment named 'z' (as GAMIC stores
    # reflectivity) must not be dropped/overwritten by the Cartesian height
    # coordinate that xradar.georeference() adds when Xradar computes gate
    # locations.
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    for sweep in dtree.match("sweep_*"):
        if "DBZ" in dtree[sweep].ds:
            dtree[sweep].ds = dtree[sweep].ds.rename({"DBZ": "z"})

    radar = pyart.xradar.Xradar(dtree)

    assert "z" in radar.fields
    # Coordinates must not leak into fields just because their dims match.
    assert "x" not in radar.fields
    assert "y" not in radar.fields
    original_z = radar.fields["z"]["data"].copy()

    # Trigger the internal xradar.georeference() call used to compute gate
    # locations.
    x, y, z = radar.get_gate_x_y_z(0)
    assert z.shape == (480, 996)

    # The 'z' moment must still be discoverable as a field with its
    # original data, not the georeferenced height coordinate.
    refreshed_fields = radar._find_fields(radar.combined_sweeps)
    assert "z" in refreshed_fields
    assert_allclose(refreshed_fields["z"]["data"], original_z)
    assert "x" not in refreshed_fields
    assert "y" not in refreshed_fields
