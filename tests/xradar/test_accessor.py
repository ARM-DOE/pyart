import numpy as np
import xarray as xr
import xradar as xd
from numpy.testing import assert_allclose, assert_almost_equal
from open_radar_data import DATASETS

import pyart

filename = DATASETS.fetch("cfrad.20080604_002217_000_SPOL_v36_SUR.nc")


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
        radar.instrument_parameters["latitude"]["data"], np.array(22.52669907)
    )
    assert_allclose(
        radar.instrument_parameters["longitude"]["data"], np.array(120.4335022)
    )


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
