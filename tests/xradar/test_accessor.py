import numpy as np
import xradar as xd
from numpy.testing import assert_almost_equal
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
    assert reflectivity.shape == (483, 996)


def test_get_gate_x_y_z(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)
    x, y, z = radar.get_gate_x_y_z(0)
    assert x.shape == (483, 996)
    assert y.shape == (483, 996)
    assert z.shape == (483, 996)


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
    assert_almost_equal(grid.x["data"], np.arange(-100_000, 120_000, 20_000))
    assert_almost_equal(grid.fields["DBZ"]["data"][0, -1, 0], np.array(0.4243435))
