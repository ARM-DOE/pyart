import xradar as xd
from open_radar_data import DATASETS

import pyart

filename = DATASETS.fetch("cfrad.20080604_002217_000_SPOL_v36_SUR.nc")


def test_get_field(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        first_dim="time",
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)
    reflectivity = radar.get_field(0, "DBZ")
    assert reflectivity.shape == (483, 996)


def test_get_gate_x_y_z(filename=filename):
    dtree = xd.io.open_cfradial1_datatree(
        filename,
        first_dim="time",
        optional=False,
    )
    radar = pyart.xradar.Xradar(dtree)
    x, y, z = radar.get_gate_x_y_z(0)
    assert x.shape == (483, 996)
    assert y.shape == (483, 996)
    assert z.shape == (483, 996)
