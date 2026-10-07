"""Unit tests for pyart.aux_io.odim_h5 module."""

import h5py
import numpy as np
import pytest
from numpy.testing import assert_allclose

import pyart

NRAYS = 4
NBINS = 5
RSCALE = 500.0


def _make_odim_file(path, conventions, rstart):
    """Create a minimal single-sweep ODIM_H5 PVOL file."""
    with h5py.File(path, "w") as hfile:
        hfile.attrs["Conventions"] = np.bytes_(conventions)

        what = hfile.create_group("what")
        what.attrs["object"] = np.bytes_("PVOL")
        what.attrs["source"] = np.bytes_("NOD:test")

        where = hfile.create_group("where")
        where.attrs["lat"] = 45.0
        where.attrs["lon"] = -90.0
        where.attrs["height"] = 300.0

        dataset1 = hfile.create_group("dataset1")
        ds_where = dataset1.create_group("where")
        ds_where.attrs["nrays"] = NRAYS
        ds_where.attrs["nbins"] = NBINS
        ds_where.attrs["elangle"] = 0.5
        ds_where.attrs["rstart"] = rstart
        ds_where.attrs["rscale"] = RSCALE

        ds_what = dataset1.create_group("what")
        ds_what.attrs["startdate"] = np.bytes_("20260101")
        ds_what.attrs["starttime"] = np.bytes_("000000")
        ds_what.attrs["enddate"] = np.bytes_("20260101")
        ds_what.attrs["endtime"] = np.bytes_("000030")

        data1 = dataset1.create_group("data1")
        d_what = data1.create_group("what")
        d_what.attrs["quantity"] = np.bytes_("DBZH")
        data1.create_dataset("data", data=np.ones((NRAYS, NBINS), dtype="uint8"))


@pytest.mark.filterwarnings("ignore:Py-ART's ODIM module is deprecated")
def test_read_odim_h5_v24_rstart_meters(tmp_path):
    # ODIM_H5 2.4 specifies rstart in metres (CfRadial2 / SI alignment)
    filename = str(tmp_path / "odim_v24.h5")
    _make_odim_file(filename, "ODIM_H5/V2_4", rstart=125.0)
    radar = pyart.aux_io.read_odim_h5(filename)

    assert radar.range["meters_to_center_of_first_gate"] == 125.0
    assert_allclose(radar.range["data"][0], 125.0 + RSCALE / 2)


@pytest.mark.filterwarnings("ignore:Py-ART's ODIM module is deprecated")
def test_read_odim_h5_legacy_rstart_km(tmp_path):
    # ODIM_H5 <= 2.3 specifies rstart in kilometres
    filename = str(tmp_path / "odim_v22.h5")
    _make_odim_file(filename, "ODIM_H5/V2_2", rstart=0.125)
    radar = pyart.aux_io.read_odim_h5(filename)

    assert radar.range["meters_to_center_of_first_gate"] == 125.0
    assert_allclose(radar.range["data"][0], 125.0 + RSCALE / 2)


@pytest.mark.filterwarnings("ignore:Py-ART's ODIM module is deprecated")
def test_read_odim_h5_zero_rstart(tmp_path):
    # rstart = 0 is unaffected by the unit convention
    filename = str(tmp_path / "odim_zero.h5")
    _make_odim_file(filename, "ODIM_H5/V2_4", rstart=0.0)
    radar = pyart.aux_io.read_odim_h5(filename)

    assert radar.range["meters_to_center_of_first_gate"] == 0.0
    assert_allclose(radar.range["data"][0], RSCALE / 2)
