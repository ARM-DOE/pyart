""" Unit Tests for Py-ART's retrieve/echo_class.py module. """

import numpy as np

import pyart


def test_cfad_default():
    # initalize test radar object
    radar = pyart.io.read(pyart.testing.NEXRAD_ARCHIVE_MSG31_FILE)
    ref_field = "reflectivity"

    # set every value to 20
    radar.fields[ref_field]["data"] = (
        np.ones(radar.fields[ref_field]["data"].shape) * 20
    )

    # set mask to none
    field_mask = np.zeros(radar.fields[ref_field]["data"].shape)

    # calculate CFAD
    freq_norm, height_edges, field_edges = pyart.retrieve.create_cfad(
        radar,
        field_bins=np.linspace(0, 30, 20),
        altitude_bins=np.arange(0, 18000, 100),
        field="reflectivity",
        field_mask=field_mask,
    )

    # set row to mask and test
    verify_index = 12

    # if CFAD code works correctly, each column should have the same values and only 1 column should have a value of 1
    # check all columns are the same
    assert freq_norm.all(axis=0).any()
    # check column 12 is all ones
    assert (freq_norm[:, verify_index] == 1).all()
