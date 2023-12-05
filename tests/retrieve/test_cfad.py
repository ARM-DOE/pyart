""" Unit Tests for Py-ART's retrieve/echo_class.py module. """

import numpy as np

import pyart


def test_cfad_default():

    # set row to mask and test
    verify_index = 5
    # create random grid of reflectivity data
    ref_random_full = np.random.random((10, 10)) * 30
    # create a mask to mask 90% data from a specific altitude
    mask = np.zeros_like(ref_random_full)
    mask[verify_index, 1:] = 1
    # mask reflectivity data
    ref_random_mask = np.ma.masked_where(mask, ref_random_full)
    # create altitude data
    z_col = np.linspace(0, 12000, 10)
    z_full = np.repeat(z_col[..., np.newaxis], 10, axis=1)
    z_mask = np.ma.masked_where(ref_random_mask.mask, z_full)
    # compute CFAD
    freq_norm, height_edges, field_edges = pyart.retrieve.create_cfad(ref_random_mask,
                                                                      z_mask,
                                                                      field_bins=np.linspace(0, 30, 20),
                                                                      altitude_bins=np.linspace(0, 12000, 10)
                                                                      )
    # if CFAD code works correctly, all values in this row should be false since this altitude has only 1 value (less
    # than the necessary fraction needed)
    assert freq_norm[verify_index, :].mask.all()
