""" Unit Tests for Py-ART's util/image_mute.py module. """

import pyart
import numpy as np

def test_image_mute_radar():
    # read in example file and image mute
    radar = pyart.io.read_nexrad_archive(pyart.testing.NEXRAD_ARCHIVE_MSG31_FILE)
    radar = pyart.util.image_mute_radar(radar, field='reflectivity',
                                        mute_field='cross_correlation_ratio', mute_threshold=0.97)

    # check that fields are added to the radar object
    assert 'nonmuted_reflectivity' in radar.fields.keys()
    assert 'muted_reflectivity' in radar.fields.keys()

    # check that the number of points in the muted and non muted reflectivity
    # fields have the same number of points as the rhoHV field
    n_rhohv = np.sum(~radar.fields['cross_correlation_ratio']['data'].mask)
    n_mutez = np.sum(~radar.fields['muted_reflectivity']['data'].mask)
    n_nonmutez = np.sum(~radar.fields['nonmuted_reflectivity']['data'].mask)

    assert n_mutez + n_nonmutez == n_rhohv

