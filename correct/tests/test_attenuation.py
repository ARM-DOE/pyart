""" Tests for the attenuation module in pyart.correct """

import os.path

import pyart
from pyart.correct import attenuation

import numpy as np
from numpy.testing import assert_array_equal

DIR = os.path.dirname(__file__)

RSLNAME = os.path.join(DIR, "sample.sigmet")
PHASENAME = os.path.join(DIR, 'reproc_phase_reference.npy')
CORZNAME = os.path.join(DIR, 'cor_z_reference.npy')
SPECATNAME = os.path.join(DIR, 'spec_at_reference.npy')


################################
# Attenuation correction tests #
################################


def test_attenuation_correction_rsl():
    """ Test calculate_attenuation on data read using RSL """

    # read in the data
    radarobj = pyart.io.py4dd.RSL_anyformat_to_radar(RSLNAME)
    radar = pyart.io.radar.Radar(radarobj)

    # add the fields created by phase_proc
    reproc_phase = np.load(PHASENAME)
    radar.fields['proc_dp_phase_shift'] = {'data': reproc_phase}

    spec_at, cor_z = attenuation.calculate_attenuation(
        radar, 8.6, debug=True, a_coef=0.17)

    ref_cor_z = np.load(CORZNAME)
    ref_spec_at = np.load(SPECATNAME)

    # XXX no recasting to float32 should be done here.
    cor_z['data'].data[cor_z['data'].mask] = -9999.0
    assert_array_equal(ref_spec_at, spec_at['data'].astype('float32'))
    assert_array_equal(ref_cor_z, cor_z['data'].astype('float32'))
