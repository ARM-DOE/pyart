""" Tests for the phase_proc module in pyart.correct """

import os.path

import pyart
from pyart.correct import phase_proc

import numpy as np
from numpy.testing import assert_array_equal

DIR = os.path.dirname(__file__)
RSLNAME = os.path.join(DIR, "sample.sigmet")
SOBNAME = os.path.join(DIR, 'sob_kdp_reference.npy')
PHASENAME = os.path.join(DIR, 'reproc_phase_reference.npy')

####################
# Phase proc tests #
####################


def test_phase_rsl():
    """ Phase Correct a file read using RSL """

    # read in the data
    radarobj = pyart.io.py4dd.RSL_anyformat_to_radar(RSLNAME)
    radar = pyart.io.radar.Radar(radarobj)

    # process phase
    gates = radar.range['data'][1] - radar.range['data'][0]
    rge = 10.0 * 1000.0
    ng = rge / gates
    lp = phase_proc.phase_proc(radar, 8.6, sys_phase=332.0,
                               overide_sys_phase=True, debug=True, nowrap=ng)
    reproc_phase, sob_kdp = lp(debug=True)

    # compare to known good data
    ref_sob_kdp = np.load(SOBNAME)
    ref_reproc_phase = np.load(PHASENAME)
    assert_array_equal(ref_reproc_phase, reproc_phase['data'])
    assert_array_equal(ref_sob_kdp, sob_kdp['data'])
