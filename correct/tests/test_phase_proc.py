""" Tests for the phase_proc module in pyart.correct """

import os.path

import pyart
from pyart.correct import phase_proc

import netCDF4
import numpy as np
from numpy.testing import assert_array_equal

DIR = os.path.dirname(__file__)
RSLNAME = os.path.join(DIR, "sample.sigmet")
SOBNAME = os.path.join(DIR, 'sob_kdp_reference.npy')
PHASENAME = os.path.join(DIR, 'reproc_phase_reference.npy')

####################
# Phase proc tests #
####################

radar = pyart.io.radar.Radar(netCDF4.Dataset('sample.nc'))


def test_det_sys_phase():
    assert round(phase_proc.det_sys_phase(radar), 2) == 126.02


def test_phase_rsl():
    """ Phase Correct a file read using RSL """

    # read in the data
    radarobj = pyart.io.py4dd.RSL_anyformat_to_radar(RSLNAME)
    radar = pyart.io.radar.Radar(radarobj)

    # process phase
    gates = radar.range['data'][1] - radar.range['data'][0]
    rge = 10.0 * 1000.0
    ng = rge / gates

    reproc_phase, sob_kdp = phase_proc.phase_proc(
        radar, 8.6, sys_phase=332.0, overide_sys_phase=True, debug=True,
        nowrap=ng)

    # compare to known good data
    ref_sob_kdp = np.load(SOBNAME)
    ref_reproc_phase = np.load(PHASENAME)
    assert_array_equal(ref_reproc_phase, reproc_phase['data'])
    assert_array_equal(ref_sob_kdp, sob_kdp['data'])


def test_phase_rsl_fast():
    """ Phase Correct a file read using RSL """

    # read in the data
    radarobj = pyart.io.py4dd.RSL_anyformat_to_radar(RSLNAME)
    v = radarobj.contents.volumes[0]
    v.h.nsweeps = 1

    radar = pyart.io.radar.Radar(radarobj)

    # process phase
    gates = radar.range['data'][1] - radar.range['data'][0]
    rge = 10.0 * 1000.0
    ng = rge / gates
    reproc_phase, sob_kdp = phase_proc.phase_proc(
        radar, 8.6, sys_phase=332.0, overide_sys_phase=True, debug=True,
        nowrap=ng)

    # compare to known good data
    ref_sob_kdp = np.load('sob_kdp_single_reference.npy')
    ref_reproc_phase = np.load('reproc_phase_single_reference.npy')
    assert_array_equal(ref_reproc_phase, reproc_phase['data'])
    assert_array_equal(ref_sob_kdp, sob_kdp['data'])
