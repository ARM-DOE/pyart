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
SOB_SINGLE = os.path.join(DIR, 'sob_kdp_single_reference.npy')
PHASE_SINGLE = os.path.join(DIR, 'reproc_phase_single_reference.npy')

####################
# Phase proc tests #
####################


def test_det_sys_phase():
    radar = pyart.io.read_netcdf(os.path.join(DIR, 'sample.nc'))
    assert round(phase_proc.det_sys_phase(radar), 2) == 126.02


""" # takes too long to run
def test_phase_rsl():

    # read in the data
    radar = pyart.io.read_rsl(RSLNAME)

    # process phase
    gates = radar.range['data'][1] - radar.range['data'][0]
    rge = 10.0 * 1000.0
    ng = rge / gates

    reproc_phase, sob_kdp = phase_proc.phase_proc_lp(
        radar, 8.6, sys_phase=332.0, overide_sys_phase=True, debug=True,
        nowrap=ng)

    # compare to known good data
    ref_sob_kdp = np.load(SOBNAME)
    ref_reproc_phase = np.load(PHASENAME)
    assert_array_equal(ref_reproc_phase, reproc_phase['data'])
    assert_array_equal(ref_sob_kdp, sob_kdp['data'])
"""


def test_phase_rsl_fast():

    # read in the data
    radar = pyart.io.read_rsl(RSLNAME)

    # hack to make the radar object appear to only have a single sweep
    radar.sweep_start_ray_index['data'] = np.array([0])
    radar.nsweeps = 1
    data = radar.fields['dp_phase_shift']['data']
    radar.fields['dp_phase_shift']['data'] = data[:360, :]
    data = radar.fields['copol_coeff']['data']
    radar.fields['copol_coeff']['data'] = data[:360, :]

    # process phase
    gates = radar.range['data'][1] - radar.range['data'][0]
    rge = 10.0 * 1000.0
    ng = rge / gates
    reproc_phase, sob_kdp = phase_proc.phase_proc_lp(
        radar, 8.6, sys_phase=332.0, overide_sys_phase=True, debug=True,
        nowrap=ng)

    # compare to known good data
    ref_sob_kdp = np.load(SOB_SINGLE)
    ref_reproc_phase = np.load(PHASE_SINGLE)
    assert_array_equal(ref_reproc_phase, reproc_phase['data'])
    assert_array_equal(ref_sob_kdp, sob_kdp['data'])
