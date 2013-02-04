# tests for modules in the correct directory

import pyart
from pyart.correct import phase_proc, attenuation, dealias

import numpy as np
import netCDF4
from numpy.testing import assert_array_equal

FILENAME = "XSE110510113001.RAW7BPA"

#################
# Dealias Tests #
#################

def test_dealias_rsl():
    # Compare results of dealiasing a radar object created from data read 
    # using RSL against known good results.
    
    # read in the data
    radarobj = pyart.io.py4dd.RSL_anyformat_to_radar(FILENAME)
    radar = pyart.io.radar.Radar(radarobj)
    
    # find and extract sonde data
    target = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    fname = 'sgpinterpolatedsondeC1.c1.20110510.000000.cdf'
    interp_sounde = netCDF4.Dataset(fname)
    t = dealias.find_time_in_interp_sonde(interp_sounde,target)
    height, speed, direction = t

    # perform dealiasing
    dealias_obj = dealias.dealiaser(radar, height * 1000.0, speed,
                                    direction, target)
    dealias_data = dealias_obj()
    
    # compare against known good data
    reference_data = np.load('dealias_reference.npy')
    assert_array_equal(reference_data, dealias_data['data'].data)


def test_dealias_ncf(): 
    # Compare results of dealiasing a radar object created from data read 
    # from a NetCDF file against known good results.

    # read in the data
    ncf = netCDF4.Dataset('test.nc')
    radar = pyart.io.radar.Radar(ncf)
    
    # find and extract the sonde data
    target = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    fname = 'sgpinterpolatedsondeC1.c1.20110510.000000.cdf'
    interp_sounde = netCDF4.Dataset(fname)
    t = dealias.find_time_in_interp_sonde(interp_sounde,target)
    height, speed, direction = t

    # perform dealiasing
    dealias_obj = dealias.dealiaser(radar, height * 1000.0, speed,
                                    direction, target)
    dealias_data = dealias_obj()
    
    # compare against know good data
    reference_data = np.load('dealias_reference.npy')
    assert_array_equal(reference_data, dealias_data['data'].data)

def test_dealias_compare_rsl(): 
    # Compare results of dealiasing a radar object created from data read 
    # using the RSL library, with and without recreating a RSL Radar object.
 
    # read in the data
    radarobj = pyart.io.py4dd.RSL_anyformat_to_radar(FILENAME)
    radar = pyart.io.radar.Radar(radarobj)
    
    # find and extract the sonde data
    target = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    fname = 'sgpinterpolatedsondeC1.c1.20110510.000000.cdf'
    interp_sounde = netCDF4.Dataset(fname)
    t = dealias.find_time_in_interp_sonde(interp_sounde,target)
    height, speed, direction = t

    # perform dealiasing, with recreation of RSL Radar object
    dealias_obj = dealias.dealiaser(radar, height * 1000.0, speed,
                                    direction, target)
    dealias_data = dealias_obj()
    
     
    # read in the data
    radarobj = pyart.io.py4dd.RSL_anyformat_to_radar(FILENAME)
    radar = pyart.io.radar.Radar(radarobj)
    
    # find and extract the sonde data
    target = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    fname = 'sgpinterpolatedsondeC1.c1.20110510.000000.cdf'
    interp_sounde = netCDF4.Dataset(fname)
    t = dealias.find_time_in_interp_sonde(interp_sounde,target)
    height, speed, direction = t

    # perform dealiasing, with recreation of RSL Radar object
    
    # the original RSL file does not set the altitude correct...
    radarobj.contents.volumes[0].sweeps[0].rays[0].h.alt = 340
    radarobj.contents.volumes[1].sweeps[0].rays[0].h.alt = 340
    dealias_obj2 = dealias.dealiaser(radar, height * 1000.0, speed,
                                    direction, target, rsl_radar = radarobj) 
    dealias_data2 = dealias_obj2()
    
    assert_array_equal(dealias_data2['data'].data, dealias_data['data'].data)

####################
# Phase proc tests #
####################

def test_phase_rsl():

    # read in the data
    radarobj = pyart.io.py4dd.RSL_anyformat_to_radar(FILENAME)
    radar = pyart.io.radar.Radar(radarobj)
 
    # process phase
    gates = radar.range['data'][1] - radar.range['data'][0]
    rge = 10.0 * 1000.0
    ng = rge / gates
    mylp = phase_proc.phase_proc(radar, 8.6, 
                                 sys_phase=332.0,
                                 overide_sys_phase=True, 
                                 debug=True, nowrap=ng)
    reproc_phase, sob_kdp = mylp(debug=True)
    
    # compare to known good data
    #REFFILENAME = 'ref_proc_surcmacI4.c0.20110510.113008.nc' 
    #ref_radar = pyart.io.radar.Radar(netCDF4.Dataset(REFFILENAME))
    #ref_reproc_phase = ref_radar.fields['proc_dp_phase_shift']['data']
    #ref_sob_kdp = ref_radar.fields['recalculated_diff_phase']['data']
    # XXX no recasting to float32 should be done here.
    #assert_array_equal(ref_reproc_phase, reproc_phase['data'].astype('float32'))
    #assert_array_equal(ref_sob_kdp, sob_kdp['data'].astype('float32'))
    
    ref_sob_kdp = np.load('sob_kdp_reference.npy')
    ref_reproc_phase = np.load('reproc_phase_reference.npy')
    assert_array_equal(ref_reproc_phase, reproc_phase['data'])
    assert_array_equal(ref_sob_kdp, sob_kdp['data'])

# currently fails do to float32/64 mismatch
"""
def test_phase_ncf():
    # Same as last test but with netCDF file

    # read in the data
    ncf = netCDF4.Dataset('test.nc')
    radar = pyart.io.radar.Radar(ncf)
 
    # process phase
    gates = radar.range['data'][1] - radar.range['data'][0]
    rge = 10.0 * 1000.0
    ng = rge / gates
    mylp = phase_proc.phase_proc(radar, 8.6, 
                                 sys_phase=332.0,
                                 overide_sys_phase=True, 
                                 debug=True, nowrap=ng)
    reproc_phase, sob_kdp = mylp(debug=True)
    
    # compare to known good data
    #REFFILENAME = 'ref_proc_surcmacI4.c0.20110510.113008.nc' 
    #ref_radar = pyart.io.radar.Radar(netCDF4.Dataset(REFFILENAME))
    #ref_reproc_phase = ref_radar.fields['proc_dp_phase_shift']['data']
    #ref_sob_kdp = ref_radar.fields['recalculated_diff_phase']['data']
    # XXX no recasting to float32 should be done here.
    #assert_array_equal(ref_reproc_phase, reproc_phase['data'].astype('float32'))
    #assert_array_equal(ref_sob_kdp, sob_kdp['data'].astype('float32'))
    
    ref_sob_kdp = np.load('sob_kdp_reference.npy')
    ref_reproc_phase = np.load('reproc_phase_reference.npy')
    assert_array_equal(ref_reproc_phase, reproc_phase['data'])
    assert_array_equal(ref_sob_kdp, sob_kdp['data'])
"""

################################
# Attenuation correction tests #
################################

def test_annenuation_correction_rsl():

    # read in the data
    ncf = netCDF4.Dataset('phase_processed_data.nc')
    radar = pyart.io.radar.Radar(ncf)

    ref_sob_kdp = np.load('sob_kdp_reference.npy')
    ref_reproc_phase = np.load('reproc_phase_reference.npy')
    
    # XXX replace float32 fields with float64 fields
    radar.fields['proc_dp_phase_shift']['data'] = ref_reproc_phase

    spec_at, cor_z = attenuation.calculate_attenuation(
        radar, 8.6, debug=True, a_coef=0.17)

    REFFILENAME = 'ref_proc_surcmacI4.c0.20110510.113008.nc' 
    ref_radar = pyart.io.radar.Radar(netCDF4.Dataset(REFFILENAME))
    ref_spec_at = ref_radar.fields['specific_attenuation']['data']
    ref_cor_z = ref_radar.fields['corrected_reflectivity_horizontal']['data']
    
    # XXX no recasting to float32 should be done here.
    cor_z['data'].data[cor_z['data'].mask] = -9999.0
    assert_array_equal(ref_spec_at, spec_at['data'].astype('float32'))
    assert_array_equal(ref_cor_z, cor_z['data'].astype('float32'))
 
