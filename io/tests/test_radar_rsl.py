# tests of the Radar object in radar.py created from a sample RSL file
# the following attributes are tested in some form or another for a
# rsl radar object:
#   azimuth
#   elevation
#   fields
#   location
#   metadata
#   naz
#   nele
#   ngates
#   nsweeps
#   range
#   scan_type
#   sweep_info
#   sweep_mode
#   sweep_number
#   time
#   inst_params
#   cal             # not defined for mdv files
#   nrays           # " "
#   tu              # " "

# The following methods are not tested
# cf2rad                    : creator, move to nc_utils.py
# extract_rsl_pointing
# get_mdv_meta
# mdv2rad                   : creator, move to py_rsl.py
# prtmode
# rsl2rad                   : creator, move to py4dd.py
# streamcf2rad              : creator, move to nc_utils.py
# ray_header_time_to_dict   :

from os.path import join, dirname

import numpy as np
from numpy.ma.core import MaskedArray

import pyart

# read in the sample file and create a a Radar object
fname = join(dirname(__file__), 'sample_uf.uf')
rsl_radar = pyart.io.radar.Radar(
    pyart.io.py4dd.RSL_anyformat_to_radar(fname))


# azimuth attribute
def test_rsl_azimuth():
    """ RSL radar, azimuth attribute """
    assert 'comment' in rsl_radar.azimuth
    assert 'long_name' in rsl_radar.azimuth
    assert 'standard_name' in rsl_radar.azimuth
    assert 'units' in rsl_radar.azimuth
    assert round(rsl_radar.azimuth['data'][0]) == 360.0
    assert round(rsl_radar.azimuth['data'][10]) == 10.0


# elevation attribute
def test_rsl_elevation():
    """ RSL radar, elevation attribute """
    assert 'comment' in rsl_radar.elevation
    assert 'long_name' in rsl_radar.elevation
    assert 'standard_name' in rsl_radar.elevation
    assert 'units' in rsl_radar.elevation
    assert rsl_radar.elevation['data'].shape == (7920, )
    assert round(rsl_radar.elevation['data'][0], 2) == 0.52


# fields attribute
def test_rsl_fields_dictionaries():
    field_names = ['norm_coherent_power',
                   'reflectivity_horizontal',
                   'dp_phase_shift',
                   'diff_reflectivity',
                   'mean_doppler_velocity',
                   'copol_coeff',
                   'diff_phase']
    for field_name in field_names:
        description = "RSL radar, fields : %s, required keys" % field_name
        check_required_field_keys.description = description
        yield check_required_field_keys, field_name


def check_required_field_keys(field_name):
    """ Check that the 7 required keys are present in a field """
    assert 'long_name' in rsl_radar.fields[field_name]
    assert 'valid_min' in rsl_radar.fields[field_name]
    assert 'least_significant_digit' in rsl_radar.fields[field_name]
    assert 'units' in rsl_radar.fields[field_name]
    assert 'valid_max' in rsl_radar.fields[field_name]
    assert 'data' in rsl_radar.fields[field_name]
    assert 'standard_name' in rsl_radar.fields[field_name]


def test_rsl_fields_nonstandard_keys():
    """ RSL radar, fields : norm_coherent_power, nonstandard keys """
    assert 'comment' in rsl_radar.fields['norm_coherent_power']


def test_rsl_fields_data_shape():
    field_names = ['norm_coherent_power',
                   'reflectivity_horizontal',
                   'dp_phase_shift',
                   'diff_reflectivity',
                   'mean_doppler_velocity',
                   'copol_coeff',
                   'diff_phase']
    for field_name in field_names:
        description = "RSL radar, fields : %s, data shape" % field_name
        check_field_data_shape.description = description
        yield check_field_data_shape, field_name


def check_field_data_shape(field_name):
    assert rsl_radar.fields[field_name]['data'].shape == (7920, 801)


def test_rsl_fields_data_type():
    field_names = ['norm_coherent_power',
                   'reflectivity_horizontal',
                   'dp_phase_shift',
                   'diff_reflectivity',
                   'mean_doppler_velocity',
                   'copol_coeff',
                   'diff_phase']
    for field_name in field_names:
        description = "RSL radar, fields : %s, data type" % field_name
        check_field_data_type.description = description
        yield check_field_data_type, field_name


def check_field_data_type(field_name):
    assert type(rsl_radar.fields[field_name]['data']) is MaskedArray


def test_rsl_fields_data_first_points():
    """ RSL radar, fields : first point in data """
    def first_point(field_name):
        return round(rsl_radar.fields[field_name]['data'][0, 0])
    # these values can be found using:
    # [round(rsl_radar.fields[f]['data'][0,0]) for f in rsl_radar.fields]
    assert first_point('norm_coherent_power') == 1.0
    assert first_point('reflectivity_horizontal') == -8.0
    assert first_point('dp_phase_shift') == 90.0
    assert first_point('diff_reflectivity') == 0.0
    assert first_point('mean_doppler_velocity') == 0.0
    assert first_point('copol_coeff') == 1.0
    assert first_point('diff_phase') == 0.0


# location attribute
def test_rsl_location():
    elements = ['latitude', 'altitude', 'longitude']
    for element in elements:
        description = "RSL radar, location : %s " % element
        check_location_element.description = description
        yield check_location_element, element


def check_location_element(element):
    """ Check that location attributes dictionaries have all required keys. """
    assert 'data' in rsl_radar.location[element]
    assert 'standard_name' in rsl_radar.location[element]
    assert 'units' in rsl_radar.location[element]


def test_rsl_location_data():
    """ RSL radar, location : data """
    assert round(rsl_radar.location['latitude']['data']) == 37.0
    assert round(rsl_radar.location['longitude']['data']) == -97.0
    assert round(rsl_radar.location['altitude']['data']) == 340.0


# metadata attribute
def test_rsl_metadata():
    """ RSL radar, metadata attribute """
    assert 'instrument_name' in rsl_radar.metadata
    assert 'country' in rsl_radar.metadata
    assert 'original_container' in rsl_radar.metadata
    assert 'project' in rsl_radar.metadata
    assert 'state' in rsl_radar.metadata


# naz attribute
def test_rsl_naz():
    """ RSL radar, naz attribute """
    assert rsl_radar.naz == 360


# nele attribute
def test_rsl_nele():
    """ RSL radar, nele attribute """
    assert rsl_radar.nele == 22


# ngates attribute
def test_rsl_ngates():
    """ RSL radar, ngates attribute """
    assert rsl_radar.ngates == 801


# nsweeps attribute
def test_rsl_nsweeps():
    """ RSL radar, nsweeps attribute """
    assert rsl_radar.nsweeps == 22


# range attribute
def test_rsl_range():
    """ RSL radar, range attribute """
    assert 'comment' in rsl_radar.range
    assert 'long_name' in rsl_radar.range
    assert 'standard_name' in rsl_radar.range
    assert 'meters_to_center_of_first_gate' in rsl_radar.range
    assert 'meters_between_gates' in rsl_radar.range
    assert 'units' in rsl_radar.range
    assert 'data' in rsl_radar.range
    assert 'spacing_is_constant' in rsl_radar.range
    assert rsl_radar.range['data'].shape == (801, )
    assert round(rsl_radar.range['data'][0]) == 0.0


# scan_type attribute
def test_rsl_scan_type():
    """ RSL radar, scan type attribute """
    assert rsl_radar.scan_type == 'ppi'


# sweep_info attribute
def test_rsl_sweep_info():
    """ RSL radar, sweep_info attribute """
    assert 'sweep_start_ray_index' in rsl_radar.sweep_info.keys()
    assert 'sweep_mode' in rsl_radar.sweep_info.keys()
    assert 'sweep_number' in rsl_radar.sweep_info.keys()
    assert 'sweep_end_ray_index' in rsl_radar.sweep_info.keys()
    assert 'fixed_angle' in rsl_radar.sweep_info.keys()


def test_rsl_sweep_info_start():
    """ RSL radar, sweep_info : sweep_start_ray_index """
    ssri = rsl_radar.sweep_info['sweep_start_ray_index']
    assert 'data' in ssri.keys()
    assert 'long_name' in ssri.keys()
    assert 'units' in ssri.keys()
    assert np.all(ssri['data'] == np.arange(0, 7561, 360))


def test_rsl_sweep_info_mode():
    """ RSL radar, sweep_info : sweep_mode """
    sm = rsl_radar.sweep_info['sweep_mode']
    assert 'units' in sm.keys()
    assert 'long_name' in sm.keys()
    assert 'data' in sm.keys()
    assert 'comment' in sm.keys()
    assert sm['data'] == ['azimuth_surveillance    '] * 22


def test_rsl_sweep_info_number():
    """ RSL radar, sweep_info : sweep_number """
    sn = rsl_radar.sweep_info['sweep_number']
    assert 'data' in sn.keys()
    assert 'long_name' in sn.keys()
    assert 'units' in sn.keys()
    assert sn['data'] == range(22)


def test_rsl_sweep_info_end():
    """ RSL radar, sweep_info : sweep_end_ray_index """
    seri = rsl_radar.sweep_info['sweep_end_ray_index']
    assert 'data' in seri.keys()
    assert 'long_name' in seri.keys()
    assert 'units' in seri.keys()
    assert np.all(seri['data'] == np.arange(359, 7920, 360))


def test_rsl_sweep_info_angle():
    """ RSL radar, sweep_info : fixed_angle """
    fa = rsl_radar.sweep_info['fixed_angle']
    assert 'data' in fa.keys()
    assert 'long_name' in fa.keys()
    assert 'units' in fa.keys()
    assert 'standard_name' in fa.keys()
    assert fa['data'].shape == (22, )


# sweep mode attribute
def test_rsl_sweep_mode():
    """ RSL radar, sweep_mode attribute """
    assert np.all(rsl_radar.sweep_mode == ['ppi'] * 22)


# sweep_number attribute
def test_rsl_sweep_number():
    """ RSL radar, sweep_number attribute """
    assert np.all(rsl_radar.sweep_number == range(22))


# time attribute
def test_rsl_time():
    """ RSL radar, time attribute """
    assert 'comment' in rsl_radar.time.keys()
    assert 'long_name' in rsl_radar.time.keys()
    assert 'standard_name' in rsl_radar.time.keys()
    assert 'units' in rsl_radar.time.keys()
    assert 'calendar' in rsl_radar.time.keys()
    assert 'data' in rsl_radar.time.keys()
    assert rsl_radar.time['units'] == 'seconds since 2011-05-23T22:42:59Z'
    assert rsl_radar.time['data'].shape == (7920, )
    assert round(rsl_radar.time['data'][600]) == 25.0


# inst_params attribute
def test_rsl_inst_params():
    """ RSL radar, inst_params attribute """
    assert 'prt' in rsl_radar.inst_params.keys()
    assert 'unambiguous_range' in rsl_radar.inst_params.keys()
    assert 'prt_mode' in rsl_radar.inst_params.keys()
    assert 'nyquist_velocity' in rsl_radar.inst_params.keys()


def test_rsl_inst_params_prt():
    """ RSL radar, inst_params : ptr """
    p = rsl_radar.inst_params['prt']
    assert 'comments' in p
    assert 'units' in p
    assert 'data' in p
    assert np.all(np.round(p['data'], 8) == 0.00046019)


def test_rsl_inst_params_unambiguous_range():
    """ RSL radar, inst_params : unambiguous_range """
    ur = rsl_radar.inst_params['unambiguous_range']
    assert 'comment' in ur
    assert 'units' in ur
    assert 'data' in ur
    assert np.all(np.round(ur['data']) == 68981.)


def test_rsl_inst_params_prt_mode():
    """ RSL radar, inst_params : prt_mode """
    pm = rsl_radar.inst_params['prt_mode']
    assert 'comments' in pm
    assert 'data' in pm
    assert np.all(pm['data'] == 'dual                    ')


def test_rsl_inst_params_nyquist_velocity():
    """ RSL radar, inst_params : nyquist_velocity """
    nv = rsl_radar.inst_params['nyquist_velocity']
    assert 'comments' in nv
    assert 'units' in nv
    assert 'data' in nv
    assert np.all(np.round(nv['data']) == 17)


###############################################
# these attributed are not found in mdv files #
###############################################


# cal attribute
def test_rsl_cal():
    """ RSL radar, cal attribute """
    assert rsl_radar.cal == 'gregorian'


def test_rsl_nrays():
    """ RSL radar, nrays attribute """
    assert rsl_radar.nrays == 360


def test_rsl_tu():
    """ RSL radar, tu attribute """
    assert rsl_radar.tu == 'seconds since 2011-05-23T22:42:59Z'
