# tests of the Radar object in radar.py created from a sample mdv file
# the following attributes are tested in some form or another for a
# mdv radar object:
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

# The following methods are not tested
# cf2rad                    : creator, move to nc_utils.py
# extract_rsl_pointing
# get_mdv_meta
# mdv2rad                   : creator, move to py_mdv.py
# prtmode
# rsl2rad                   : creator, move to py4dd.py
# streamcf2rad              : creator, move to nc_utils.py
# ray_header_time_to_dict   :

from os.path import join, dirname

import numpy as np
from numpy.ma.core import MaskedArray

import pyart

# read in the sample file and create a a Radar object
fname = join(dirname(__file__), '110635.mdv')
mdv_radar = pyart.io.read_mdv(fname)
#mdv_radar = pyart.io.radar.Radar(pyart.io.py_mdv.read_mdv(fname))


# azimuth attribute
def test_mdv_azimuth():
    """ MDV radar, azimuth attribute """
    assert 'comment' in mdv_radar.azimuth
    assert 'long_name' in mdv_radar.azimuth
    assert 'standard_name' in mdv_radar.azimuth
    assert 'units' in mdv_radar.azimuth
    assert mdv_radar.azimuth['data'][0] == 0
    assert mdv_radar.azimuth['data'][10] == 10


# elevation attribute
def test_mdv_elevation():
    """ MDV radar, elevation attribute """
    assert 'comment' in mdv_radar.elevation
    assert 'long_name' in mdv_radar.elevation
    assert 'standard_name' in mdv_radar.elevation
    assert 'units' in mdv_radar.elevation
    assert mdv_radar.elevation['data'].shape == (6120, )
    assert mdv_radar.elevation['data'][0] == 0.75


# fields attribute
def test_mdv_fields_dictionaries():
    field_names = ['norm_coherent_power',
                   'reflectivity_horizontal',
                   'dp_phase_shift',
                   'doppler_spectral_width',
                   'diff_reflectivity',
                   'mean_doppler_velocity',
                   'copol_coeff',
                   'diff_phase']
    for field_name in field_names:
        description = "MDV radar, fields : %s, required keys" % field_name
        check_required_field_keys.description = description
        yield check_required_field_keys, field_name


def check_required_field_keys(field_name):
    """ Check that the 7 required keys are present in a field """
    assert 'long_name' in mdv_radar.fields[field_name]
    assert 'valid_min' in mdv_radar.fields[field_name]
    assert 'least_significant_digit' in mdv_radar.fields[field_name]
    assert 'units' in mdv_radar.fields[field_name]
    assert 'valid_max' in mdv_radar.fields[field_name]
    assert 'data' in mdv_radar.fields[field_name]
    assert 'standard_name' in mdv_radar.fields[field_name]


def test_mdv_fields_nonstandard_keys():
    """ MDV radar, fields : norm_coherent_power, nonstandard keys """
    assert 'comment' in mdv_radar.fields['norm_coherent_power']


def test_mdv_fields_data_shape():
    field_names = ['norm_coherent_power',
                   'reflectivity_horizontal',
                   'dp_phase_shift',
                   'doppler_spectral_width',
                   'diff_reflectivity',
                   'mean_doppler_velocity',
                   'copol_coeff',
                   'diff_phase']
    for field_name in field_names:
        description = "MDV radar, fields : %s, data shape" % field_name
        check_field_data_shape.description = description
        yield check_field_data_shape, field_name


def check_field_data_shape(field_name):
    assert mdv_radar.fields[field_name]['data'].shape == (6120, 983)


def test_mdv_fields_data_type():
    field_names = ['norm_coherent_power',
                   'reflectivity_horizontal',
                   'dp_phase_shift',
                   'doppler_spectral_width',
                   'diff_reflectivity',
                   'mean_doppler_velocity',
                   'copol_coeff',
                   'diff_phase']
    for field_name in field_names:
        description = "MDV radar, fields : %s, data type" % field_name
        check_field_data_type.description = description
        yield check_field_data_type, field_name


def check_field_data_type(field_name):
    assert type(mdv_radar.fields[field_name]['data']) is MaskedArray


def test_mdv_fields_data_first_points():
    """ MDV radar, fields : first point in data """
    def first_point(field_name):
        return round(mdv_radar.fields[field_name]['data'][0, 0])
    # these values can be found using:
    # [round(mdv_radar.fields[f]['data'][0,0]) for f in mdv_radar.fields]
    assert first_point('norm_coherent_power') == 0.0
    assert first_point('reflectivity_horizontal') == 24.0
    assert first_point('dp_phase_shift') == 144.0
    assert first_point('doppler_spectral_width') == 12.0
    assert first_point('diff_reflectivity') == -0.0
    assert first_point('mean_doppler_velocity') == 9.0
    assert first_point('copol_coeff') == 1.0
    assert first_point('diff_phase') == 0.0


# location attribute
def test_mdv_location():
    elements = ['latitude', 'altitude', 'longitude']
    for element in elements:
        description = "MDV radar, location : %s " % element
        check_location_element.description = description
        yield check_location_element, element


def check_location_element(element):
    """ Check that location attributes dictionaries have all required keys. """
    assert 'data' in mdv_radar.location[element]
    assert 'standard_name' in mdv_radar.location[element]
    assert 'units' in mdv_radar.location[element]


def test_mdv_location_data():
    """ MDV radar, location : data """
    assert round(mdv_radar.location['latitude']['data']) == 37.0
    assert round(mdv_radar.location['longitude']['data']) == -97.0
    assert round(mdv_radar.location['altitude']['data']) == 328.0


# metadata attribute
def test_mdv_metadata():
    """ MDV radar, metadata attribute """
    assert 'instrument_name' in mdv_radar.metadata
    assert 'source' in mdv_radar.metadata


# naz attribute
def test_mdv_naz():
    """ MDV radar, naz attribute """
    assert mdv_radar.naz == 360


# nele attribute
def test_mdv_nele():
    """ MDV radar, nele attribute """
    assert mdv_radar.nele == 17


# ngates attribute
def test_mdv_ngates():
    """ MDV radar, ngates attribute """
    assert mdv_radar.ngates == 983


# nsweeps attribute
def test_mdv_nsweeps():
    """ MDV radar, nsweeps attribute """
    assert mdv_radar.nsweeps == 17


# range attribute
def test_mdv_range():
    """ MDV radar, range attribute """
    assert 'comment' in mdv_radar.range
    assert 'long_name' in mdv_radar.range
    assert 'standard_name' in mdv_radar.range
    assert 'meters_to_center_of_first_gate' in mdv_radar.range
    assert 'meters_between_gates' in mdv_radar.range
    assert 'units' in mdv_radar.range
    assert 'data' in mdv_radar.range
    assert 'spacing_is_constant' in mdv_radar.range
    assert mdv_radar.range['data'].shape == (983, )
    assert round(mdv_radar.range['data'][0]) == 118.0


# scan_type attribute
def test_mdv_scan_type():
    """ MDV radar, scan type attribute """
    assert mdv_radar.scan_type == 'ppi'


# sweep_info attribute
def test_mdv_sweep_info():
    """ MDV radar, sweep_info attribute """
    assert 'sweep_start_ray_index' in mdv_radar.sweep_info.keys()
    assert 'sweep_mode' in mdv_radar.sweep_info.keys()
    assert 'sweep_number' in mdv_radar.sweep_info.keys()
    assert 'sweep_end_ray_index' in mdv_radar.sweep_info.keys()
    assert 'fixed_angle' in mdv_radar.sweep_info.keys()


def test_mdv_sweep_info_start():
    """ MDV radar, sweep_info : sweep_start_ray_index """
    ssri = mdv_radar.sweep_info['sweep_start_ray_index']
    assert 'data' in ssri.keys()
    assert 'long_name' in ssri.keys()
    assert 'units' in ssri.keys()
    assert np.all(ssri['data'] == np.arange(0, 5761, 360))


def test_mdv_sweep_info_mode():
    """ MDV radar, sweep_info : sweep_mode """
    sm = mdv_radar.sweep_info['sweep_mode']
    assert 'units' in sm.keys()
    assert 'long_name' in sm.keys()
    assert 'data' in sm.keys()
    assert 'comment' in sm.keys()
    assert np.all(sm['data'] == ['azimuth_surveillance    '])


def test_mdv_sweep_info_number():
    """ MDV radar, sweep_info : sweep_number """
    sn = mdv_radar.sweep_info['sweep_number']
    assert 'data' in sn.keys()
    assert 'long_name' in sn.keys()
    assert 'units' in sn.keys()
    assert np.all(sn['data'] == range(17))


def test_mdv_sweep_info_end():
    """ MDV radar, sweep_info : sweep_end_ray_index """
    seri = mdv_radar.sweep_info['sweep_end_ray_index']
    assert 'data' in seri.keys()
    assert 'long_name' in seri.keys()
    assert 'units' in seri.keys()
    assert np.all(seri['data'] == np.arange(359, 6120, 360))


def test_mdv_sweep_info_angle():
    """ MDV radar, sweep_info : fixed_angle """
    fa = mdv_radar.sweep_info['fixed_angle']
    assert 'data' in fa.keys()
    assert 'long_name' in fa.keys()
    assert 'units' in fa.keys()
    assert 'standard_name' in fa.keys()
    assert fa['data'].shape == (17, )


# sweep mode attribute
def test_mdv_sweep_mode():
    """ MDV radar, sweep_mode attribute """
    assert np.all(mdv_radar.sweep_mode == ['ppi'] * 17)


# sweep_number attribute
def test_mdv_sweep_number():
    """ MDV radar, sweep_number attribute """
    assert np.all(mdv_radar.sweep_number == range(17))


# time attribute
def test_mdv_time():
    """ MDV radar, time attribute """
    assert 'comment' in mdv_radar.time.keys()
    assert 'long_name' in mdv_radar.time.keys()
    assert 'standard_name' in mdv_radar.time.keys()
    assert 'units' in mdv_radar.time.keys()
    assert 'calendar' in mdv_radar.time.keys()
    assert 'data' in mdv_radar.time.keys()
    assert mdv_radar.time['units'] == 'seconds since 2011-05-20T11:01:00.0Z'
    assert mdv_radar.time['data'].shape == (6120, )
    assert round(mdv_radar.time['data'][600]) == 33.


# inst_params attribute
def test_mdv_inst_params():
    """ MDV radar, inst_params attribute """
    assert 'prt' in mdv_radar.inst_params.keys()
    assert 'unambiguous_range' in mdv_radar.inst_params.keys()
    assert 'prt_mode' in mdv_radar.inst_params.keys()
    assert 'nyquist_velocity' in mdv_radar.inst_params.keys()


def test_mdv_inst_params_prt():
    """ MDV radar, inst_params : ptr """
    p = mdv_radar.inst_params['prt']
    assert 'comments' in p
    assert 'units' in p
    assert 'data' in p
    assert np.all(np.round(p['data'], 8) == 0.00080645)


def test_mdv_inst_params_unambiguous_range():
    """ MDV radar, inst_params : unambiguous_range """
    ur = mdv_radar.inst_params['unambiguous_range']
    assert 'comment' in ur
    assert 'units' in ur
    assert 'data' in ur
    assert np.all(np.round(ur['data']) == 117996)


def test_mdv_inst_params_prt_mode():
    """ MDV radar, inst_params : prt_mode """
    pm = mdv_radar.inst_params['prt_mode']
    assert 'comments' in pm
    assert 'data' in pm
    assert np.all(pm['data'] == 'fixed                   ')


def test_mdv_inst_params_nyquist_velocity():
    """ MDV radar, inst_params : nyquist_velocity """
    nv = mdv_radar.inst_params['nyquist_velocity']
    assert 'comments' in nv
    assert 'units' in nv
    assert 'data' in nv
    assert np.all(np.round(nv['data']) == 17)
