""" Unit Tests for Py-ART's io/mdv.py module. """

from os.path import join, dirname

import numpy as np
from numpy.ma.core import MaskedArray
import netCDF4

import pyart

# read in the sample file and create a a Radar object
fname = join(dirname(__file__), '110635.mdv')
radar = pyart.io.read_mdv(fname)


# time attribute
def test_time():
    assert 'comment' in radar.time.keys()
    assert 'long_name' in radar.time.keys()
    assert 'standard_name' in radar.time.keys()
    assert 'units' in radar.time.keys()
    assert 'calendar' in radar.time.keys()
    assert 'data' in radar.time.keys()
    assert radar.time['units'] == 'seconds since 2011-05-20T11:01:00Z'
    assert radar.time['data'].shape == (6120, )
    assert round(radar.time['data'][600]) == 33.


# range attribute
def test_range():
    assert 'long_name' in radar.range
    assert 'standard_name' in radar.range
    assert 'meters_to_center_of_first_gate' in radar.range
    assert 'meters_between_gates' in radar.range
    assert 'units' in radar.range
    assert 'data' in radar.range
    assert 'spacing_is_constant' in radar.range
    assert radar.range['data'].shape == (983, )
    assert round(radar.range['data'][0]) == 118.0


# fields attribute is tested later


# metadata attribute
def test_metadata():
    assert 'instrument_name' in radar.metadata
    assert 'source' in radar.metadata


# scan_type attribute
def test_scan_type():
    assert radar.scan_type == 'ppi'


# latitude attribute
def test_latitude():
    assert 'data' in radar.latitude
    assert 'standard_name' in radar.latitude
    assert 'units' in radar.latitude
    assert round(radar.latitude['data']) == 37.0


# longitude attribute
def test_longitude():
    assert 'data' in radar.longitude
    assert 'standard_name' in radar.longitude
    assert 'units' in radar.longitude
    assert round(radar.longitude['data']) == -97.0


# altitude attribute
def test_altitude():
    assert 'data' in radar.altitude
    assert 'standard_name' in radar.altitude
    assert 'units' in radar.altitude
    assert 'positive' in radar.altitude
    assert round(radar.altitude['data']) == 328.0


# altitude_agl attribute
def test_altitude_agl():
    assert radar.altitude_agl is None


# sweep_number attribute
def test_sweep_number():
    assert 'standard_name' in radar.sweep_number
    assert np.all(radar.sweep_number['data'] == range(17))


# sweep_mode attribute
def test_sweep_mode():
    assert 'standard_name' in radar.sweep_mode
    assert radar.sweep_mode['data'].shape == (17, )
    assert np.all(radar.sweep_mode['data'] == ['azimuth_surveillance'])


# fixed_angle attribute
def test_fixed_angle():
    assert 'standard_name' in radar.fixed_angle
    assert 'units' in radar.fixed_angle
    assert radar.fixed_angle['data'].shape == (17, )
    assert round(radar.fixed_angle['data'][0], 2) == 0.75


# sweep_start_ray_index attribute
def test_sweep_start_ray_index():
    assert 'long_name' in radar.sweep_start_ray_index
    assert radar.sweep_start_ray_index['data'].shape == (17, )
    assert round(radar.sweep_start_ray_index['data'][1]) == 360.0


# sweep_end_ray_index attribute
def test_sweep_end_ray_index():
    assert 'long_name' in radar.sweep_end_ray_index
    assert radar.sweep_end_ray_index['data'].shape == (17, )
    assert round(radar.sweep_end_ray_index['data'][0]) == 359.0


# target_scan_rate attribute
def test_target_scan_rate():
    assert radar.target_scan_rate is None


# azimuth attribute
def test_azimuth():
    assert 'standard_name' in radar.azimuth
    assert 'long_name' in radar.azimuth
    assert 'units' in radar.azimuth
    assert 'axis' in radar.azimuth
    assert round(radar.azimuth['data'][0]) == 0.0
    assert round(radar.azimuth['data'][10]) == 10.0


# elevation attribute
def test_elevation():
    assert 'standard_name' in radar.elevation
    assert 'long_name' in radar.azimuth
    assert 'units' in radar.elevation
    assert 'axis' in radar.elevation
    assert radar.elevation['data'].shape == (6120, )
    assert round(radar.elevation['data'][0], 2) == 0.75


# scan_rate attribute
def test_scan_rate():
    assert radar.scan_rate is None


# antenna_transition attribute
def test_antenna_transition():
    assert radar.antenna_transition is None


# instrument_parameters attribute
def test_instument_parameters():
    # instrument_parameter sub-convention
    keys = ['prt', 'unambiguous_range', 'prt_mode', 'nyquist_velocity']
    for k in keys:
        description = 'instrument_parameters: %s' % k
        check_instrument_parameter.description = description
        yield check_instrument_parameter, k


def check_instrument_parameter(param):
    assert param in radar.instrument_parameters
    param_dic = radar.instrument_parameters[param]
    assert param_dic['meta_group'] == 'instrument_parameters'


# radar_calibration attribute
def test_radar_calibration():
    assert radar.radar_calibration is None


# ngates attribute
def test_ngates():
    assert radar.ngates == 983


# nrays attribute
def test_nrays():
    assert radar.nrays == 6120


# nsweeps attribute
def test_nsweeps():
    assert radar.nsweeps == 17


####################
# fields attribute #
####################


def test_field_dics():
    fields = ['norm_coherent_power',
              'reflectivity_horizontal',
              'dp_phase_shift',
              'doppler_spectral_width',
              'diff_reflectivity',
              'mean_doppler_velocity',
              'copol_coeff',
              'diff_phase']
    for field in fields:
        description = "fields : %s, dictionary" % field
        check_field_dic.description = description
        yield check_field_dic, field


def check_field_dic(field):
    """ Check that the required keys are present in a field dictionary. """
    assert 'standard_name' in radar.fields[field]
    assert 'units' in radar.fields[field]
    assert '_FillValue' in radar.fields[field]
    assert 'coordinates' in radar.fields[field]


def test_field_shapes():
    fields = ['norm_coherent_power',
              'reflectivity_horizontal',
              'dp_phase_shift',
              'doppler_spectral_width',
              'diff_reflectivity',
              'mean_doppler_velocity',
              'copol_coeff',
              'diff_phase']
    for field in fields:
        description = "field : %s, shape" % field
        check_field_shape.description = description
        yield check_field_shape, field


def check_field_shape(field):
    assert radar.fields[field]['data'].shape == (6120, 983)


def test_field_types():
    fields = {'norm_coherent_power': MaskedArray,
              'reflectivity_horizontal': MaskedArray,
              'dp_phase_shift': MaskedArray,
              'doppler_spectral_width': MaskedArray,
              'diff_reflectivity': MaskedArray,
              'mean_doppler_velocity': MaskedArray,
              'copol_coeff': MaskedArray,
              'diff_phase': MaskedArray}
    for field, field_type in fields.iteritems():
        description = "field : %s, type" % field
        check_field_type.description = description
        yield check_field_type, field, field_type


def check_field_type(field, field_type):
    assert type(radar.fields[field]['data']) is field_type


def test_field_first_points():
    # these values can be found using:
    # [round(radar.fields[f]['data'][0,0]) for f in radar.fields]
    fields = {'norm_coherent_power': 0.0,
              'reflectivity_horizontal': 24.0,
              'dp_phase_shift': 144.0,
              'doppler_spectral_width': 12.0,
              'diff_reflectivity': -0.0,
              'mean_doppler_velocity': 9.0,
              'copol_coeff': 1.0,
              'diff_phase': 0.0}
    for field, field_value in fields.iteritems():
        description = "field : %s, first point" % field
        check_field_first_point.description = description
        yield check_field_first_point, field, field_value


def check_field_first_point(field, value):
    assert round(radar.fields[field]['data'][0, 0]) == value
