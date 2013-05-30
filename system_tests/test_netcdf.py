""" Unit Tests for Py-ART's io/netcdf.py module. """

from os.path import join, dirname

import numpy as np
from numpy.ma.core import MaskedArray
import netCDF4

import pyart

# read in the sample file and create a a Radar object
fname = join(dirname(__file__), 'sample_nc.nc')
radar = pyart.io.read_netcdf(fname)


# time attribute
def test_time():
    assert 'comment' in radar.time.keys()
    assert 'long_name' in radar.time.keys()
    assert 'standard_name' in radar.time.keys()
    assert 'units' in radar.time.keys()
    assert 'calendar' in radar.time.keys()
    assert 'data' in radar.time.keys()
    assert radar.time['units'] == 'seconds since 2012-12-11T22:48:21Z'
    assert radar.time['data'].shape == (3524, )
    assert round(radar.time['data'][600]) == 51.


# range attribute
def test_range():
    assert 'long_name' in radar.range
    assert 'standard_name' in radar.range
    assert 'meters_to_center_of_first_gate' in radar.range
    assert 'meters_between_gates' in radar.range
    assert 'units' in radar.range
    assert 'data' in radar.range
    assert 'spacing_is_constant' in radar.range
    assert radar.range['data'].shape == (1333, )
    assert round(radar.range['data'][0]) == 0.0


# fields attribute is tested later


# metadata attribute
def test_metadata():
    assert 'comment' in radar.metadata
    assert 'scan_id' in radar.metadata
    assert 'platform_is_mobile' in radar.metadata
    assert 'site_name' in radar.metadata
    assert 'title' in radar.metadata
    assert 'scan_name' in radar.metadata
    assert 'n_gates_vary' in radar.metadata
    assert 'Conventions' in radar.metadata
    assert 'source' in radar.metadata
    assert 'version' in radar.metadata
    assert 'references' in radar.metadata
    assert 'instrument_name' in radar.metadata
    assert 'ray_times_increase' in radar.metadata
    assert 'institution' in radar.metadata
    assert 'history' in radar.metadata


# scan_type attribute
def test_scan_type():
    assert radar.scan_type == 'ppi'


# latitude attribute
def test_latitude():
    assert 'data' in radar.latitude
    assert 'standard_name' in radar.latitude
    assert 'units' in radar.latitude
    assert round(radar.latitude['data']) == -2.0


# longitude attribute
def test_longitude():
    assert 'data' in radar.longitude
    assert 'standard_name' in radar.longitude
    assert 'units' in radar.longitude
    assert round(radar.longitude['data']) == 147.0


# altitude attribute
def test_altitude():
    assert 'data' in radar.altitude
    assert 'standard_name' in radar.altitude
    assert 'units' in radar.altitude
    assert 'positive' in radar.altitude
    assert round(radar.altitude['data']) == 126.0


# altitude_agl attribute
def test_altitude_agl():
    assert radar.altitude_agl is None


# sweep_number attribute
def test_sweep_number():
    assert 'standard_name' in radar.sweep_number
    assert '_FillValue' in radar.sweep_number
    assert np.all(radar.sweep_number['data'] == range(9, 19))


# sweep_mode attribute
def test_sweep_mode():
    assert 'standard_name' in radar.sweep_mode
    assert radar.sweep_mode['data'].shape == (10, 32)
    str_array = netCDF4.chartostring(radar.sweep_mode['data'])
    assert np.all(str_array == ['azimuth_surveillance'])


# fixed_angle attribute
def test_fixed_angle():
    assert 'standard_name' in radar.fixed_angle
    assert 'units' in radar.fixed_angle
    assert '_FillValue' in radar.fixed_angle
    assert radar.fixed_angle['data'].shape == (10, )
    assert round(radar.fixed_angle['data'][0], 2) == 6.4


# sweep_start_ray_index attribute
def test_sweep_start_ray_index():
    assert 'standard_name' in radar.sweep_start_ray_index
    assert '_FillValue' in radar.sweep_start_ray_index
    assert radar.sweep_start_ray_index['data'].shape == (10, )
    assert round(radar.sweep_start_ray_index['data'][1]) == 271.0


# sweep_end_ray_index attribute
def test_sweep_end_ray_index():
    assert 'standard_name' in radar.sweep_end_ray_index
    assert '_FillValue' in radar.sweep_end_ray_index
    assert radar.sweep_end_ray_index['data'].shape == (10, )
    assert round(radar.sweep_end_ray_index['data'][0]) == 270.0


# target_scan_rate attribute
def test_target_scan_rate():
    assert radar.target_scan_rate is None


# azimuth attribute
def test_azimuth():
    assert '_FillValue' in radar.azimuth
    assert 'standard_name' in radar.azimuth
    assert 'units' in radar.azimuth
    assert round(radar.azimuth['data'][0]) == 6.0
    assert round(radar.azimuth['data'][10]) == 16.0


# elevation attribute
def test_elevation():
    assert '_FillValue' in radar.elevation
    assert 'positive' in radar.elevation
    assert 'standard_name' in radar.elevation
    assert 'units' in radar.elevation
    assert radar.elevation['data'].shape == (3524, )
    assert round(radar.elevation['data'][0], 2) == 6.4


# scan_rate attribute
def test_scan_rate():
    assert 'standard_name' in radar.scan_rate
    #assert 'meta_group' in radar.scan_rate  # shouldn't be in this group
    assert 'units' in radar.scan_rate
    assert '_FillValue' in radar.scan_rate
    assert np.ma.is_masked(radar.scan_rate['data'][0])


# antenna_transition attribute
def test_antenna_transition():
    assert 'standard_name' in radar.antenna_transition
    assert 'comment' in radar.antenna_transition
    assert np.all(radar.antenna_transition['data'] == 0)


# instrument_parameters attribute
def test_instument_parameters():
    # instrument_parameter sub-convention
    keys = ['frequency', 'follow_mode', 'pulse_width', 'prt_mode', 'prt',
            'prt_ratio', 'polarization_mode', 'nyquist_velocity',
            'unambiguous_range', 'n_samples']
    # radar_parameters sub-convention
    keys += ['radar_antenna_gain_h', 'radar_antenna_gain_v',
             'radar_beam_width_h', 'radar_beam_width_v',
             'measured_transmit_power_h', 'measured_transmit_power_v',
             'radar_rx_bandwidth']
    for k in keys:
        description = 'instrument_parameters: %s' % k
        check_instrument_parameter.description = description
        yield check_instrument_parameter, k


def check_instrument_parameter(param):
    assert param in radar.instrument_parameters
    # the meta_group is not reliably set in CF/Radial files
    #param_dic = radar.instrument_parameters[param]
    #assert param_dic['meta_group'] == 'instrument_parameters'


# radar_calibration attribute
def test_radar_calibration():
    keys = ['r_calib_index',
            'r_calib_time',
            'r_calib_pulse_width',
            'r_calib_antenna_gain_h',   # should be r_calib_ant_gain_h
            'r_calib_antenna_gain_v',   # should be r_calib_ant_gain_v
            'r_calib_xmit_power_h',
            'r_calib_xmit_power_v',
            'r_calib_two_way_waveguide_loss_h',
            'r_calib_two_way_waveguide_loss_v',
            'r_calib_two_way_radome_loss_h',
            'r_calib_two_way_radome_loss_v',
            'r_calib_receiver_mismatch_loss',
            'r_calib_radar_constant_h',
            'r_calib_radar_constant_v',
            'r_calib_noise_hc',
            'r_calib_noise_vc',
            'r_calib_noise_hx',
            'r_calib_noise_vx',
            'r_calib_receiver_gain_hc',
            'r_calib_receiver_gain_vc',
            'r_calib_receiver_gain_hx',
            'r_calib_receiver_gain_vx',
            'r_calib_base_dbz_1km_hc',
            'r_calib_base_dbz_1km_vc',
            'r_calib_base_dbz_1km_hx',
            'r_calib_base_dbz_1km_vx',
            'r_calib_sun_power_hc',
            'r_calib_sun_power_vc',
            'r_calib_sun_power_hx',
            'r_calib_sun_power_vx',
            'r_calib_noise_source_power_h',
            'r_calib_noise_source_power_v',
            'r_calib_power_measure_loss_h',
            'r_calib_power_measure_loss_v',
            'r_calib_coupler_forward_loss_h',
            'r_calib_coupler_forward_loss_v',
            'r_calib_zdr_correction',
            'r_calib_ldr_correction_h',
            'r_calib_ldr_correction_v',
            'r_calib_system_phidp',
            'r_calib_test_power_h',
            'r_calib_test_power_v',
            'r_calib_receiver_slope_hc',
            'r_calib_receiver_slope_vc',
            'r_calib_receiver_slope_hx',
            'r_calib_receiver_slope_vx', ]
    for k in keys:
        description = 'radar_calibration: %s' % k
        check_radar_calibration.description = description
        yield check_radar_calibration, k


def check_radar_calibration(param):
    assert param in radar.radar_calibration
    param_dic = radar.radar_calibration[param]
    assert param_dic['meta_group'] == 'radar_calibration'


# ngates attribute
def test_ngates():
    assert radar.ngates == 1333


# nrays attribute
def test_nrays():
    assert radar.nrays == 3524


# nsweeps attribute
def test_nsweeps():
    assert radar.nsweeps == 10


####################
# fields attribute #
####################


def test_field_dics():
    fields = ['PHIDP', 'WIDTH', 'DBZ', 'DBM', 'RHOHV', 'SNR', 'VEL', 'ZDR']
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
    fields = ['PHIDP', 'WIDTH', 'DBZ', 'DBM', 'RHOHV', 'SNR', 'VEL', 'ZDR']
    for field in fields:
        description = "field : %s, shape" % field
        check_field_shape.description = description
        yield check_field_shape, field


def check_field_shape(field):
    assert radar.fields[field]['data'].shape == (3524, 1333)


def test_field_types():
    fields = {'PHIDP': np.ndarray,
              'WIDTH': np.ndarray,
              'DBZ': MaskedArray,
              'DBM': np.ndarray,
              'RHOHV': np.ndarray,
              'SNR': MaskedArray,
              'VEL': np.ndarray,
              'ZDR': MaskedArray}
    for field, field_type in fields.iteritems():
        description = "field : %s, type" % field
        check_field_type.description = description
        yield check_field_type, field, field_type


def check_field_type(field, field_type):
    assert type(radar.fields[field]['data']) is field_type


def test_field_first_points():
    # these values can be found using:
    # [round(radar.fields[f]['data'][0,0]) for f in radar.fields]
    fields = {'PHIDP': -111.0,
              'WIDTH': 5.0,
              'DBZ': 20.0,
              'DBM': -50.0,
              'RHOHV': 1.0,
              'SNR': 63.0,
              'VEL': -8.0,
              'ZDR': 0.0}
    for field, field_value in fields.iteritems():
        description = "field : %s, first point" % field
        check_field_first_point.description = description
        yield check_field_first_point, field, field_value


def check_field_first_point(field, value):
    assert round(radar.fields[field]['data'][0, 0]) == value
