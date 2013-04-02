# tests of the Radar object in radar.py created from a sample nc file
# the following attributes are tested in some form or another for a
# nc radar object:
#   azimuth
#   elevation
#   fields
#   location
#   metadata
#   naz
#   nele
#   ngates
#   range
#   scan_type
#   sweep_info
#   time
#   inst_params

#   nsweeps     # not defined
#   sweep_mode  # not defined
#   sweep_number # not defined

# The following methods are not tested
# cf2rad                    : creator, move to nc_utils.py
# extract_rsl_pointing
# get_mdv_meta
# mdv2rad                   : creator, move to py_nc.py
# prtmode
# rsl2rad                   : creator, move to py4dd.py
# streamcf2rad              : creator, move to nc_utils.py
# ray_header_time_to_dict   :

from os.path import join, dirname

import numpy as np
from numpy.ma.core import MaskedArray
import netCDF4

import pyart

# read in the sample file and create a a Radar object
fname = join(dirname(__file__), 'sample_nc.nc')
nc_radar = pyart.io.read_netcdf(fname)
#nc_radar = pyart.io.radar.Radar(netCDF4.Dataset(fname))


# azimuth attribute
def test_nc_azimuth():
    """ NC radar, azimuth attribute """
    assert '_FillValue' in nc_radar.azimuth
    assert 'standard_name' in nc_radar.azimuth
    assert 'units' in nc_radar.azimuth
    assert nc_radar.azimuth['data'][0] == 6.0
    assert nc_radar.azimuth['data'][10] == 16.0


# elevation attribute
def test_nc_elevation():
    """ NC radar, elevation attribute """
    assert '_FillValue' in nc_radar.elevation
    assert 'positive' in nc_radar.elevation
    assert 'standard_name' in nc_radar.elevation
    assert 'units' in nc_radar.elevation
    assert nc_radar.elevation['data'].shape == (3524, )
    assert round(nc_radar.elevation['data'][0], 2) == 6.4


# fields attribute
def test_nc_fields_dictionaries():
    field_names = [u'PHIDP', u'WIDTH', u'DBZ', u'DBM', u'RHOHV', u'SNR',
                   u'VEL', u'ZDR']
    for field_name in field_names:
        description = "NC radar, fields : %s, required keys" % field_name
        check_required_field_keys.description = description
        yield check_required_field_keys, field_name


def check_required_field_keys(field_name):
    """ Check that the 11 standard keys are present in a field """
    assert u'threshold_value' in nc_radar.fields[field_name]
    assert u'_FillValue' in nc_radar.fields[field_name]
    assert u'grid_mapping' in nc_radar.fields[field_name]
    assert u'scale_factor' in nc_radar.fields[field_name]
    assert u'sampling_ratio' in nc_radar.fields[field_name]
    assert u'coordinates' in nc_radar.fields[field_name]
    assert u'add_offset' in nc_radar.fields[field_name]
    assert u'standard_name' in nc_radar.fields[field_name]
    assert u'threshold_field_name' in nc_radar.fields[field_name]
    assert u'units' in nc_radar.fields[field_name]
    assert 'data' in nc_radar.fields[field_name]


def test_nc_fields_data_shape():
    field_names = [u'PHIDP', u'WIDTH', u'DBZ', u'DBM', u'RHOHV', u'SNR',
                   u'VEL', u'ZDR']
    for field_name in field_names:
        description = "NC radar, fields : %s, data shape" % field_name
        check_field_data_shape.description = description
        yield check_field_data_shape, field_name


def check_field_data_shape(field_name):
    assert nc_radar.fields[field_name]['data'].shape == (3524, 1333)


def test_nc_fields_data_type():
    field_names = {'PHIDP': np.ndarray,
                   'WIDTH': np.ndarray,
                   'DBZ': MaskedArray,
                   'DBM': np.ndarray,
                   'RHOHV': np.ndarray,
                   'SNR': MaskedArray,
                   'VEL': np.ndarray,
                   'ZDR': MaskedArray}
    for field_name, field_type in field_names.iteritems():
        description = "NC radar, fields : %s, data type" % field_name
        check_field_data_type.description = description
        yield check_field_data_type, field_name, field_type


def check_field_data_type(field_name, field_type):
    assert type(nc_radar.fields[field_name]['data']) is field_type


def test_nc_fields_data_first_points():
    """ NC radar, fields : first point in data """
    def first_point(field_name):
        return round(nc_radar.fields[field_name]['data'][0, 0])
    # these values can be found using:
    # [round(nc_radar.fields[f]['data'][0,0]) for f in nc_radar.fields]
    assert first_point('PHIDP') == -111.0
    assert first_point('WIDTH') == 5.0
    assert first_point('DBZ') == 20.0
    assert first_point('DBM') == -50.0
    assert first_point('RHOHV') == 1.0
    assert first_point('SNR') == 63.0
    assert first_point('VEL') == -8.0
    assert first_point('ZDR') == 0.0


# location attribute
def test_nc_location():
    elements = ['latitude', 'altitude', 'longitude']
    for element in elements:
        description = "NC radar, location : %s " % element
        check_location_element.description = description
        yield check_location_element, element


def check_location_element(element):
    """ Check that location attributes dictionaries have all required keys. """
    assert 'data' in nc_radar.location[element]
    assert 'standard_name' in nc_radar.location[element]
    assert 'units' in nc_radar.location[element]


def test_nc_location_data():
    """ NC radar, location : elevation extra key """
    assert 'positive' in nc_radar.location['altitude']


def test_nc_location_data():
    """ NC radar, location : data """
    assert round(nc_radar.location['latitude']['data']) == -2.0
    assert round(nc_radar.location['longitude']['data']) == 147.0
    assert round(nc_radar.location['altitude']['data']) == 126.0


# metadata attribute
def test_nc_metadata():
    """ NC radar, metadata attribute """
    assert 'comment' in nc_radar.metadata
    assert 'scan_id' in nc_radar.metadata
    assert 'platform_is_mobile' in nc_radar.metadata
    assert 'site_name' in nc_radar.metadata
    assert 'title' in nc_radar.metadata
    assert 'scan_name' in nc_radar.metadata
    assert 'n_gates_vary' in nc_radar.metadata
    assert 'Conventions' in nc_radar.metadata
    assert 'source' in nc_radar.metadata
    assert 'version' in nc_radar.metadata
    assert 'references' in nc_radar.metadata
    assert 'instrument_name' in nc_radar.metadata
    assert 'ray_times_increase' in nc_radar.metadata
    assert 'institution' in nc_radar.metadata
    assert 'history' in nc_radar.metadata


# naz attribute
def test_nc_naz():
    """ NC radar, naz attribute """
    assert nc_radar.naz == 271


# nele attribute
def test_nc_nele():
    """ NC radar, nele attribute """
    assert nc_radar.nele == 10


# ngates attribute
def test_nc_ngates():
    """ NC radar, ngates attribute """
    assert nc_radar.ngates == 1333


# nsweeps attribute
#def test_nc_nsweeps():
#    """ NC radar, nsweeps attribute """
#    assert nc_radar.nsweeps == 17


# range attribute
def test_nc_range():
    """ NC radar, range attribute """
    assert 'long_name' in nc_radar.range
    assert 'standard_name' in nc_radar.range
    assert 'meters_to_center_of_first_gate' in nc_radar.range
    assert 'meters_between_gates' in nc_radar.range
    assert 'units' in nc_radar.range
    assert 'data' in nc_radar.range
    assert 'spacing_is_constant' in nc_radar.range
    assert nc_radar.range['data'].shape == (1333, )
    assert round(nc_radar.range['data'][0]) == 0.0


# scan_type attribute
def test_nc_scan_type():
    """ NC radar, scan type attribute """
    assert nc_radar.scan_type == 'ppi'


# sweep_info attribute
def test_nc_sweep_info():
    """ NC radar, sweep_info attribute """
    assert 'sweep_start_ray_index' in nc_radar.sweep_info.keys()
    assert 'sweep_mode' in nc_radar.sweep_info.keys()
    assert 'sweep_number' in nc_radar.sweep_info.keys()
    assert 'sweep_end_ray_index' in nc_radar.sweep_info.keys()
    assert 'fixed_angle' in nc_radar.sweep_info.keys()


def test_nc_sweep_info_start():
    """ NC radar, sweep_info : sweep_start_ray_index """
    ssri = nc_radar.sweep_info['sweep_start_ray_index']
    assert 'data' in ssri.keys()
    assert 'standard_name' in ssri.keys()
    assert '_FillValue' in ssri.keys()
    assert ssri['data'].shape == (10, )


def test_nc_sweep_info_mode():
    """ NC radar, sweep_info : sweep_mode """
    sm = nc_radar.sweep_info['sweep_mode']
    assert 'options' in sm.keys()
    assert 'standard_name' in sm.keys()
    assert 'data' in sm.keys()
    assert sm['data'].shape == (10, 32)


def test_nc_sweep_info_number():
    """ NC radar, sweep_info : sweep_number """
    sn = nc_radar.sweep_info['sweep_number']
    assert 'data' in sn.keys()
    assert 'standard_name' in sn.keys()
    assert '_FillValue' in sn.keys()
    assert np.all(sn['data'] == range(9, 19))


def test_nc_sweep_info_end():
    """ NC radar, sweep_info : sweep_end_ray_index """
    seri = nc_radar.sweep_info['sweep_end_ray_index']
    assert 'data' in seri.keys()
    assert 'standard_name' in seri.keys()
    assert '_FillValue' in seri.keys()
    assert seri['data'].shape == (10, )


def test_nc_sweep_info_angle():
    """ NC radar, sweep_info : fixed_angle """
    fa = nc_radar.sweep_info['fixed_angle']
    assert 'data' in fa.keys()
    assert 'standard_name' in fa.keys()
    assert 'units' in fa.keys()
    assert '_FillValue' in fa.keys()
    assert fa['data'].shape == (10, )


# sweep mode attribute
#def test_nc_sweep_mode():
#    """ NC radar, sweep_mode attribute """
#    assert np.all(nc_radar.sweep_mode == ['ppi'] * 17)


# sweep_number attribute
#def test_nc_sweep_number():
#    """ NC radar, sweep_number attribute """
#    assert np.all(nc_radar.sweep_number == range(17))


# time attribute
def test_nc_time():
    """ NC radar, time attribute """
    assert 'comment' in nc_radar.time.keys()
    assert 'long_name' in nc_radar.time.keys()
    assert 'standard_name' in nc_radar.time.keys()
    assert 'units' in nc_radar.time.keys()
    #assert 'calendar' in nc_radar.time.keys()
    assert 'data' in nc_radar.time.keys()
    assert nc_radar.time['units'] == 'seconds since 2012-12-11T22:48:21Z'
    assert nc_radar.time['data'].shape == (3524, )
    assert round(nc_radar.time['data'][600]) == 51.


# inst_params attribute
def test_nc_inst_params():
    """ NC radar, inst_params attribute """
    assert 'prt_ratio' in nc_radar.inst_params.keys()
    assert 'follow_mode' in nc_radar.inst_params.keys()
    assert 'frequency' in nc_radar.inst_params.keys()
    assert 'prt_mode' in nc_radar.inst_params.keys()
    assert 'prt' in nc_radar.inst_params.keys()
    assert 'pulse_width' in nc_radar.inst_params.keys()
    assert 'polarization_mode' in nc_radar.inst_params.keys()
    assert 'unambiguous_range' in nc_radar.inst_params.keys()
    assert 'nyquist_velocity' in nc_radar.inst_params.keys()
    assert 'n_samples' in nc_radar.inst_params.keys()


def test_nc_inst_params_prt_ratio():
    """ NC radar, inst_params : ptr_ratio """
    pr = nc_radar.inst_params['prt_ratio']
    assert '_FillValue' in pr
    assert 'data' in pr
    assert 'meta_group' in pr
    assert 'standard_name' in pr
    assert 'units' in pr
    assert np.all(pr['data'] == 1)


def test_nc_inst_params_follow_mode():
    """ NC radar, inst_params : follow_mode """
    fm = nc_radar.inst_params['follow_mode']
    assert 'data' in fm
    assert 'meta_group' in fm
    assert 'standard_name' in fm
    assert 'options' in fm
    assert fm['data'].shape == (10, 32)


def test_nc_inst_params_frequency():
    """ NC radar, inst_params : frequency """
    f = nc_radar.inst_params['frequency']
    assert '_FillValue' in f
    assert 'data' in f
    assert 'meta_group' in f
    assert 'standard_name' in f
    assert 'units' in f
    assert f['data'].shape == (50, )


def test_nc_inst_params_prt_mode():
    """ NC radar, inst_params : prt_mode """
    pm = nc_radar.inst_params['prt_mode']
    assert 'data' in pm
    assert 'meta_group' in pm
    assert 'standard_name' in pm
    assert 'options' in pm
    assert pm['data'].shape == (10, 32)


def test_nc_inst_params_prt():
    """ NC radar, inst_params : prt """
    p = nc_radar.inst_params['prt']
    assert '_FillValue' in p
    assert 'data' in p
    assert 'meta_group' in p
    assert 'standard_name' in p
    assert 'units' in p
    assert np.all(p['data'] == 0.0016665)


def test_nc_inst_params_pulse_width():
    """ NC radar, inst_params : pulse_width """
    pw = nc_radar.inst_params['pulse_width']
    assert '_FillValue' in pw
    assert 'data' in pw
    assert 'meta_group' in pw
    assert 'standard_name' in pw
    assert 'units' in pw
    assert np.all(pw['data'] == 1.99999999e-06)


def test_nc_inst_params_polarization_mode():
    """ NC radar, inst_params : polarization_mode """
    pm = nc_radar.inst_params['prt_mode']
    assert 'data' in pm
    assert 'meta_group' in pm
    assert 'standard_name' in pm
    assert 'options' in pm
    assert pm['data'].shape == (10, 32)


def test_nc_inst_params_unambiguous_range():
    """ NC radar, inst_params : unambiguous_range """
    ur = nc_radar.inst_params['unambiguous_range']
    assert '_FillValue' in ur
    assert 'data' in ur
    assert 'meta_group' in ur
    assert 'standard_name' in ur
    assert 'units' in ur
    assert np.all(ur['data'] == 249802.0625)


def test_nc_inst_params_nyquist_velocity():
    """ NC radar, inst_params : nyquist_velocity """
    nv = nc_radar.inst_params['nyquist_velocity']
    assert '_FillValue' in nv
    assert 'data' in nv
    assert 'meta_group' in nv
    assert 'standard_name' in nv
    assert 'units' in nv
    assert np.all(np.round(nv['data']) == 16)


def test_nc_inst_params_n_samples():
    """ NC radar, inst_params : n_samples """
    ns = nc_radar.inst_params['n_samples']
    assert '_FillValue' in ns
    assert 'data' in ns
    assert 'meta_group' in ns
    assert 'standard_name' in ns
    assert np.all(ns['data'] == 128)
