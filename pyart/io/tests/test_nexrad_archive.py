""" Unit Tests for Py-ART's io/nexrad_archive.py module. """

import numpy as np
from numpy.ma.core import MaskedArray

import pyart

#######################################################
# read_nexrad_archive tests (verify radar attributes) #
#######################################################

# read in the sample file and create the four radar objects
radars = pyart.io.read_nexrad_archive(pyart.testing.NEXRAD_ARCHIVE_FILE)
radar = radars[0]


def test_refl_super_resolution():
    r = radars[0]
    assert 'reflectivity' in r.fields.keys()
    assert r.fields['reflectivity']['data'].shape == (2880, 1832)


def test_doppler_super_resolution():
    r = radars[1]
    assert 'velocity' in r.fields.keys()
    assert r.fields['velocity']['data'].shape == (1440, 1192)
    assert 'spectrum_width' in r.fields.keys()
    assert r.fields['spectrum_width']['data'].shape == (1440, 1192)
    assert 'correlation_coefficient' in r.fields.keys()
    assert r.fields['correlation_coefficient']['data'].shape == (1440, 1192)
    assert 'differential_reflectivity' in r.fields.keys()
    assert r.fields['differential_reflectivity']['data'].shape == (1440, 1192)
    assert 'differential_phase' in r.fields.keys()
    assert r.fields['differential_phase']['data'].shape == (1440, 1192)


def test_refl_standard_resolution():
    r = radars[2]
    assert 'reflectivity' in r.fields.keys()
    assert r.fields['reflectivity']['data'].shape == (4320, 1352)


def test_doppler_standard_resolution():
    r = radars[3]
    assert 'velocity' in r.fields.keys()
    assert r.fields['velocity']['data'].shape == (4320, 1192)
    assert 'spectrum_width' in r.fields.keys()
    assert r.fields['spectrum_width']['data'].shape == (4320, 1192)
    assert 'correlation_coefficient' in r.fields.keys()
    assert r.fields['correlation_coefficient']['data'].shape == (4320, 1192)
    assert 'differential_reflectivity' in r.fields.keys()
    assert r.fields['differential_reflectivity']['data'].shape == (4320, 1192)
    assert 'differential_phase' in r.fields.keys()
    assert r.fields['differential_phase']['data'].shape == (4320, 1192)


# time attribute
def test_time():
    assert 'comment' in radar.time.keys()
    assert 'long_name' in radar.time.keys()
    assert 'standard_name' in radar.time.keys()
    assert 'units' in radar.time.keys()
    assert 'calendar' in radar.time.keys()
    assert 'data' in radar.time.keys()
    assert radar.time['units'] == 'seconds since 2013-07-17T19:50:21Z'
    assert radar.time['data'].shape == (2880, )
    assert round(radar.time['data'][1]) == 1.


# range attribute
def test_range():
    assert 'long_name' in radar.range
    assert 'standard_name' in radar.range
    assert 'meters_to_center_of_first_gate' in radar.range
    assert 'meters_between_gates' in radar.range
    assert 'units' in radar.range
    assert 'data' in radar.range
    assert 'spacing_is_constant' in radar.range
    assert radar.range['data'].shape == (1832, )
    assert round(radar.range['data'][1]) == 2375.


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
    assert radar.latitude['data'].shape == (1, )
    assert round(radar.latitude['data']) == 48.0


# longitude attribute
def test_longitude():
    assert 'data' in radar.longitude
    assert 'standard_name' in radar.longitude
    assert 'units' in radar.longitude
    assert radar.longitude['data'].shape == (1, )
    assert round(radar.longitude['data']) == -122.0


# altitude attribute
def test_altitude():
    assert 'data' in radar.altitude
    assert 'standard_name' in radar.altitude
    assert 'units' in radar.altitude
    assert 'positive' in radar.altitude
    assert radar.altitude['data'].shape == (1, )
    assert round(radar.altitude['data']) == 161.0


# altitude_agl attribute
def test_altitude_agl():
    assert radar.altitude_agl is None


# sweep_number attribute
def test_sweep_number():
    assert 'standard_name' in radar.sweep_number
    assert np.all(radar.sweep_number['data'] == range(4))


# sweep_mode attribute
def test_sweep_mode():
    assert 'standard_name' in radar.sweep_mode
    assert radar.sweep_mode['data'].shape == (4, )
    assert np.all(radar.sweep_mode['data'] == ['azimuth_surveillance'])


# fixed_angle attribute
def test_fixed_angle():
    assert 'standard_name' in radar.fixed_angle
    assert 'units' in radar.fixed_angle
    assert radar.fixed_angle['data'].shape == (2880, )
    assert round(radar.fixed_angle['data'][0], 2) == 0.75


# sweep_start_ray_index attribute
def test_sweep_start_ray_index():
    assert 'long_name' in radar.sweep_start_ray_index
    assert radar.sweep_start_ray_index['data'].shape == (4, )
    assert round(radar.sweep_start_ray_index['data'][0]) == 0.0


# sweep_end_ray_index attribute
def test_sweep_end_ray_index():
    assert 'long_name' in radar.sweep_end_ray_index
    assert radar.sweep_end_ray_index['data'].shape == (4, )
    assert round(radar.sweep_end_ray_index['data'][0]) == 719.0


# target_scan_rate attribute
def test_target_scan_rate():
    assert radar.target_scan_rate is None


# azimuth attribute
def test_azimuth():
    assert 'standard_name' in radar.azimuth
    assert 'long_name' in radar.azimuth
    assert 'units' in radar.azimuth
    assert 'axis' in radar.azimuth
    assert round(radar.azimuth['data'][0]) == 350.0
    assert round(radar.azimuth['data'][10]) == 355.0


# elevation attribute
def test_elevation():
    assert 'standard_name' in radar.elevation
    assert 'long_name' in radar.azimuth
    assert 'units' in radar.elevation
    assert 'axis' in radar.elevation
    assert radar.elevation['data'].shape == (2880, )
    assert round(radar.elevation['data'][0], 2) == 0.75


# scan_rate attribute
def test_scan_rate():
    assert radar.scan_rate is None


# antenna_transition attribute
def test_antenna_transition():
    assert radar.antenna_transition is None


# instrument_parameters attribute
def test_instument_parameters():
    assert radar.instrument_parameters is None


# radar_calibration attribute
def test_radar_calibration():
    assert radar.radar_calibration is None


# ngates attribute
def test_ngates():
    assert radar.ngates == 1832


# nrays attribute
def test_nrays():
    assert radar.nrays == 2880


# nsweeps attribute
def test_nsweeps():
    assert radar.nsweeps == 4


####################
# fields attribute #
####################


def test_field_dics():
    fields = ['reflectivity', ]
    for field in fields:
        description = "field : %s, dictionary" % field
        check_field_dic.description = description
        yield check_field_dic, field


def check_field_dic(field):
    """ Check that the required keys are present in a field dictionary. """
    assert 'standard_name' in radar.fields[field]
    assert 'units' in radar.fields[field]
    assert '_FillValue' in radar.fields[field]
    assert 'coordinates' in radar.fields[field]


def test_field_shapes():
    fields = ['reflectivity', ]
    for field in fields:
        description = "field : %s, shape" % field
        check_field_shape.description = description
        yield check_field_shape, field


def check_field_shape(field):
    assert radar.fields[field]['data'].shape == (2880, 1832)


def test_field_types():
    fields = {'reflectivity': MaskedArray, }
    for field, field_type in fields.iteritems():
        description = "field : %s, type" % field
        check_field_type.description = description
        yield check_field_type, field, field_type


def check_field_type(field, field_type):
    assert type(radar.fields[field]['data']) is field_type


def test_field_first_points():
    # these values can be found using:
    # [round(radar.fields[f]['data'][0,0]) for f in radar.fields]
    fields = {'reflectivity': -32.0}
    for field, field_value in fields.iteritems():
        description = "field : %s, first point" % field
        check_field_first_point.description = description
        yield check_field_first_point, field, field_value


def check_field_first_point(field, value):
    assert round(radar.fields[field]['data'][0, 0]) == value


def test_compressed_archive():
    # the compressed archive only contains the first 120 radials
    cradars = pyart.io.read_nexrad_archive(
        pyart.testing.NEXRAD_ARCHIVE_COMPRESSED_FILE)
    r = cradars[0]
    assert 'reflectivity' in r.fields.keys()
    rdata = r.fields['reflectivity']['data']
    assert rdata.shape == (120, 1832)
    assert round(rdata[0, 0], 1) == 10.5
    assert round(rdata[2, 0], 1) == 9.5
    assert cradars[1] is None
    assert cradars[2] is None
    assert cradars[3] is None
