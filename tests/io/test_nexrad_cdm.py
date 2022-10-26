""" Unit Tests for Py-ART's io/nexrad_cdm.py module. """

import bz2

import numpy as np
from numpy.testing import assert_almost_equal
from numpy.ma.core import MaskedArray
import pytest

import pyart

###################################################
# read_nexrad_cdm tests (verify radar attributes) #
###################################################

# read in the sample file and create the radar objects
# We need to decompress the bz2 file which contains this data
with pyart.testing.InTemporaryDirectory():
    tmpfile = 'tmp_nexrad.nc'
    with open(tmpfile, 'wb') as f:
        f.write(bz2.BZ2File(pyart.testing.NEXRAD_CDM_FILE).read())
    radar = pyart.io.read_nexrad_cdm(tmpfile)


# time attribute
def test_time():
    assert 'comment' in radar.time.keys()
    assert 'long_name' in radar.time.keys()
    assert 'standard_name' in radar.time.keys()
    assert 'units' in radar.time.keys()
    assert 'calendar' in radar.time.keys()
    assert 'data' in radar.time.keys()
    assert radar.time['units'] == 'seconds since 2013-07-17T19:50:21Z'
    assert radar.time['data'].shape == (7200, )
    assert_almost_equal(radar.time['data'][1], 0.677, 3)


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
    assert_almost_equal(radar.range['data'][1], 2375, 0)


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
    assert_almost_equal(radar.latitude['data'], 48, 0)


# longitude attribute
def test_longitude():
    assert 'data' in radar.longitude
    assert 'standard_name' in radar.longitude
    assert 'units' in radar.longitude
    assert radar.longitude['data'].shape == (1, )
    assert_almost_equal(radar.longitude['data'], -122, 0)


# altitude attribute
def test_altitude():
    assert 'data' in radar.altitude
    assert 'standard_name' in radar.altitude
    assert 'units' in radar.altitude
    assert 'positive' in radar.altitude
    assert radar.altitude['data'].shape == (1, )
    assert_almost_equal(radar.altitude['data'], 151, 0)     # 10 m different


# altitude_agl attribute
def test_altitude_agl():
    assert radar.altitude_agl is None


# sweep_number attribute
def test_sweep_number():
    assert 'standard_name' in radar.sweep_number
    assert np.all(radar.sweep_number['data'] == range(16))


# sweep_mode attribute
def test_sweep_mode():
    assert 'standard_name' in radar.sweep_mode
    assert radar.sweep_mode['data'].shape == (16, )
    assert radar.sweep_mode['data'].dtype.char == 'S'
    assert np.all(radar.sweep_mode['data'] == [b'azimuth_surveillance'])


# fixed_angle attribute
def test_fixed_angle():
    assert 'standard_name' in radar.fixed_angle
    assert 'units' in radar.fixed_angle
    assert radar.fixed_angle['data'].shape == (16, )
    assert_almost_equal(radar.fixed_angle['data'][0], 0.53, 2)


# sweep_start_ray_index attribute
def test_sweep_start_ray_index():
    assert 'long_name' in radar.sweep_start_ray_index
    assert radar.sweep_start_ray_index['data'].shape == (16, )
    assert_almost_equal(radar.sweep_start_ray_index['data'][0], 0, 0)


# sweep_end_ray_index attribute
def test_sweep_end_ray_index():
    assert 'long_name' in radar.sweep_end_ray_index
    assert radar.sweep_end_ray_index['data'].shape == (16, )
    assert_almost_equal(radar.sweep_end_ray_index['data'][0], 719, 0)


# target_scan_rate attribute
def test_target_scan_rate():
    assert radar.target_scan_rate is None


# azimuth attribute
def test_azimuth():
    assert 'standard_name' in radar.azimuth
    assert 'long_name' in radar.azimuth
    assert 'units' in radar.azimuth
    assert 'axis' in radar.azimuth
    assert_almost_equal(radar.azimuth['data'][0], 350, 0)
    assert_almost_equal(radar.azimuth['data'][10], 355, 0)


# elevation attribute
def test_elevation():
    assert 'standard_name' in radar.elevation
    assert 'long_name' in radar.azimuth
    assert 'units' in radar.elevation
    assert 'axis' in radar.elevation
    assert radar.elevation['data'].shape == (7200, )
    assert_almost_equal(radar.elevation['data'][0], 0.75, 2)


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
    assert radar.nrays == 7200


# nsweeps attribute
def test_nsweeps():
    assert radar.nsweeps == 16


####################
# fields attribute #
####################

fields = [
    'differential_phase', 'spectrum_width', 'cross_correlation_ratio',
    'reflectivity', 'differential_reflectivity', 'velocity']


@pytest.mark.parametrize("field", fields)
def test_field_dics(field):
    description = "field : %s, dictionary" % field
    check_field_dic.description = description
    check_field_dic(field)


def check_field_dic(field):
    """ Check that the required keys are present in a field dictionary. """
    assert 'standard_name' in radar.fields[field]
    assert 'units' in radar.fields[field]
    assert '_FillValue' in radar.fields[field]
    assert 'coordinates' in radar.fields[field]


@pytest.mark.parametrize("field", fields)
def test_field_shapes(field):
    description = "field : %s, shape" % field
    check_field_shape.description = description
    check_field_shape(field)


def check_field_shape(field):
    assert radar.fields[field]['data'].shape == (7200, 1832)


fields = {'differential_phase': MaskedArray,
          'spectrum_width': MaskedArray,
          'cross_correlation_ratio': MaskedArray,
          'reflectivity': MaskedArray,
          'differential_reflectivity': MaskedArray,
          'velocity': MaskedArray}
@pytest.mark.parametrize(
    "field, field_type", fields.items(), ids=list(fields.keys()))
def test_field_types(field, field_type):
    description = "field : %s, type" % field
    check_field_type.description = description
    check_field_type(field, field_type)


def check_field_type(field, field_type):
    assert type(radar.fields[field]['data']) is field_type


fields = {'differential_phase': 181.0,
          'spectrum_width': np.ma.masked,
          'cross_correlation_ratio': 0.0,
          'reflectivity': -32.0,
          'differential_reflectivity': -8.0,
          'velocity': np.ma.masked}
@pytest.mark.parametrize(
    "field, field_value", fields.items(), ids=list(fields.keys()))
def test_field_first_points(field, field_value):
    # these values can be found using:
    # [round(radar.fields[f]['data'][0,0]) for f in radar.fields]
    description = "field : %s, first point" % field
    check_field_first_point.description = description
    check_field_first_point(field, field_value)


def check_field_first_point(field, value):
    if np.ma.is_masked(value):
        assert np.ma.is_masked(radar.fields[field]['data'][0, 0])
    else:
        assert_almost_equal(radar.fields[field]['data'][0, 0], value, 0)
