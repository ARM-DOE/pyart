""" Unit Tests for Py-ART's io/nexrad_archive.py module using a MSG1 file. """

import warnings

import numpy as np
from numpy.testing import assert_almost_equal
from numpy.ma.core import MaskedArray
import pytest

import pyart

#######################################################
# read_nexrad_archive tests (verify radar attributes) #
#######################################################


NRAYS = 2567
NGATES = 1840


# read in the sample file, ignore warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', category=UserWarning)
    radar = pyart.io.read_nexrad_archive(
        pyart.testing.NEXRAD_ARCHIVE_MSG1_FILE, station='KLOT')


# time attribute
def test_time():
    assert 'comment' in radar.time.keys()
    assert 'long_name' in radar.time.keys()
    assert 'standard_name' in radar.time.keys()
    assert 'units' in radar.time.keys()
    assert 'calendar' in radar.time.keys()
    assert 'data' in radar.time.keys()
    assert radar.time['units'] == 'seconds since 2003-01-01T00:09:21Z'
    assert radar.time['data'].shape == (NRAYS, )
    assert_almost_equal(radar.time['data'][1], .506, 3)


# range attribute
def test_range():
    assert 'long_name' in radar.range
    assert 'standard_name' in radar.range
    assert 'meters_to_center_of_first_gate' in radar.range
    assert 'meters_between_gates' in radar.range
    assert 'units' in radar.range
    assert 'data' in radar.range
    assert 'spacing_is_constant' in radar.range
    assert radar.range['data'].shape == (NGATES, )
    assert_almost_equal(radar.range['data'][0], -375, 0)
    assert_almost_equal(radar.range['data'][1], -125, 0)
    assert_almost_equal(radar.range['data'][-1], 459375, 0)


# metadata attribute
def test_metadata():
    assert 'instrument_name' in radar.metadata
    assert 'source' in radar.metadata
    assert 'vcp_pattern' not in radar.metadata


# scan_type attribute
def test_scan_type():
    assert radar.scan_type == 'ppi'


# latitude attribute
def test_latitude():
    assert 'data' in radar.latitude
    assert 'standard_name' in radar.latitude
    assert 'units' in radar.latitude
    assert radar.latitude['data'].shape == (1, )
    assert_almost_equal(radar.latitude['data'], 41.6, 1)


# longitude attribute
def test_longitude():
    assert 'data' in radar.longitude
    assert 'standard_name' in radar.longitude
    assert 'units' in radar.longitude
    assert radar.longitude['data'].shape == (1, )
    assert_almost_equal(radar.longitude['data'], -88.1, 1)


# altitude attribute
def test_altitude():
    assert 'data' in radar.altitude
    assert 'standard_name' in radar.altitude
    assert 'units' in radar.altitude
    assert 'positive' in radar.altitude
    assert radar.altitude['data'].shape == (1, )
    assert_almost_equal(radar.altitude['data'], 663, 0)


# altitude_agl attribute
def test_altitude_agl():
    assert radar.altitude_agl is None


# sweep_number attribute
def test_sweep_number():
    assert 'standard_name' in radar.sweep_number
    assert np.all(radar.sweep_number['data'] == range(7))


# sweep_mode attribute
def test_sweep_mode():
    assert 'standard_name' in radar.sweep_mode
    assert radar.sweep_mode['data'].shape == (7, )
    assert np.all(radar.sweep_mode['data'] == [b'azimuth_surveillance'])


# fixed_angle attribute
def test_fixed_angle():
    assert 'standard_name' in radar.fixed_angle
    assert 'units' in radar.fixed_angle
    assert radar.fixed_angle['data'].shape == (7, )
    assert_almost_equal(radar.fixed_angle['data'][0], 0.5, 1)


# sweep_start_ray_index attribute
def test_sweep_start_ray_index():
    assert 'long_name' in radar.sweep_start_ray_index
    assert radar.sweep_start_ray_index['data'].shape == (7, )
    assert_almost_equal(radar.sweep_start_ray_index['data'][1], 367, 0)


# sweep_end_ray_index attribute
def test_sweep_end_ray_index():
    assert 'long_name' in radar.sweep_end_ray_index
    assert radar.sweep_end_ray_index['data'].shape == (7, )
    assert_almost_equal(radar.sweep_end_ray_index['data'][1], 733, 0)


# target_scan_rate attribute
def test_target_scan_rate():
    assert radar.target_scan_rate is None


# azimuth attribute
def test_azimuth():
    assert 'standard_name' in radar.azimuth
    assert 'long_name' in radar.azimuth
    assert 'units' in radar.azimuth
    assert 'axis' in radar.azimuth
    assert_almost_equal(radar.azimuth['data'][0], 245.874, 3)
    assert_almost_equal(radar.azimuth['data'][10], 255.71777, 3)


# elevation attribute
def test_elevation():
    assert 'standard_name' in radar.elevation
    assert 'long_name' in radar.azimuth
    assert 'units' in radar.elevation
    assert 'axis' in radar.elevation
    assert radar.elevation['data'].shape == (NRAYS, )
    assert_almost_equal(radar.elevation['data'][0], 0.48, 2)


# scan_rate attribute
def test_scan_rate():
    assert radar.scan_rate is None


# antenna_transition attribute
def test_antenna_transition():
    assert radar.antenna_transition is None


# instrument_parameters attribute
def test_instument_parameters():
    assert 'nyquist_velocity' in radar.instrument_parameters
    nyq = radar.instrument_parameters['nyquist_velocity']['data']
    assert_almost_equal(nyq[0], 0, 0)
    assert_almost_equal(nyq[400], 28.34, 2)

    assert 'unambiguous_range' in radar.instrument_parameters
    unamb = radar.instrument_parameters['unambiguous_range']['data']
    assert_almost_equal(unamb[0], 466000., 0)
    assert_almost_equal(unamb[400], 137000., 0)


# radar_calibration attribute
def test_radar_calibration():
    assert radar.radar_calibration is None


# ngates attribute
def test_ngates():
    assert radar.ngates == NGATES


# nrays attribute
def test_nrays():
    assert radar.nrays == NRAYS


# nsweeps attribute
def test_nsweeps():
    assert radar.nsweeps == 7


####################
# fields attribute #
####################


FIELDS = ['reflectivity', 'spectrum_width', 'velocity']


@pytest.mark.parametrize("field", FIELDS)
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


@pytest.mark.parametrize("field", FIELDS)
def test_field_shapes(field):
    description = "field : %s, shape" % field
    check_field_shape.description = description
    check_field_shape(field)


def check_field_shape(field):
    assert radar.fields[field]['data'].shape == (NRAYS, NGATES)


fields = {'spectrum_width': MaskedArray,
          'reflectivity': MaskedArray,
          'velocity': MaskedArray}
@pytest.mark.parametrize(
    "field, field_type", fields.items(), ids=list(fields.keys()))
def test_field_types(field, field_type):
    description = "field : %s, type" % field
    check_field_type.description = description
    check_field_type(field, field_type)


def check_field_type(field, field_type):
    assert type(radar.fields[field]['data']) is field_type


def test_field_data_spectrum_width():
    assert_almost_equal(radar.fields['spectrum_width']['data'][367, 16:19],
                        [0.0, 15.5, 16.5])


def test_field_data_velocity():
    assert_almost_equal(radar.fields['velocity']['data'][367, 16:19],
                        [-17.5, -13.0, 7.5])


def test_field_data_reflectivity():
    assert_almost_equal(
        radar.fields['reflectivity']['data'][0, 8:16],
        [1.0, 1.0, 0.4375, -0.6875, -1.8125, -2.9375, -3.5, -3.5])


def test_field_data_nearest_neighbor():
    # read in the sample file, ignore warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=UserWarning)
        radar = pyart.io.read_nexrad_archive(
            pyart.testing.NEXRAD_ARCHIVE_MSG1_FILE,
            station='KLOT', linear_interp=False)
    assert_almost_equal(
        radar.fields['reflectivity']['data'][0, 8:16],
        [1.0, 1.0, 1.0, 1.0, -3.5, -3.5, -3.5, -3.5])


def test_single_scan():
    radar = pyart.io.read_nexrad_archive(
        pyart.testing.NEXRAD_ARCHIVE_MSG1_FILE,
        station='KLOT', scans=[0])
    assert radar.nrays == 367
    assert radar.ngates == 460
    assert radar.nsweeps == 1
    assert_almost_equal(radar.fields['reflectivity']['data'][0, 6:10],
                        [19, 26.5, 14.5, 21.5])
    assert_almost_equal(radar.range['data'][0], 0, 0)
    assert_almost_equal(radar.range['data'][1], 1000, 0)


def test_only_ref():
    radar = pyart.io.read_nexrad_archive(
        pyart.testing.NEXRAD_ARCHIVE_MSG1_FILE,
        station='KLOT', scans=[0, 2, 4, 5, 6],
        exclude_fields=['velocity', 'spectrum_width'])
    assert_almost_equal(radar.range['data'][0], 0, 0)
    assert_almost_equal(radar.range['data'][1], 1000, 0)
    assert_almost_equal(radar.fields['reflectivity']['data'][0, 6:10],
                        [19, 26.5, 14.5, 21.5])
    assert radar.nrays == 367 + 368 + 366 + 366 + 366
    assert radar.ngates == 460
    assert radar.nsweeps == 5
