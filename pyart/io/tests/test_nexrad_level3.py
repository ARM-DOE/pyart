""" Unit Tests for Py-ART's io/nexrad_level3.py module. """

import numpy as np
from numpy.ma.core import MaskedArray

import pyart


def test_nexrad_level3_msg19():
    radar = pyart.io.read_nexrad_level3(pyart.testing.NEXRAD_LEVEL3_MSG19)

    assert radar.time['units'] == 'seconds since 2015-01-02T02:05:28Z'
    assert radar.time['data'].shape == (360, )
    assert round(radar.time['data'][0]) == 0.

    assert radar.range['data'].shape == (230, )
    assert round(radar.range['data'][100]) == 99900

    assert radar.scan_type == 'ppi'

    assert radar.latitude['data'].shape == (1, )
    assert round(radar.latitude['data'][0]) == 33.0

    assert radar.longitude['data'].shape == (1, )
    assert round(radar.longitude['data'][0]) == -87.0

    assert radar.altitude['data'].shape == (1, )
    assert round(radar.altitude['data'][0]) == 759.

    assert radar.altitude_agl is None

    assert radar.sweep_number['data'].shape == (1, )
    assert radar.sweep_number['data'][0] == 0

    assert radar.sweep_mode['data'].shape == (1, )
    assert radar.sweep_mode['data'].dtype.char == 'S'
    assert np.all(radar.sweep_mode['data'] == [b'azimuth_surveillance'])

    assert radar.fixed_angle['data'].shape == (1, )
    assert round(radar.fixed_angle['data'][0], 2) == 0.50

    assert radar.sweep_start_ray_index['data'].shape == (1, )
    assert round(radar.sweep_start_ray_index['data'][0]) == 0.0

    assert radar.sweep_end_ray_index['data'].shape == (1, )
    assert round(radar.sweep_end_ray_index['data'][0]) == 359.0

    assert radar.target_scan_rate is None

    assert round(radar.azimuth['data'][0]) == 320.0
    assert round(radar.azimuth['data'][10]) == 330.0

    assert radar.elevation['data'].shape == (360, )
    assert round(radar.elevation['data'][0], 2) == 0.50

    assert radar.scan_rate is None
    assert radar.antenna_transition is None
    assert radar.instrument_parameters is None
    assert radar.radar_calibration is None

    assert radar.ngates == 230
    assert radar.nrays == 360
    assert radar.nsweeps == 1

    assert 'reflectivity' in radar.fields.keys()
    assert radar.fields['reflectivity']['data'].shape == (360, 230)
    assert type(radar.fields['reflectivity']['data']) is MaskedArray
    assert round(radar.fields['reflectivity']['data'][10, 10]) == 25.


def test_nexrad_level3_msg161():
    radar = pyart.io.read_nexrad_level3(pyart.testing.NEXRAD_LEVEL3_MSG163)

    assert radar.time['units'] == 'seconds since 2015-01-02T02:05:28Z'
    assert radar.time['data'].shape == (360, )
    assert round(radar.time['data'][0]) == 0.

    assert radar.range['data'].shape == (1200, )
    assert round(radar.range['data'][100]) == 24975.

    assert radar.scan_type == 'ppi'

    assert radar.latitude['data'].shape == (1, )
    assert round(radar.latitude['data'][0]) == 33.0

    assert radar.longitude['data'].shape == (1, )
    assert round(radar.longitude['data'][0]) == -87.0

    assert radar.altitude['data'].shape == (1, )
    assert round(radar.altitude['data'][0]) == 759.

    assert radar.altitude_agl is None

    assert radar.sweep_number['data'].shape == (1, )
    assert radar.sweep_number['data'][0] == 0

    assert radar.sweep_mode['data'].shape == (1, )
    assert np.all(radar.sweep_mode['data'] == [b'azimuth_surveillance'])

    assert radar.fixed_angle['data'].shape == (1, )
    assert round(radar.fixed_angle['data'][0], 2) == 0.50

    assert radar.sweep_start_ray_index['data'].shape == (1, )
    assert round(radar.sweep_start_ray_index['data'][0]) == 0.0

    assert radar.sweep_end_ray_index['data'].shape == (1, )
    assert round(radar.sweep_end_ray_index['data'][0]) == 359.0

    assert radar.target_scan_rate is None

    assert round(radar.azimuth['data'][0]) == 329.0
    assert round(radar.azimuth['data'][10]) == 339.0

    assert radar.elevation['data'].shape == (360, )
    assert round(radar.elevation['data'][0], 2) == 0.50

    assert radar.scan_rate is None
    assert radar.antenna_transition is None
    assert radar.instrument_parameters is None
    assert radar.radar_calibration is None

    assert radar.ngates == 1200
    assert radar.nrays == 360
    assert radar.nsweeps == 1

    field_name = 'specific_differential_phase'
    assert field_name in radar.fields.keys()
    assert radar.fields[field_name]['data'].shape == (360, 1200)
    assert type(radar.fields[field_name]['data']) is MaskedArray
    assert round(radar.fields[field_name]['data'][103, 170]) == 2.


def test_nexrad_level3_msg161_fileobj():
    fh = open(pyart.testing.NEXRAD_LEVEL3_MSG163, 'rb')
    radar = pyart.io.read_nexrad_level3(fh)
    fh.close()
    field_name = 'specific_differential_phase'
    assert field_name in radar.fields.keys()
    assert radar.fields[field_name]['data'].shape == (360, 1200)
    assert type(radar.fields[field_name]['data']) is MaskedArray
    assert round(radar.fields[field_name]['data'][103, 170]) == 2.
