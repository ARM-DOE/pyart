""" Unit Tests for Py-ART's io/radar.py module. """

import sys
from io import BytesIO

import numpy as np
from numpy.testing import assert_raises
import pyart


def test_extract_sweeps():
    radar = pyart.testing.make_empty_ppi_radar(100, 360, 3)
    radar.fields['reflectivity'] = {'data': np.zeros((1080, 100))}
    radar.fields['velocity'] = {'data': np.zeros((1080, 100))}

    eradar = radar.extract_sweeps([0, 2])

    # extracted radar should have 720 rays, 2 sweeps, 100 gates
    assert eradar.time['data'].shape == (720, )
    assert eradar.range['data'].shape == (100, )

    assert eradar.metadata['instrument_name'] == 'fake_radar'
    assert eradar.scan_type == 'ppi'

    assert eradar.latitude['data'].shape == (1, )
    assert eradar.longitude['data'].shape == (1, )
    assert eradar.altitude['data'].shape == (1, )
    assert eradar.altitude_agl is None

    assert eradar.sweep_number['data'].shape == (2, )
    assert eradar.sweep_mode['data'].shape == (2, )
    assert eradar.fixed_angle['data'].shape == (2, )
    assert eradar.sweep_start_ray_index['data'].shape == (2, )
    assert eradar.sweep_end_ray_index['data'].shape == (2, )
    assert eradar.target_scan_rate is None

    assert eradar.azimuth['data'].shape == (720, )
    assert eradar.elevation['data'].shape == (720, )
    assert eradar.scan_rate is None
    assert eradar.antenna_transition is None

    assert eradar.instrument_parameters is None
    assert eradar.radar_calibration is None

    assert eradar.ngates == 100
    assert eradar.nrays == 720
    assert eradar.nsweeps == 2

    assert eradar.fields['reflectivity']['data'].shape == (720, 100)
    assert eradar.fields['velocity']['data'].shape == (720, 100)


def test_extract_sweeps_extra():
    radar = pyart.testing.make_empty_ppi_radar(10, 36, 3)
    radar.instrument_parameters = {
        'prt': {'data': np.zeros((108, ))},
        'prt_mode': {'data': np.array(['fixed'] * 3)},
        'radar_antenna_gain_h': {'data': np.array(0)},
    }

    radar.radar_calibration = {
        'r_calib_index': {'data': np.zeros((108, ))},
        'r_calib_time': {'data': np.zeros((8, ))}
    }

    eradar = radar.extract_sweeps([0, 2])

    instr = eradar.instrument_parameters
    assert instr['prt']['data'].shape == (72, )
    assert instr['prt_mode']['data'].shape == (2, )
    assert instr['radar_antenna_gain_h']['data'].shape == ()

    calib = eradar.radar_calibration
    assert calib['r_calib_index']['data'].shape == (72, )
    assert calib['r_calib_time']['data'].shape == (8, )


def test_extract_sweeps_errors():
    radar = pyart.testing.make_empty_ppi_radar(10, 36, 2)
    assert_raises(ValueError, radar.extract_sweeps, [0, 2])
    assert_raises(ValueError, radar.extract_sweeps, [-1, 1])


def test_radar_creation():
    radar = pyart.testing.make_target_radar()
    assert isinstance(radar, pyart.io.Radar)


def test_add_field():
    radar = pyart.testing.make_target_radar()
    dic = {'data': np.zeros((360, 50)), 'standard_name': 'test'}
    radar.add_field('test', dic)
    assert 'test' in radar.fields
    assert 'data' in radar.fields['test']
    assert radar.fields['test']['standard_name'] == 'test'


def test_add_field_errors():
    radar = pyart.testing.make_target_radar()

    assert_raises(ValueError, radar.add_field, 'reflectivity', {})

    dic = {'dat': np.zeros((360, 50)), 'standard_name': 'test'}
    assert_raises(KeyError, radar.add_field, 'test', dic)

    dic = {'data': np.zeros((360, 49)), 'standard_name': 'test'}
    assert_raises(ValueError, radar.add_field, 'test', dic)


def test_add_field_like():
    radar = pyart.testing.make_target_radar()
    data = np.zeros((360, 50))
    radar.add_field_like('reflectivity', 'test', data)
    assert 'test' in radar.fields
    assert 'data' in radar.fields['test']
    assert radar.fields['test']['units'] == 'dBZ'


def test_add_field_like_errors():
    radar = pyart.testing.make_target_radar()
    assert_raises(ValueError, radar.add_field_like, 'foo', 'bar', [])


def test_info_levels():
    for level in ['standard', 's', 'compact', 'c', 'full', 'f']:
        yield check_info, level


def test_info_nonstandard():
    # build a non-standard radar object for testing all paths in info
    radar = pyart.testing.make_target_radar()
    radar.fields['reflectivity']['data'] = [1, 2, 3, 4]
    radar.instrument_parameters = {'foobar': {'data': [1, 2], 'bar': 'foo'}}
    radar.radar_calibration = {'foobar': {'data': [1, 2], 'bar': 'foo'}}
    check_info('standard', radar)


def check_info(level, radar=None):
    out = BytesIO()
    get_info(level, out, radar)
    # don't check the output, just that something was printed.
    assert len(out.getvalue()) != 0


def get_info(level='standard', out=sys.stdout, radar=None):
    if radar is None:
        radar = pyart.testing.make_target_radar()
    radar.info(level, out)


def test_info_errors():
    assert_raises(ValueError, check_info, 'foo')


def test_is_vpt():
    radar = pyart.testing.make_empty_ppi_radar(10, 36, 3)
    assert not pyart.io.is_vpt(radar)
    pyart.io.to_vpt(radar)
    assert pyart.io.is_vpt(radar)


def test_to_vpt():
    # single scan
    radar = pyart.testing.make_empty_ppi_radar(10, 36, 3)
    radar.instrument_parameters = {
        'prt_mode': {'data': np.array(['fixed'] * 3)}
    }
    pyart.io.to_vpt(radar)
    assert pyart.io.is_vpt(radar)
    assert radar.nsweeps == 1
    assert radar.azimuth['data'][10] == 0.0
    assert radar.elevation['data'][0] == 90.0
    assert len(radar.instrument_parameters['prt_mode']['data']) == 1

    # multiple scans
    radar = pyart.testing.make_empty_ppi_radar(10, 36, 3)
    radar.instrument_parameters = {
        'prt_mode': {'data': np.array(['fixed'] * 3)}
    }
    pyart.io.to_vpt(radar, False)
    assert pyart.io.is_vpt(radar)
    assert radar.nsweeps == 108
    assert radar.azimuth['data'][10] == 10.0
    assert radar.elevation['data'][0] == 90.0
    assert len(radar.instrument_parameters['prt_mode']['data']) == 108
