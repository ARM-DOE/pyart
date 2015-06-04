""" Unit Tests for Py-ART's core/radar.py module. """

import sys
# we need a class which excepts str for writing in Python 2 and 3
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import inspect

import numpy as np
from numpy.testing import assert_raises
import pyart


def test_rays_per_sweep_attribute():
    radar = pyart.testing.make_target_radar()
    rays_per_sweep = radar.rays_per_sweep
    assert isinstance(rays_per_sweep, dict)
    assert rays_per_sweep['data'].shape == (1, )
    assert rays_per_sweep['data'][0] == 360


def test_iterators():
    radar = pyart.testing.make_empty_ppi_radar(30, 20, 5)
    radar.fields['reflectivity'] = {
        'data': np.zeros((100, 30), dtype=np.float32)}

    starts = [0, 20, 40, 60, 80]
    ends = [19, 39, 59, 79, 99]
    starts_ends = [(s, e) for s, e in zip(starts, ends)]

    assert inspect.isgenerator(radar.iter_start())
    assert [s for s in radar.iter_start()] == starts

    assert inspect.isgenerator(radar.iter_end())
    assert [s for s in radar.iter_end()] == ends

    assert inspect.isgenerator(radar.iter_start_end())
    assert [s for s in radar.iter_start_end()] == starts_ends

    assert inspect.isgenerator(radar.iter_slice())
    for s, start, end in zip(radar.iter_slice(), starts, ends):
        assert s.start == start
        assert s.stop == end + 1
        assert s.step is None

    assert inspect.isgenerator(radar.iter_field('reflectivity'))
    for d in radar.iter_field('reflectivity'):
        assert d.shape == (20, 30)
        assert d.dtype == np.float32
    assert_raises(KeyError, radar.iter_field, 'foobar')

    assert inspect.isgenerator(radar.iter_azimuth())
    for d in radar.iter_azimuth():
        assert d.shape == (20, )

    assert inspect.isgenerator(radar.iter_elevation())
    for d in radar.iter_elevation():
        assert d.shape == (20, )


def test_get_methods():
    radar = pyart.testing.make_empty_ppi_radar(30, 20, 5)
    radar.fields['reflectivity'] = {
        'data': np.zeros((100, 30), dtype=np.float32)}

    assert radar.get_start(0) == 0
    assert radar.get_start(1) == 20
    assert_raises(IndexError, radar.get_start, -1)
    assert_raises(IndexError, radar.get_start, 20)

    assert radar.get_end(0) == 19
    assert radar.get_end(1) == 39
    assert_raises(IndexError, radar.get_end, -1)
    assert_raises(IndexError, radar.get_end, 20)

    assert radar.get_start_end(0) == (0, 19)
    assert radar.get_start_end(1) == (20, 39)
    assert_raises(IndexError, radar.get_start_end, -1)
    assert_raises(IndexError, radar.get_start_end, 20)

    assert radar.get_slice(0) == slice(0, 20)
    assert radar.get_slice(1) == slice(20, 40)
    assert_raises(IndexError, radar.get_slice, -1)
    assert_raises(IndexError, radar.get_slice, 20)

    data = radar.get_field(0, 'reflectivity')
    assert data.shape == (20, 30)
    assert data.dtype == np.float32
    data = radar.get_field(1, 'reflectivity')
    assert data.shape == (20, 30)
    assert data.dtype == np.float32
    assert_raises(KeyError, radar.get_field, 0, 'foobar')
    assert_raises(IndexError, radar.get_field, -1, 'reflectivity')
    assert_raises(IndexError, radar.get_field, 20, 'reflectivity')

    assert radar.get_azimuth(0).shape == (20, )
    assert_raises(IndexError, radar.get_azimuth, -1)
    assert_raises(IndexError, radar.get_azimuth, 20)

    assert radar.get_elevation(0).shape == (20, )
    assert_raises(IndexError, radar.get_elevation, -1)
    assert_raises(IndexError, radar.get_elevation, 20)

    assert_raises(LookupError, radar.get_nyquist_vel, 0)
    radar.instrument_parameters = {
        'nyquist_velocity': {'data': np.ones((100,))}
    }
    assert round(radar.get_nyquist_vel(0)) == 1
    assert_raises(IndexError, radar.get_nyquist_vel, -1)
    radar.instrument_parameters['nyquist_velocity']['data'][0] = 2
    assert_raises(Exception, radar.get_nyquist_vel, 0)


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
    assert isinstance(radar, pyart.core.Radar)


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


def test_add_field_like_bug():
    # tests for bug where adding a field over-writes 'like' field
    # data/metadata.
    radar = pyart.testing.make_target_radar()
    data = np.ones((360, 50))
    radar.add_field_like('reflectivity', 'test', data)
    radar.fields['test']['units'] = 'fake'

    # check field added
    assert radar.fields['test']['units'] == 'fake'
    assert radar.fields['test']['data'][0, 0] == 1

    # check original field
    assert radar.fields['reflectivity']['units'] == 'dBZ'
    assert radar.fields['reflectivity']['data'][0, 0] == 0


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
    out = StringIO()
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
    assert not pyart.core.is_vpt(radar)
    pyart.core.to_vpt(radar)
    assert pyart.core.is_vpt(radar)


def test_to_vpt():
    # single scan
    radar = pyart.testing.make_empty_ppi_radar(10, 36, 3)
    radar.instrument_parameters = {
        'prt_mode': {'data': np.array(['fixed'] * 3)}
    }
    pyart.core.to_vpt(radar)
    assert pyart.core.is_vpt(radar)
    assert radar.nsweeps == 1
    assert radar.azimuth['data'][10] == 0.0
    assert radar.elevation['data'][0] == 90.0
    assert len(radar.instrument_parameters['prt_mode']['data']) == 1

    # multiple scans
    radar = pyart.testing.make_empty_ppi_radar(10, 36, 3)
    radar.instrument_parameters = {
        'prt_mode': {'data': np.array(['fixed'] * 3)}
    }
    pyart.core.to_vpt(radar, False)
    assert pyart.core.is_vpt(radar)
    assert radar.nsweeps == 108
    assert radar.azimuth['data'][10] == 10.0
    assert radar.elevation['data'][0] == 90.0
    assert len(radar.instrument_parameters['prt_mode']['data']) == 108
