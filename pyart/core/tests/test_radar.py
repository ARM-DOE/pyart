""" Unit Tests for Py-ART's core/radar.py module. """

import pickle
import sys
# we need a class which excepts str for writing in Python 2 and 3
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import inspect

import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal
import pytest

import pyart
from pyart.lazydict import LazyLoadDict


def test_radar_picklable():
    # verify that Radar instances are picklable
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    picklestring = pickle.dumps(radar)
    radar_new = pickle.loads(picklestring)
    assert 'data' in radar.gate_x
    assert 'data' in radar_new.gate_x


def test_gate_longitude_latitude():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    radar.azimuth['data'][:] = [0, 90, 180, 270, 0, 90, 180, 270]
    radar.elevation['data'][:] = [0, 0, 0, 0, 10, 10, 10, 10]
    radar.range['data'][:] = [5, 15, 25, 35, 45]

    assert radar.gate_longitude['data'].shape == (8, 5)
    assert radar.gate_latitude['data'].shape == (8, 5)
    assert_almost_equal(radar.gate_longitude['data'][0, 0], -97.5, 1)
    assert_almost_equal(radar.gate_latitude['data'][0, 0], 36.5, 1)

    # reset and try again with a non-default lat_0/lon_0
    radar.init_gate_longitude_latitude()
    radar.projection.pop('_include_lon_0_lat_0')
    radar.projection['lat_0'] = 20.0
    radar.projection['lon_0'] = 60.0

    assert_almost_equal(radar.gate_longitude['data'][0, 0], 60.0, 1)
    assert_almost_equal(radar.gate_latitude['data'][0, 0], 20.0, 1)


def test_gate_altitude():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    radar.azimuth['data'][:] = [0, 90, 180, 270, 0, 90, 180, 270]
    radar.elevation['data'][:] = [0, 0, 0, 0, 10, 10, 10, 10]
    radar.range['data'][:] = [5, 15, 25, 35, 45]

    assert radar.gate_altitude['data'].shape == (8, 5)
    assert_almost_equal(radar.gate_altitude['data'][0, 0], 200.0, 1)

    radar.init_gate_altitude()
    radar.altitude['data'][0] = 150.
    assert_almost_equal(radar.gate_altitude['data'][0, 0], 150.0, 1)


def test_gate_x_y_z():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    radar.azimuth['data'][:] = [0, 90, 180, 270, 0, 90, 180, 270]
    radar.elevation['data'][:] = [0, 0, 0, 0, 10, 10, 10, 10]
    radar.range['data'][:] = [5, 15, 25, 35, 45]

    assert radar.gate_x['data'].shape == (8, 5)
    assert_allclose(radar.gate_x['data'][0], [0, 0, 0, 0, 0], atol=1e-14)
    assert_allclose(radar.gate_x['data'][1], [5, 15, 25, 35, 45], atol=1e-14)
    assert_allclose(radar.gate_x['data'][2], [0, 0, 0, 0, 0], atol=1e-5)
    assert_allclose(
        radar.gate_x['data'][3], [-5, -15, -25, -35, -45], atol=1e-14)

    assert radar.gate_y['data'].shape == (8, 5)
    assert_allclose(radar.gate_y['data'][0], [5, 15, 25, 35, 45], atol=1e-14)
    assert_allclose(radar.gate_y['data'][1], [0, 0, 0, 0, 0], atol=1e-5)
    assert_allclose(
        radar.gate_y['data'][2], [-5, -15, -25, -35, -45], atol=1e-14)
    assert_allclose(radar.gate_y['data'][3], [0, 0, 0, 0, 0], atol=1e-6)

    assert radar.gate_z['data'].shape == (8, 5)
    z_sweep0 = np.array([1.47e-6, 1.324e-5, 3.679e-5, 7.210e-5, 1.1919e-4])
    assert_allclose(radar.gate_z['data'][0], z_sweep0, atol=1e-3)
    assert_allclose(radar.gate_z['data'][1], z_sweep0, atol=1e-3)
    assert_allclose(radar.gate_z['data'][2], z_sweep0, atol=1e-3)
    assert_allclose(radar.gate_z['data'][3], z_sweep0, atol=1e-3)


def test_get_gate_x_y_z():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    radar.azimuth['data'][:] = [0, 90, 180, 270, 0, 90, 180, 270]
    radar.elevation['data'][:] = [0, 0, 0, 0, 10, 10, 10, 10]
    radar.range['data'][:] = [5, 15, 25, 35, 45]

    gate_x, gate_y, gate_z = radar.get_gate_x_y_z(0)
    assert gate_x.shape == (4, 5)
    assert_allclose(gate_x[0], [0, 0, 0, 0, 0], atol=1e-5)
    assert_allclose(gate_x[1], [5, 15, 25, 35, 45], atol=1e-14)
    assert_allclose(gate_x[2], [0, 0, 0, 0, 0], atol=1e-5)
    assert_allclose(gate_x[3], [-5, -15, -25, -35, -45], atol=1e-14)

    assert gate_y.shape == (4, 5)
    assert_allclose(gate_y[0], [5, 15, 25, 35, 45], atol=1e-14)
    assert_allclose(gate_y[1], [0, 0, 0, 0, 0], atol=1e-5)
    assert_allclose(gate_y[2], [-5, -15, -25, -35, -45], atol=1e-14)
    assert_allclose(gate_y[3], [0, 0, 0, 0, 0], atol=1e-5)

    assert gate_z.shape == (4, 5)
    z_sweep0 = np.array([1.47e-6, 1.324e-5, 3.679e-5, 7.210e-5, 1.1919e-4])
    assert_allclose(gate_z[0], z_sweep0, atol=1e-3)
    assert_allclose(gate_z[1], z_sweep0, atol=1e-3)
    assert_allclose(gate_z[2], z_sweep0, atol=1e-3)
    assert_allclose(gate_z[3], z_sweep0, atol=1e-3)


def test_get_gate_x_y_z_edges():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 1)
    radar.azimuth['data'][:] = [315., 45, 135., 225.]
    radar.elevation['data'][:] = [0, 0, 0, 0]
    radar.range['data'][:] = [5, 15, 25, 35, 45]

    gate_x, gate_y, gate_z = radar.get_gate_x_y_z(0, edges=True)

    zeros = np.array([0, 0, 0, 0, 0, 0])
    even = np.array([0, 10, 20, 30, 40, 50])
    atol = 1e-4

    assert gate_x.shape == (5, 6)
    assert_allclose(gate_x[0], -even, atol=atol)
    assert_allclose(gate_x[1], zeros, atol=atol)
    assert_allclose(gate_x[2], even, atol=atol)
    assert_allclose(gate_x[3], zeros, atol=atol)
    assert_allclose(gate_x[4], -even, atol=atol)

    assert gate_y.shape == (5, 6)
    assert_allclose(gate_y[0], zeros, atol=atol)
    assert_allclose(gate_y[1], even, atol=atol)
    assert_allclose(gate_y[2], zeros, atol=atol)
    assert_allclose(gate_y[3], -even, atol=atol)
    assert_allclose(gate_y[4], zeros, atol=atol)

    assert gate_z.shape == (5, 6)
    z_sweep0 = np.array([0, 5.89e-6, 2.354e-5, 5.297e-5, 9.418e-5, 1.4715e-4])
    assert_allclose(gate_z[0], z_sweep0, atol=1e-3)
    assert_allclose(gate_z[1], z_sweep0, atol=1e-3)
    assert_allclose(gate_z[2], z_sweep0, atol=1e-3)
    assert_allclose(gate_z[3], z_sweep0, atol=1e-3)
    assert_allclose(gate_z[4], z_sweep0, atol=1e-3)


def test_get_gate_x_y_z_transitions():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    radar.azimuth['data'][:] = [0, 90, 180, 270, 0, 90, 180, 270]
    radar.elevation['data'][:] = [0, 0, 0, 0, 10, 10, 10, 10]
    radar.range['data'][:] = [5, 15, 25, 35, 45]
    radar.antenna_transition = {'data': np.array([0, 0, 1, 0, 0, 0, 0, 0])}

    gate_x, gate_y, gate_z = radar.get_gate_x_y_z(0, filter_transitions=True)
    assert gate_x.shape == (3, 5)
    assert_allclose(gate_x[0], [0, 0, 0, 0, 0], atol=1e-14)
    assert_allclose(gate_x[1], [5, 15, 25, 35, 45], atol=1e-14)
    assert_allclose(gate_x[2], [-5, -15, -25, -35, -45], atol=1e-14)

    assert gate_y.shape == (3, 5)
    assert_allclose(gate_y[0], [5, 15, 25, 35, 45], atol=1e-14)
    assert_allclose(gate_y[1], [0, 0, 0, 0, 0], atol=1e-5)
    assert_allclose(gate_y[2], [0, 0, 0, 0, 0], atol=1e-5)

    assert gate_z.shape == (3, 5)
    z_sweep0 = np.array([1.47e-6, 1.324e-5, 3.679e-5, 7.210e-5, 1.1919e-4])
    assert_allclose(gate_z[0], z_sweep0, atol=1e-3)
    assert_allclose(gate_z[1], z_sweep0, atol=1e-3)
    assert_allclose(gate_z[2], z_sweep0, atol=1e-3)


def test_get_gate_lat_lon_alt():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    lat, lon, alt = radar.get_gate_lat_lon_alt(0)
    assert lat.shape == (4, 5)
    assert_allclose(lat[0], [36.5, 36.502243, 36.50449, 36.506744, 36.50899],
                    atol=1e-3)
    assert_allclose(lat[1], [36.5, 36.502243, 36.50449, 36.506737, 36.508984],
                    atol=1e-3)
    assert_allclose(lat[2], [36.5, 36.502243, 36.50449, 36.506737, 36.508984],
                    atol=1e-3)
    assert_allclose(lat[3], [36.5, 36.50224, 36.504486, 36.50673 , 36.508976],
                    atol=1e-3)

    assert lon.shape == (4, 5)
    assert_allclose(lon[0], [-97.49999, -97.49999, -97.49999, -97.49999, -97.49999],
                    atol=1e-3)
    assert_allclose(lon[1], [-97.49999, -97.49995, -97.4999, -97.499855, -97.499794],
                    atol=1e-3)
    assert_allclose(lon[2], [-97.49999, -97.4999, -97.499794, -97.4997, -97.4996],
                    atol=1e-3)
    assert_allclose(lon[3], [-97.49999, -97.499855, -97.4997, -97.49956, -97.499405],
                    atol=1e-3)

    assert lat.shape == (4, 5)
    alt_sweep0 = np.array([200., 203., 206., 209., 213.])
    assert_allclose(alt[0], alt_sweep0, atol=1e-3)
    assert_allclose(alt[1], alt_sweep0, atol=1e-3)
    assert_allclose(alt[2], alt_sweep0, atol=1e-3)
    assert_allclose(alt[3], alt_sweep0, atol=1e-3)


def test_get_gate_lat_lon_alt_transitions():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 2)
    radar.antenna_transition = {'data': np.array([0, 0, 1, 0, 0, 0, 0, 0])}
    lat, lon, alt = radar.get_gate_lat_lon_alt(0, filter_transitions=True)
    assert lat.shape == (3, 5)
    assert_allclose(lat[0], [36.5, 36.502243, 36.50449, 36.506744, 36.50899],
                    atol=1e-3)
    assert_allclose(lat[1], [36.5, 36.502243, 36.50449, 36.506737, 36.508984],
                    atol=1e-3)
    assert_allclose(lat[2], [36.5, 36.50224, 36.504486, 36.50673, 36.508976],
                    atol=1e-3)

    assert lon.shape == (3, 5)
    assert_allclose(lon[0], [-97.49999, -97.49999, -97.49999, -97.49999, -97.49999],
                    atol=1e-3)
    assert_allclose(lon[1], [-97.49999, -97.49995, -97.4999, -97.499855, -97.499794],
                    atol=1e-3)
    assert_allclose(lon[2], [-97.49999, -97.499855, -97.4997, -97.49956, -97.499405],
                    atol=1e-3)

    assert lat.shape == (3, 5)
    alt_sweep0 = np.array([200., 203., 206., 209., 213.])
    assert_allclose(alt[0], alt_sweep0, atol=1e-3)
    assert_allclose(alt[1], alt_sweep0, atol=1e-3)
    assert_allclose(alt[2], alt_sweep0, atol=1e-3)


def test_init_gate_x_y_z():
    radar = pyart.testing.make_empty_ppi_radar(5, 4, 1)
    radar.azimuth['data'][:] = [0, 90, 180, 270]
    radar.elevation['data'][:] = [0, 0, 0, 0]
    radar.range['data'][:] = [5, 15, 25, 35, 45]

    # access and check initial gate locations
    assert_allclose(radar.gate_x['data'][1], [5, 15, 25, 35, 45], atol=1e-14)
    assert_allclose(radar.gate_y['data'][0], [5, 15, 25, 35, 45], atol=1e-14)
    z_sweep0 = np.array([1.47e-6, 1.324e-5, 3.679e-5, 7.210e-5, 1.1919e-4])
    assert_allclose(radar.gate_z['data'][0], z_sweep0, atol=1e-3)

    # change range, gate_x, y, z are not updated
    radar.range['data'][:] = [15, 25, 35, 45, 55]
    assert_allclose(radar.gate_x['data'][1], [5, 15, 25, 35, 45], atol=1e-14)
    assert_allclose(radar.gate_y['data'][0], [5, 15, 25, 35, 45], atol=1e-14)
    z_sweep0 = np.array([1.47e-6, 1.324e-5, 3.679e-5, 7.210e-5, 1.1919e-4])
    assert_allclose(radar.gate_z['data'][0], z_sweep0, atol=1e-3)

    # call init_gate_x_y_z, now the attributes are updated
    radar.init_gate_x_y_z()
    assert_allclose(radar.gate_x['data'][1], [15, 25, 35, 45, 55], atol=1e-14)
    assert_allclose(radar.gate_y['data'][0], [15, 25, 35, 45, 55], atol=1e-14)
    z_sweep0 = np.array([1.324e-5, 3.679e-5, 7.210e-5, 1.1919e-4, 1.7805e-4])
    assert_allclose(radar.gate_z['data'][0], z_sweep0, atol=1e-3)


def test_rays_per_sweep_attribute():
    radar = pyart.testing.make_target_radar()
    rays_per_sweep = radar.rays_per_sweep
    assert isinstance(rays_per_sweep, LazyLoadDict)
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
    pytest.raises(KeyError, radar.iter_field, 'foobar')

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
    pytest.raises(IndexError, radar.get_start, -1)
    pytest.raises(IndexError, radar.get_start, 20)

    assert radar.get_end(0) == 19
    assert radar.get_end(1) == 39
    pytest.raises(IndexError, radar.get_end, -1)
    pytest.raises(IndexError, radar.get_end, 20)

    assert radar.get_start_end(0) == (0, 19)
    assert radar.get_start_end(1) == (20, 39)
    pytest.raises(IndexError, radar.get_start_end, -1)
    pytest.raises(IndexError, radar.get_start_end, 20)

    assert radar.get_slice(0) == slice(0, 20)
    assert radar.get_slice(1) == slice(20, 40)
    pytest.raises(IndexError, radar.get_slice, -1)
    pytest.raises(IndexError, radar.get_slice, 20)

    data = radar.get_field(0, 'reflectivity')
    assert data.shape == (20, 30)
    assert data.dtype == np.float32
    data = radar.get_field(1, 'reflectivity')
    assert data.shape == (20, 30)
    assert data.dtype == np.float32
    pytest.raises(KeyError, radar.get_field, 0, 'foobar')
    pytest.raises(IndexError, radar.get_field, -1, 'reflectivity')
    pytest.raises(IndexError, radar.get_field, 20, 'reflectivity')

    assert radar.get_azimuth(0).shape == (20, )
    pytest.raises(IndexError, radar.get_azimuth, -1)
    pytest.raises(IndexError, radar.get_azimuth, 20)

    assert radar.get_elevation(0).shape == (20, )
    pytest.raises(IndexError, radar.get_elevation, -1)
    pytest.raises(IndexError, radar.get_elevation, 20)

    pytest.raises(LookupError, radar.get_nyquist_vel, 0)
    radar.instrument_parameters = {
        'nyquist_velocity': {'data': np.ones((100,))}
    }
    assert round(radar.get_nyquist_vel(0)) == 1
    pytest.raises(IndexError, radar.get_nyquist_vel, -1)
    radar.instrument_parameters['nyquist_velocity']['data'][0] = 2
    pytest.raises(Exception, radar.get_nyquist_vel, 0)


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
    pytest.raises(ValueError, radar.extract_sweeps, [0, 2])
    pytest.raises(ValueError, radar.extract_sweeps, [-1, 1])


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

    pytest.raises(ValueError, radar.add_field, 'reflectivity', {})

    dic = {'dat': np.zeros((360, 50)), 'standard_name': 'test'}
    pytest.raises(KeyError, radar.add_field, 'test', dic)

    dic = {'data': np.zeros((360, 49)), 'standard_name': 'test'}
    pytest.raises(ValueError, radar.add_field, 'test', dic)


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
    pytest.raises(ValueError, radar.add_field_like, 'foo', 'bar', [])


@pytest.mark.parametrize(
    "level", ['standard', 's', 'compact', 'c', 'full', 'f'])
def test_info_levels(level):
    check_info(level)


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
    pytest.raises(ValueError, check_info, 'foo')
