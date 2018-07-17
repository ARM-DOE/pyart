""" Unit Tests for Py-ART's io/uf.py and io/uffile.py modules. """

from datetime import datetime
try:
    from StringIO import StringIO
except ImportError:
    from io import BytesIO as StringIO

import numpy as np
from numpy.testing import assert_raises, assert_almost_equal, assert_warns
import pytest

import pyart
from pyart.io.uffile import UFFile, UFRay

radar = pyart.io.read_uf(pyart.testing.UF_FILE, file_field_names=True)


def test_time():
    assert radar.time['units'] == 'seconds since 2011-05-20T10:54:16Z'
    assert_almost_equal(radar.time['data'][0], 0)


def test_range():
    assert np.allclose(radar.range['data'], np.arange(667) * 60 + 30)
    assert radar.range['meters_to_center_of_first_gate'] == 30
    assert radar.range['meters_between_gates'] == 60


def test_lat_lon_alt():
    assert_almost_equal(radar.latitude['data'], 36.49, 2)
    assert_almost_equal(radar.longitude['data'], -97.59, 2)
    assert_almost_equal(radar.altitude['data'], 214)


def test_sweep_start_ray_index():
    assert np.allclose(radar.sweep_start_ray_index['data'], np.array([0]))


def test_sweep_end_ray_index():
    assert np.allclose(radar.sweep_end_ray_index['data'], np.array([0]))


def test_sweep_number():
    assert np.allclose(radar.sweep_number['data'], np.array([0]))


def test_scan_type():
    assert radar.scan_type == 'ppi'


def test_sweep_mode():
    assert radar.sweep_mode['data'].dtype.char == 'S'
    assert radar.sweep_mode['data'] == np.array([b'azimuth_surveillance'])


def test_fixed_angle():
    assert np.allclose(radar.fixed_angle['data'], np.array([0.5]))


def test_elevation():
    assert np.allclose(radar.elevation['data'], np.array([0.484375]))


def test_azimuth():
    assert np.allclose(radar.azimuth['data'], np.array([359.9375]))


fields = {'CZ': 0.0, 'DR': 0.0, 'DZ': -6.0499999999999998,
          'HC': 1.0, 'KD': 0.0, 'PH': 90.0, 'RH': 1.0,
          'SQ': 1.0, 'SW': 0.01, 'VR': 0.0, 'ZD': 0.0,
          'ZT': -6.0499999999999998}
@pytest.mark.parametrize(
    "field, field_value", fields.items(), ids=list(fields.keys()))
def test_fields(field, field_value):
    check_field(field, field_value)


def check_field(field, value):
    assert_almost_equal(radar.fields[field]['data'][0, 0], value)


def test_raises_ioerror():
    fake_bad_file = StringIO(b'XXXXXXXX')
    pytest.raises(IOError, pyart.io.read_uf, fake_bad_file)


def test_read_fileobj():
    fh = open(pyart.testing.UF_FILE, 'rb')
    radar = pyart.io.read_uf(fh)
    fh.close()


def test_instrument_parameters():
    # test that instrument parameters exist
    assert radar.scan_rate is not None
    assert 'pulse_width' in radar.instrument_parameters
    assert 'radar_beam_width_h' in radar.instrument_parameters
    assert 'radar_beam_width_v' in radar.instrument_parameters
    assert 'radar_receiver_bandwidth' in radar.instrument_parameters
    assert 'polarization_mode' in radar.instrument_parameters
    assert 'frequency' in radar.instrument_parameters
    assert 'prt' in radar.instrument_parameters
    assert 'nyquist_velocity' in radar.instrument_parameters


def test_nyquist_vel():
    ufile = UFFile(pyart.testing.UF_FILE)
    ufile.rays[0].field_headers[1].pop('nyquist')
    assert ufile.get_nyquists() is None

    ufile = UFFile(pyart.testing.UF_FILE)
    ray = UFRay(ufile.rays[0]._buf)
    ray.field_headers[1].pop('nyquist')
    ufile.rays.append(ray)
    ufile.nrays = 2
    assert ufile.get_nyquists() is None


def test_datetime():
    ufile = UFFile(pyart.testing.UF_FILE)
    ray = ufile.rays[0]
    assert ray.get_datetime() == datetime(2011, 5, 20, 10, 54, 16)

    # test case where midnight is incorrectly represented as 24:00:00
    ray.mandatory_header['hour'] = 24
    ray.mandatory_header['minute'] = 0
    ray.mandatory_header['second'] = 0
    assert ray.get_datetime() == datetime(2011, 5, 21, 0, 0, 0)


def test_polarization():
    ufile = UFFile(pyart.testing.UF_FILE)
    ufile.rays[0].field_headers[0]['polarization'] = 99
    assert ufile.get_sweep_polarizations()[0] == 'elliptical'


def test_frequency():
    ufile = UFFile(pyart.testing.UF_FILE)
    filemetadata = pyart.config.FileMetadata('uf')

    ip = pyart.io.uf._get_instrument_parameters(ufile, filemetadata)
    assert 'frequency' in ip
    assert_almost_equal(ip['frequency']['data'][0], 9.69026150e+09, -3)

    # An invalid wavelength should throw a warning
    ufile.rays[0].field_headers[0]['wavelength_cm'] = 0
    assert_warns(UserWarning,
                 pyart.io.uf._get_instrument_parameters, ufile, filemetadata)


def test_scan_type():
    ufile = UFFile(pyart.testing.UF_FILE)
    ufray = ufile.rays[0]

    scan_type = pyart.io.uf._get_scan_type(ufray)
    assert scan_type == 'ppi'

    # An invalid scan mode should throw a warning
    ufray.mandatory_header['sweep_mode'] = 99
    assert_warns(UserWarning, pyart.io.uf._get_scan_type, ufray)


def test_skip_field():
    test_radar = pyart.io.read_uf(
        pyart.testing.UF_FILE, exclude_fields=['DZ'], file_field_names=True)
    assert 'DZ' not in test_radar.fields.keys()


def test_range_start_nonzero():
    radar = pyart.io.read_uf(pyart.testing.UF_FILE)
    radar.range['meters_to_center_of_first_gate'] += 1500.
    in_mem = StringIO()
    pyart.io.write_uf(in_mem, radar)

    in_mem.seek(0)
    radar2 = pyart.io.read_uf(in_mem)
    assert_almost_equal(radar2.range['meters_to_center_of_first_gate'], 1530)
    assert_almost_equal(radar2.range['data'][0], 1530)
