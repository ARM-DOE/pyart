""" Unit tests for rainfall rate estimation module. """

import numpy as np
from numpy.testing import assert_allclose

import pyart


def test_est_rain_rate_zpoly():
    radar = pyart.io.read(pyart.testing.UF_FILE)
    rain = pyart.retrieve.est_rain_rate_zpoly(radar)

    assert 'units' in rain.keys()
    assert 'long_name' in rain.keys()
    assert 'coordinates' in rain.keys()
    assert rain['units'] == 'mm/hr'

    assert_allclose(
        rain['data'][0][0:5],
        [2.89949111e-04, 1.26031121e-02, 9.54675995e-06,
         2.14264704e-01, 9.47527914e-01], atol=1e7)

    assert_allclose(
        rain['data'][0][-5:],
        [0.20666919, 0.10129632, 0.13916402, 0.12526962,
         0.12681702], atol=1e7)
    del radar


def test_est_rain_rate_z():
    radar = pyart.io.read(pyart.testing.UF_FILE)
    rain = pyart.retrieve.est_rain_rate_z(radar)

    assert 'units' in rain.keys()
    assert 'long_name' in rain.keys()
    assert 'coordinates' in rain.keys()
    assert rain['units'] == 'mm/hr'

    assert_allclose(
        rain['data'][0][0:5],
        [0.01604766, 0.05375671, 0.00767614,
         0.27197879, 1.04875935], atol=1e7)

    assert_allclose(
        rain['data'][0][-5:],
        [0.26443019, 0.15909914, 0.1973247,
         0.18339988, 0.18495507], atol=1e7)
    del radar


def test_est_rain_rate_kdp():
    radar = pyart.io.read(pyart.testing.UF_FILE)
    rain = pyart.retrieve.est_rain_rate_kdp(radar)

    assert 'units' in rain.keys()
    assert 'long_name' in rain.keys()
    assert 'coordinates' in rain.keys()
    assert rain['units'] == 'mm/hr'

    assert_allclose(
        rain['data'][0][0:5],
        [0., 0., 0., 0., 0.], atol=1e7)

    assert_allclose(
        rain['data'][0][-5:],
        [57.4029591, 56.8539739, 56.57898183,
         56.30365354, 56.11991333], atol=1e7)

    radar.instrument_parameters.pop('frequency')
    rain_no_freq = pyart.retrieve.est_rain_rate_kdp(radar)

    assert_allclose(
        rain_no_freq['data'][0][0:5],
        [0., 0., 0., 0., 0.], atol=1e7)

    assert_allclose(
        rain_no_freq['data'][0][-5:],
        [117.04546973, 115.85529152, 115.25939369,
         114.66295169, 114.26501948], atol=1e7)
    del radar


def test_est_rain_rate_a():
    radar = pyart.io.read(pyart.testing.UF_FILE)
    spec_a, _ = pyart.correct.calculate_attenuation(radar, z_offset=1.0)
    radar.add_field('specific_attenuation', spec_a, replace_existing=True)
    rain = pyart.retrieve.est_rain_rate_a(radar)

    assert 'units' in rain.keys()
    assert 'long_name' in rain.keys()
    assert 'coordinates' in rain.keys()
    assert rain['units'] == 'mm/hr'

    assert_allclose(
        rain['data'][0][0:5],
        [0.03885328, 0.02587043, 0.02171088,
         0.13686746, 1.2150548], atol=1e7)

    assert_allclose(
        rain['data'][0][-5:],
        [0., 0., 0., 0., 0.], atol=1e7)

    radar.instrument_parameters.pop('frequency')
    rain_no_freq = pyart.retrieve.est_rain_rate_kdp(radar)

    assert_allclose(
        rain_no_freq['data'][0][0:5],
        [0.10804109, 0.06917369, 0.05707911,
         0.42970675, 4.708362], atol=1e7)

    assert_allclose(
        rain_no_freq['data'][0][-5:],
        [0., 0., 0., 0., 0.], atol=1e7)
    del radar


def test_coeff_ra_table():
    coeff_ra_dict = pyart.retrieve.qpe._coeff_ra_table()
    test_dict = {'S': (3100.0, 1.03),
                 'C': (250.0, 0.91),
                 'X': (45.5, 0.83)}

    assert coeff_ra_dict == test_dict


def test_coeff_rkdp_table():
    coeff_rkdp_dict = pyart.retrieve.qpe._coeff_rkdp_table()
    test_dict = {'S': (50.7, 0.85),
                 'C': (29.7, 0.85),
                 'X': (15.81, 0.7992)}

    assert coeff_rkdp_dict == test_dict


def test_get_coeff_ra():
    coeff_ra_good = pyart.retrieve.qpe._get_coeff_ra(3e9)
    assert coeff_ra_good == (3100.0, 1.03)

    coeff_ra_use_s = pyart.retrieve.qpe._get_coeff_ra(1e9)
    assert coeff_ra_use_s == (3100.0, 1.03)

    coeff_ra_use_x = pyart.retrieve.qpe._get_coeff_ra(13e9)
    assert coeff_ra_use_x == (45.5, 0.83)


def test_get_coeff_rkdp():
    coeff_rkdp_good = pyart.retrieve.qpe._get_coeff_rkdp(3e9)
    assert coeff_rkdp_good == (50.7, 0.85)

    coeff_rkdp_use_s = pyart.retrieve.qpe._get_coeff_rkdp(1e9)
    assert coeff_rkdp_use_s == (50.7, 0.85)

    coeff_rkdp_use_x = pyart.retrieve.qpe._get_coeff_rkdp(13e9)
    assert coeff_rkdp_use_x == (15.81, 0.7992)
