""" Unit Tests for Py-ART's io/arm_sonde.py module. """

import datetime

import numpy as np
from numpy.testing import assert_almost_equal
from pytest import raises

import pyart


def test_read_arm_sonde_vap_target_datetime():
    target_datetime = datetime.datetime(2011, 5, 10, 11, 30, 5)
    profile_datetime, hprofile = pyart.io.read_arm_sonde_vap(
        pyart.testing.INTERP_SOUNDE_FILE, target_datetime=target_datetime)

    assert profile_datetime == datetime.datetime(2011, 5, 10, 11, 30)
    assert_almost_equal(hprofile.height[:5], [318, 338, 358, 378, 398], 1)
    assert_almost_equal(hprofile.speed[:5], [8.4, 7.4, 9.2, 10.8, 12.3], 1)
    assert_almost_equal(
        hprofile.direction[:5], [182.6, 181.9, 182.3, 183.3, 184.2], 1)


def test_read_arm_sonde_vap_radar():
    radar = pyart.testing.make_empty_ppi_radar(1, 1, 1)
    radar.time['units'] = 'seconds since 2011-05-10T11:30:05Z'
    profile_datetime, hprofile = pyart.io.read_arm_sonde_vap(
        pyart.testing.INTERP_SOUNDE_FILE, radar=radar)

    assert profile_datetime == datetime.datetime(2011, 5, 10, 11, 30)
    assert_almost_equal(hprofile.height[:5], [318, 338, 358, 378, 398], 1)
    assert_almost_equal(hprofile.speed[:5], [8.4, 7.4, 9.2, 10.8, 12.3], 1)
    assert_almost_equal(
        hprofile.direction[:5], [182.6, 181.9, 182.3, 183.3, 184.2], 1)


def test_read_arm_sonde_vap_errors():

    # radar or target_datetime must be specified
    raises(
        ValueError, pyart.io.read_arm_sonde_vap,
        pyart.testing.INTERP_SOUNDE_FILE)

    # only one of radar or target_datetime can be specified
    radar = pyart.testing.make_empty_ppi_radar(1, 1, 1)
    target_datetime = datetime.datetime(2011, 5, 10, 11, 30, 5)
    raises(
        ValueError, pyart.io.read_arm_sonde_vap,
        pyart.testing.INTERP_SOUNDE_FILE,
        radar=radar, target_datetime=target_datetime)


def test_read_arm_sonde():
    profile_dt, hprofile = pyart.io.read_arm_sonde(pyart.testing.SONDE_FILE)

    assert profile_dt == datetime.datetime(2011, 5, 20, 8, 28)
    assert_almost_equal(hprofile.height[:5], [315, 321, 328, 336, 344], 0)
    assert_almost_equal(hprofile.speed[:5], [5, 3.2, 3.7, 4.3, 4.8], 1)
    assert_almost_equal(
        hprofile.direction[:5], [215., 193., 191., 191., 189.], 1)
