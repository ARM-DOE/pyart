# check that the dummy NEXRAD file is simlar to non-dummy file.

from __future__ import print_function
import numpy as np
from numpy.testing import assert_array_equal, assert_allclose

import pyart

NEXRAD_FILE = 'KATX20130717_195021_V06'
OUTPUT_FILE = 'KATX20130717_195021_V06_DUMMY'


def test_dummy_similar():
    radar1 = pyart.io.read_nexrad_archive(NEXRAD_FILE)
    radar2 = pyart.io.read_nexrad_archive(OUTPUT_FILE)
    assert radars_similar(radar1, radar2)


def radars_similar(r1, r2):

    ###########################
    # Attribute that are None #
    ###########################
    assert dics_similar(r1.altitude_agl, r2.altitude_agl)
    assert dics_similar(r1.target_scan_rate, r2.target_scan_rate)
    assert dics_similar(r1.scan_rate, r2.scan_rate)
    assert dics_similar(r1.antenna_transition, r2.antenna_transition)
    assert dics_similar(r1.radar_calibration, r2.radar_calibration)

    #########################
    # Dictionary attributes #
    #########################
    #assert dics_similar(r1.time, r2.time)          # start time mismatch
    assert dics_similar(r1.range, r2.range)
    #assert dics_similar(r1.metadata, r2.metadata)  # DO not match

    assert dics_similar(r1.latitude, r2.latitude)
    assert dics_similar(r1.longitude, r2.longitude)
    #assert dics_similar(r1.altitude, r2.altitude)  # differ by 10 meters

    assert dics_similar(r1.sweep_number, r2.sweep_number)
    assert dics_similar(r1.sweep_mode, r2.sweep_mode)
    assert dics_similar(r1.fixed_angle, r2.fixed_angle)
    assert dics_similar(r1.sweep_start_ray_index, r2.sweep_start_ray_index)
    assert dics_similar(r1.sweep_end_ray_index, r2.sweep_end_ray_index)

    assert dics_similar(r1.azimuth, r2.azimuth)
    assert dics_similar(r1.elevation, r2.elevation)

    ###########
    # scalars #
    ###########

    assert r1.ngates == r2.ngates
    assert r1.nrays == r2.nrays
    assert r1.nsweeps == r2.nsweeps
    assert r1.scan_type == r2.scan_type

    ##########
    # fields #
    ##########

    print(r1.fields.keys())
    print(r2.fields.keys())
    assert set(r1.fields.keys()).difference(r2.fields.keys()) == set()

    for field in r1.fields:
        print(field)
        assert dics_similar(r1.fields[field], r2.fields[field])

    #radar1.fields
    return True


def dics_similar(dic1, dic2):
    """ Determine if two dictionaries are similar. """
    if dic1 is None:
        if dic2 is None:
            return True
        else:
            return False

    # all keys in both dictionaries
    print(dic1.keys())
    print(dic2.keys())
    assert set(dic1.keys()).difference(dic2.keys()) == set()

    for key in dic1.keys():
        if key == 'data':
            continue
        print(key, dic1[key], dic2[key])
        assert dic1[key] == dic2[key]

    # do not check data, since it should not match
    return True
