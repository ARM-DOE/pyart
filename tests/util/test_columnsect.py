""" Unit Tests for Py-ART's util/columnsect.py module. """

import numpy as np

import pyart

# read in example file
radar = pyart.io.read_nexrad_archive(pyart.testing.NEXRAD_ARCHIVE_MSG31_FILE)


def test_get_azimuth():
    # test to make sure azimuth is correct to Everett, WA
    azimuth = pyart.util.columnsect.for_azimuth(
        radar.latitude["data"][0], 47.97, radar.longitude["data"][0], -122.20
    )
    test_azi = abs(np.around(azimuth, 2) - 138.57)
    assert test_azi < 0.001


def test_sphere_distance():
    # test to make sure sphere distance is correct to Everett, WA
    distance = pyart.util.columnsect.sphere_distance(
        radar.latitude["data"][0], 47.97, (radar.longitude["data"][0]), -122.20
    )
    test_dis = abs(np.around(distance / 1000.0, 2) - 33.27)
    assert test_dis < 0.001


def test_get_field_location():
    # test to make sure column above location is pulled correctly
    column = pyart.util.columnsect.get_field_location(radar, 47.97, -122.20)

    # check to make sure z-gate is pulled correctly.
    test_height = abs(column.height[0] - 564)
    assert test_height < 0.001

    # check to make sure reflectivity value is minimum
    test_z = abs(column.reflectivity[0] + 32)
    assert test_z < 0.001
