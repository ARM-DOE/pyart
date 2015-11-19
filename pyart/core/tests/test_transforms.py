""" Unit Tests for Py-ART's core/transforms.py module. """

import warnings

from pyart.core import transforms
from numpy.testing import assert_almost_equal


def test_corner_to_point():
    corner = (36.5, -97.5)
    point = (36.4, -97.6)
    x, y = transforms.corner_to_point(corner, point)
    assert round(x) == -8950.
    assert round(y) == -11119.0


def test_ax_radius_degrees():
    R = transforms._ax_radius(36.5, units='degrees')
    # answer found from: 6,371,000 * cos(36.5 / 180. * pi)
    assert round(R) == 5121372.


def test_ax_radius_radians():
    R = transforms._ax_radius(0.637045177, units='radians')
    assert round(R) == 5121372.


def test_cartesian_to_geographic_aeqd():
    # Example taken from:
    # Snyder, J.P. Map Projections A Working Manual, 1987, page 338.
    R = 3.0
    lon_0 = -100.        # 100 degrees West longitude
    lat_0 = 40.0        # 40 degrees North latitude
    x = -5.8311398
    y = 5.5444634

    lon, lat = transforms.cartesian_to_geographic_aeqd(x, y, lon_0, lat_0, R)
    assert_almost_equal(lat, -20.0, 3)  # 20.0 S latitude
    assert_almost_equal(lon, 100.0, 3)  # 100.0 E longitude

    # edge case, distance from projection center is zero
    x = y = 0
    with warnings.catch_warnings():  # invalid divide is handled by code
        warnings.simplefilter("ignore")
        lon, lat = transforms.cartesian_to_geographic_aeqd(
            x, y, lon_0, lat_0, R)
    assert_almost_equal(lon, -100.0, 3)
    assert_almost_equal(lat, 40.0, 3)
