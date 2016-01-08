""" Unit Tests for Py-ART's core/transforms.py module. """

import warnings

from pyart.core import transforms
from numpy.testing import assert_almost_equal
from numpy.testing.decorators import skipif


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


def test_geographic_to_cartesian_aeqd():
    # Example taken from:
    # Snyder, J.P. Map Projections A Working Manual, 1987, page 338.
    R = 3.0
    lat_0 = 40.0        # 40 degrees North latitude
    lon_0 = -100.       # 100 degrees West longitude
    lat = -20.0         # 20 degrees S latitude
    lon = 100.0         # 100.0 E longitude
    x = -5.8311398
    y = 5.5444634

    with warnings.catch_warnings():  # invalid divide is handled by code
        warnings.simplefilter("ignore", category=RuntimeWarning)
        x, y = transforms.geographic_to_cartesian_aeqd(
            lon, lat, lon_0, lat_0, R)
    assert_almost_equal(x, -5.8311398, 7)
    assert_almost_equal(y, 5.5444634, 7)

    # edge case, distance from projection center is zero
    lat = 40.0
    lon = -100.
    with warnings.catch_warnings():  # invalid divide is handled by code
        # ignore division runtime warning
        warnings.simplefilter("ignore", category=RuntimeWarning)
        x, y = transforms.geographic_to_cartesian_aeqd(
            lon, lat, lon_0, lat_0, R)
    assert_almost_equal(x, 0.0, 5)
    assert_almost_equal(y, 0.0, 5)

    # edge case, sin(c) is zero
    with warnings.catch_warnings():  # invalid divide is handled by code
        warnings.simplefilter("ignore", category=RuntimeWarning)
        x, y = transforms.geographic_to_cartesian_aeqd(
            10.0, 90.0, 20.0, 90.0, 3.0)
    assert_almost_equal(x, 0.0, 5)
    assert_almost_equal(y, 0.0, 5)


def test_cartesian_to_geographic():
    # Example taken from:
    # Snyder, J.P. Map Projections A Working Manual, 1987, page 338.
    R = 3.0
    lon_0 = -100.        # 100 degrees West longitude
    lat_0 = 40.0        # 40 degrees North latitude
    x = -5.8311398
    y = 5.5444634

    projparams = {
        'proj': 'pyart_aeqd',
        'lat_0': lat_0,
        'lon_0': lon_0,
        'R': R,
    }
    lon, lat = transforms.cartesian_to_geographic(x, y, projparams)
    assert_almost_equal(lat, -20.0, 3)  # 20.0 S latitude
    assert_almost_equal(lon, 100.0, 3)  # 100.0 E longitude

    # Use the default R value
    x *= 6370997. / 3.
    y *= 6370997. / 3.
    projparams.pop('R')
    lon, lat = transforms.cartesian_to_geographic(x, y, projparams)
    assert_almost_equal(lat, -20.0, 3)  # 20.0 S latitude
    assert_almost_equal(lon, 100.0, 3)  # 100.0 E longitude


@skipif(not transforms._PYPROJ_AVAILABLE)
def test_cartesian_to_geographic_pyproj():
    R = 3.0
    lon_0 = -100.        # 100 degrees West longitude
    lat_0 = 40.0        # 40 degrees North latitude
    x = -5.8311398
    y = 5.5444634
    projparams = {
        'proj': 'aeqd',
        'lat_0': lat_0,
        'lon_0': lon_0,
        'R': R,
    }
    projparams['proj'] = 'aeqd'
    lon, lat = transforms.cartesian_to_geographic(x, y, projparams)
    assert_almost_equal(lat, -20.0, 3)  # 20.0 S latitude
    assert_almost_equal(lon, 100.0, 3)  # 100.0 E longitude


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
