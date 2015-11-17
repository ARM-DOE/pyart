""" Unit Tests for Py-ART's util/transforms.py module. """

from pyart.util import transforms


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
