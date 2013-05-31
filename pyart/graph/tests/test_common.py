""" Unit Tests for Py-ART's graph/common.py module. """

import pyart


def test_corner_to_point():
    corner = (36.5, -97.5)
    point = (36.4, -97.6)
    x, y = pyart.graph.common.corner_to_point(corner, point)
    assert round(x) == -8950.
    assert round(y) == -11119.0


def test_ax_radius_degrees():
    R = pyart.graph.common.ax_radius(36.5, units='degrees')
    # answer found from: 6,371,000 * cos(36.5 / 180. * pi)
    assert round(R) == 5121372.


def test_ax_radius_radians():
    R = pyart.graph.common.ax_radius(0.637045177, units='radians')
    assert round(R) == 5121372.
