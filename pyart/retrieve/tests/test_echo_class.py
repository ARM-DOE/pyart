""" Unit Tests for Py-ART's retrieve/echo_class.py module. """

import numpy as np
import pytest

import pyart


def test_steiner_conv_strat():
    grid = pyart.testing.make_storm_grid()
    eclass = pyart.retrieve.steiner_conv_strat(grid)
    assert np.all(eclass['data'][25] == np.array(
        [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
         2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]))
