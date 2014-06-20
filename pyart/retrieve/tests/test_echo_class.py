""" Unit Tests for Py-ART's io/mdv.py module. """

import numpy as np
from numpy.testing.decorators import skipif

import pyart


@skipif(not pyart.retrieve._F90_EXTENSIONS_AVAILABLE)
def test_steiner_conv_strat():
    grid = pyart.testing.make_storm_grid()
    eclass = pyart.retrieve.steiner_conv_strat(grid)
    assert np.all(eclass['data'][25] == np.array(
        [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
         2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]))
