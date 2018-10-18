""" Unit Tests for Py-ART's io/mdv.py module. """

import numpy as np
import pytest

import pyart


@pytest.mark.skipif(not pyart.retrieve.echo_class._F90_EXTENSIONS_AVAILABLE,
                    reason=("Py-ART was not built on a system with a Fortran "
                            "Compiler."))
def test_steiner_conv_strat():
    grid = pyart.testing.make_storm_grid()
    eclass = pyart.retrieve.steiner_conv_strat(grid)
    assert np.all(eclass['data'][25] == np.array(
        [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
         2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]))
