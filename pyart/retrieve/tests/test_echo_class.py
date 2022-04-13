""" Unit Tests for Py-ART's retrieve/echo_class.py module. """

import numpy as np
import pytest

import pyart


def test_steiner_conv_strat_default():
    grid = pyart.testing.make_storm_grid()
    eclass_default = pyart.retrieve.steiner_conv_strat(grid)
    assert np.all(eclass_default['data'][25] == np.array(
        [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
         2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]))

@pytest.mark.parametrize('area_relation', ['small', 'medium', 'large', 'sgp'])
def test_steiner_conv_strat_modify_area(area_relation):
    grid = pyart.testing.make_storm_grid()
    eclass = pyart.retrieve.steiner_conv_strat(grid,
                                               area_relation=area_relation)
    assert eclass['data'].min() == 0
    assert eclass['data'].max() == 2
