""" Unit Tests for Py-ART's retrieve/precip_rate.py module. """

import pyart
import math


def test_precip_rate():
    grid = pyart.testing.make_storm_grid()
    dict = pyart.retrieve.ZtoR(grid)

    # check that field is in grid object
    assert "NWS_primary_prate" in grid.fields.keys()

    # check calculations are within 10^-4 orders of magnitude
    assert math.floor(
        math.log10(
            grid.fields['NWS_primary_prate']['data'][0,10,10] -
            (((10**(grid.fields['reflectivity']['data'][0,10,10]/10)) / 300) ** (1 / 1.4)))) < -4




