""" Unit Tests for Py-ART's retrieve/advection.py module. """

import numpy as np
from numpy.testing.decorators import skipif

import pyart

def test_grid_displacememt_pc():
    test_storm = pyart.testing.make_storm_grid()
    test_storm_2.fields['reflectivity']['data'] =\
            test_storm_2.fields['reflectivity']['data'][:, 9:-11, 9:-1]
    test_storm.fields['reflectivity']['data'] =\
            test_storm.fields['reflectivity']['data'][:, 5:-5, :]
    test_storm_2 = pyart.testing.make_storm_grid()
    test_storm_2.fields['reflectivity']['data'] =\
            test_storm_2.fields['reflectivity']['data'][:, 0:-10, :]
    test_storm_2.fields['reflectivity']['data'][:,:,4:-3] =\
            test_storm_2.fields['reflectivity']['data'][:, :, 1:-6]
    test_storm_2.fields['reflectivity']['valid_min']=0.0
    test_storm.fields['reflectivity']['valid_min']=0.0
    assert pyart.retrieve.grid_displacememt(test_storm,
            test_storm_2, 'reflectivity', 0) == (3,5)


