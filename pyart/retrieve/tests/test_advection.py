""" Unit Tests for Py-ART's retrieve/advection.py module. """

import numpy as np
from numpy.testing.decorators import skipif

import pyart

@skipif(not pyart.retrieve._ADVECTION_AVAILABLE)
def test_grid_displacement_pc():
    test_storm = pyart.testing.make_storm_grid()
    test_storm.fields['reflectivity']['data'] =\
            test_storm.fields['reflectivity']['data'][:, 5:-5, :]
    test_storm_2 = pyart.testing.make_storm_grid()
    test_storm_2.fields['reflectivity']['data'] =\
            test_storm_2.fields['reflectivity']['data'][:, 0:-10, :]
    test_storm_2.fields['reflectivity']['data'][:,:,4:-3] =\
            test_storm_2.fields['reflectivity']['data'][:, :, 1:-6]
    test_storm_2.fields['reflectivity']['valid_min']=0.0
    test_storm.fields['reflectivity']['valid_min']=0.0
    assert pyart.retrieve.grid_displacement_pc(test_storm,
            test_storm_2, 'reflectivity', 0) == (3,5)

@skipif(not pyart.retrieve._ADVECTION_AVAILABLE)
def test_grid_shift():
    #create two guassian storms
    tgrid0 = pyart.testing.make_normal_storm(10.0, [0.0,0.0])
    tgrid1 = pyart.testing.make_normal_storm(10.0, [5.0,5.0])
    #trim one, trim and shift the other
    tgrid1_reduced = pyart.retrieve.grid_shift(tgrid1, [0.0, 0.0],
            trim_edges = 10)
    tgrid0_shifted = pyart.retrieve.grid_shift(tgrid0, [5.0, 5.0],
            trim_edges = 10)
    #take the difference
    diff = tgrid0_shifted.fields['reflectivity']['data'][0,:,:]\
            - tgrid1_reduced.fields['reflectivity']['data'][0,:,:]
    #Assert that the difference is basically digital noise.
    #This actual value on my MB Air is -1.65633898546e-21
    assert diff.mean() < 1.0e-10

@skipif(not pyart.retrieve._ADVECTION_AVAILABLE)
def test_add_grids():
    tgrid0 = pyart.testing.make_normal_storm(10.00, [0.0,0.0])
    tgrid1 = pyart.testing.make_normal_storm(10.00, [0.0,0.0])
    rgrid = pyart.testing.make_normal_storm(10.0, [0.0,0.0])
    sgrid = pyart.retrieve.add_grids([tgrid0, tgrid1], [1.0, 1.0])
    image_resultant = rgrid.fields['reflectivity']['data']
    image_test = sgrid.fields['reflectivity']['data']
    assert np.abs(image_resultant - image_test).mean() < 1.0e-6




