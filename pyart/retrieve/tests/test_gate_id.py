""" Unit Tests for Py-ART's retrieve/gate_id.py module. """

import numpy as np
import pyart

def test_map_profile_to_gates():
    test_radar = pyart.testing.make_empty_ppi_radar(100, 360, 5)
    foo_field = {'data' : np.zeros([360*5, 100])}
    test_radar.add_field('foo' , foo_field)
    z_dict, temp_dict = pyart.retrieve.map_profile_to_gates(np.ones(100),
                                         np.linspace(0,1000,100),
                                         test_radar)
    assert temp_dict['data'].mean() == 1.0
