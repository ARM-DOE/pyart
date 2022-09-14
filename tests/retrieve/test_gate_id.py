""" Unit Tests for Py-ART's retrieve/gate_id.py module. """

import numpy as np
import netCDF4

import pyart


def test_map_profile_to_gates():
    test_radar = pyart.testing.make_empty_ppi_radar(100, 360, 5)
    foo_field = {'data': np.zeros([360 * 5, 100])}
    test_radar.add_field('foo', foo_field)
    temp_dict = pyart.retrieve.map_profile_to_gates(
        np.ones(100), np.linspace(0, 1000, 100), test_radar)[1]
    assert temp_dict['data'].mean() == 1.0


def test_fetch_radar_time_profile():
    test_radar = pyart.testing.make_empty_ppi_radar(100, 360, 5)
    test_radar.time['units'] = 'seconds since 2011-05-10T00:00:01Z'
    test_radar.time['data'][0] = 41220.     # 2nd time in interpolated sonde

    sonde_dset = netCDF4.Dataset(pyart.testing.INTERP_SOUNDE_FILE)

    dic = pyart.retrieve.fetch_radar_time_profile(sonde_dset, test_radar)
    assert 'wdir' in dic
    assert 'wspd' in dic
    assert 'height' in dic
    assert round(dic['wdir'][0]) == 185
