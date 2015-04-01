""" Unit Tests for Py-ART's retrieve/simple_moment_calculation.py module. """

import numpy as np
import pyart


def test_calculate_snr_from_reflectivity():
    test_radar = pyart.testing.make_empty_ppi_radar(100, 360, 5)
    test_radar.range['data'] = test_radar.range['data'] * 100.0
    range_grid = np.meshgrid(
        test_radar.range['data'],
        np.ma.ones(test_radar.time['data'].shape))[0] + 1.0
    foo_field = {'data': np.zeros(
        [360 * 5, 100]) + 20.0 * np.log10(range_grid / 1000.)}
    test_radar.add_field('reflectivity', foo_field)
    snr = pyart.retrieve.calculate_snr_from_reflectivity(test_radar, toa=500)
    assert snr['data'].mean() < 1e-6
