""" Unit Tests for Py-ART's io/mdv.py module. """


import numpy as np
import pytest

import pyart


@pytest.mark.skipif(not pyart.bridge.wradlib_bridge._WRADLIB_AVAILABLE,
                    reason="Wradlib is not installed.")
def test_texture_of_complex_phase():
    test_radar = pyart.testing.make_empty_ppi_radar(100, 360, 5)
    foo_field = {'data': np.zeros([360*5, 100])}
    test_radar.add_field('differential_phase', foo_field)
    test_text = pyart.retrieve.texture_of_complex_phase(test_radar)
    assert test_text['data'].mean() == 0.0
