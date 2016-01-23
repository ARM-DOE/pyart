""" Unit Tests for Py-ART's util/radar_utils.py module. """

import numpy as np
import pyart


def test_is_vpt():
    radar = pyart.testing.make_empty_ppi_radar(10, 36, 3)
    assert not pyart.util.is_vpt(radar)
    pyart.util.to_vpt(radar)
    assert pyart.util.is_vpt(radar)


def test_to_vpt():
    # single scan
    radar = pyart.testing.make_empty_ppi_radar(10, 36, 3)
    radar.instrument_parameters = {
        'prt_mode': {'data': np.array(['fixed'] * 3)}
    }
    pyart.util.to_vpt(radar)
    assert pyart.util.is_vpt(radar)
    assert radar.nsweeps == 1
    assert radar.azimuth['data'][10] == 0.0
    assert radar.elevation['data'][0] == 90.0
    assert len(radar.instrument_parameters['prt_mode']['data']) == 1

    # multiple scans
    radar = pyart.testing.make_empty_ppi_radar(10, 36, 3)
    radar.instrument_parameters = {
        'prt_mode': {'data': np.array(['fixed'] * 3)}
    }
    pyart.util.to_vpt(radar, False)
    assert pyart.util.is_vpt(radar)
    assert radar.nsweeps == 108
    assert radar.azimuth['data'][10] == 10.0
    assert radar.elevation['data'][0] == 90.0
    assert len(radar.instrument_parameters['prt_mode']['data']) == 108