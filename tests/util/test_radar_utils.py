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


# read in example file
radar = pyart.io.read_nexrad_archive(pyart.testing.NEXRAD_ARCHIVE_MSG31_FILE)
def test_image_mute_radar():

    # image mute example file
    mute_radar = pyart.util.image_mute_radar(
        radar, field='reflectivity', mute_field='cross_correlation_ratio',
        mute_threshold=0.97)

    # check that fields are added to the radar object
    assert 'nonmuted_reflectivity' in mute_radar.fields.keys()
    assert 'muted_reflectivity' in mute_radar.fields.keys()

    # check that the number of points in the muted and non muted reflectivity
    # fields have the same number of points as the rhoHV field
    n_rhohv = np.sum(~mute_radar.fields['cross_correlation_ratio']['data'].mask)
    n_mutez = np.sum(~mute_radar.fields['muted_reflectivity']['data'].mask)
    n_nonmutez = np.sum(~mute_radar.fields['nonmuted_reflectivity']['data'].mask)

    assert n_mutez + n_nonmutez == n_rhohv
