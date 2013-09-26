""" Unit Tests for Py-ART's io/radar.py module. """

import numpy as np
from numpy.testing import assert_raises
import pyart


def test_radar_creation():
    radar = pyart.testing.make_target_radar()
    assert isinstance(radar, pyart.io.Radar)


def test_add_field():
    radar = pyart.testing.make_target_radar()
    dic = {'data': np.zeros((360, 50)), 'standard_name': 'test'}
    radar.add_field('test', dic)
    assert 'test' in radar.fields
    assert 'data' in radar.fields['test']
    assert radar.fields['test']['standard_name'] == 'test'


def test_add_field_errors():
    radar = pyart.testing.make_target_radar()

    assert_raises(ValueError, radar.add_field, 'reflectivity_horizontal', {})

    dic = {'dat': np.zeros((360, 50)), 'standard_name': 'test'}
    assert_raises(KeyError, radar.add_field, 'test', dic)

    dic = {'data': np.zeros((360, 49)), 'standard_name': 'test'}
    assert_raises(ValueError, radar.add_field, 'test', dic)


def test_add_field_like():
    radar = pyart.testing.make_target_radar()
    data = np.zeros((360, 50))
    radar.add_field_like('reflectivity_horizontal', 'test', data)
    assert 'test' in radar.fields
    assert 'data' in radar.fields['test']
    assert radar.fields['test']['units'] == 'dBZ'


def test_add_field_like_errors():
    radar = pyart.testing.make_target_radar()
    assert_raises(ValueError, radar.add_field_like, 'foo', 'bar', [])
