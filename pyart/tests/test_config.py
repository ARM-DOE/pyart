""" Unit Tests for Py-ART's config.py module. """

import os
import warnings

import pyart

try:
    from importlib import reload
except ImportError:
    from imp import reload

dirname = os.path.dirname(__file__)
CUSTOM_CONFIG_FILE = os.path.join(dirname, 'custom_config.py')


def test_config_functions():
    pyart.load_config(CUSTOM_CONFIG_FILE)

    assert pyart.config.get_field_name('reflectivity') == 'REF'
    assert pyart.config.get_fillvalue() == -5555.0

    metadata = pyart.config.get_metadata('azimuth')
    assert isinstance(metadata, dict)
    assert metadata['units'] == 'foo'

    assert pyart.config.get_metadata('foobar') == {}

    pyart.load_config()     # load default
    assert pyart.config.get_field_name('reflectivity') == 'reflectivity'


def test_filemetadata_custom():

    pyart.load_config(CUSTOM_CONFIG_FILE)

    # The MDV reflectivity field is mapped to velocity
    filemetadata = pyart.config.FileMetadata('mdv')
    assert filemetadata.get_field_name('DBZ_F') == 'velocity'

    # there is no foobar field
    assert filemetadata.get_field_name('foobar') is None

    # and has a custom time metadata
    time = filemetadata('time')
    assert isinstance(time, dict)
    assert 'foo' in time
    assert time['foo'] == 'bar'

    # and standard elevation metadata
    elev = filemetadata('elevation')
    assert isinstance(elev, dict)
    assert 'units' in elev
    assert elev['units'] == 'degrees'

    assert filemetadata('foobar') == {}

    # invalid filetype does no field name mapping
    filemetadata = pyart.config.FileMetadata('foo')
    assert filemetadata.get_field_name('foobar') == 'foobar'

    # additional metadata
    filemetadata = pyart.config.FileMetadata(
        'sigmet', additional_metadata={'baz': {'units': 'baz_unit'}})
    assert filemetadata('baz')['units'] == 'baz_unit'

    # file_file_names
    filemetadata = pyart.config.FileMetadata('sigmet', file_field_names=True)
    assert filemetadata.get_field_name('DBT') == 'DBT'

    # exclude fields
    filemetadata = pyart.config.FileMetadata(
        'sigmet', exclude_fields=['spectrum_width'])
    assert filemetadata.get_field_name('WIDTH2') is None


def test_init_load():

    # load the custom config and verify
    os.environ['PYART_CONFIG'] = CUSTOM_CONFIG_FILE
    reload(pyart.config)
    assert pyart.config.get_field_name('reflectivity') == 'REF'

    # starting up with an invalid config file raises a warning and
    # loads the default config
    os.environ['PYART_CONFIG'] = 'nullnullnull'
    with warnings.catch_warnings(record=True) as w:
        reload(pyart.config)
        assert len(w) > 0
    assert pyart.config.get_field_name('reflectivity') == 'reflectivity'

    # load the custom config and verify
    os.environ['PYART_CONFIG'] = CUSTOM_CONFIG_FILE
    reload(pyart.config)
    assert pyart.config.get_field_name('reflectivity') == 'REF'

    # no PYART_CONFIG loads the default
    os.environ.pop('PYART_CONFIG')
    reload(pyart.config)
    assert pyart.config.get_field_name('reflectivity') == 'reflectivity'


def test_intergration():

    pyart.load_config()     # load default
    radar = pyart.io.read_mdv(pyart.testing.MDV_PPI_FILE)
    assert 'reflectivity' in radar.fields
    assert 'velocity' not in radar.fields

    # MDV custom maps reflectivity to velocity, and has custom time metadata
    pyart.load_config(CUSTOM_CONFIG_FILE)
    radar = pyart.io.read_mdv(pyart.testing.MDV_PPI_FILE)
    assert 'reflectivity' not in radar.fields
    assert 'velocity' in radar.fields
    assert radar.time['foo'] == 'bar'
