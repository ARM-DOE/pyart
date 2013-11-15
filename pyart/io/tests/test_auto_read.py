""" Unit Tests for Py-ART's io/mdv.py module. """

import tempfile
import os
import bz2
from StringIO import StringIO

from numpy.testing.decorators import skipif
from numpy.testing import assert_raises

import pyart


def test_autoread_mdv():
    radar = pyart.io.read(pyart.testing.MDV_PPI_FILE)
    source = radar.metadata['source']
    assert source == 'MDV radar volume file created by Dsr2Vol.'

    radar = pyart.io.read(pyart.testing.MDV_RHI_FILE)
    source = radar.metadata['source']
    assert source == 'MDV radar volume file created by Dsr2Vol.'


def test_autoread_sigmet():
    radar = pyart.io.read(pyart.testing.SIGMET_PPI_FILE, use_rsl=False)
    assert radar.metadata['original_container'] == 'sigmet'

    radar = pyart.io.read(pyart.testing.SIGMET_RHI_FILE, use_rsl=False)
    assert radar.metadata['original_container'] == 'sigmet'


@skipif(not pyart.io._RSL_AVAILABLE)
def test_autoread_sigmet_rsl():
    radar = pyart.io.read(pyart.testing.SIGMET_PPI_FILE, use_rsl=True)
    assert radar.metadata['original_container'] == 'rsl'

    radar = pyart.io.read(pyart.testing.SIGMET_RHI_FILE, use_rsl=True)
    assert radar.metadata['original_container'] == 'rsl'


def test_autoread_cfradial():
    radar = pyart.io.read(pyart.testing.CFRADIAL_PPI_FILE)
    assert radar.metadata['comment'] == 'none'

    radar = pyart.io.read(pyart.testing.CFRADIAL_RHI_FILE)
    assert radar.metadata['comment'] == 'none'


def test_autoread_nexrad_archive():
    radar = pyart.io.read(pyart.testing.NEXRAD_ARCHIVE_COMPRESSED_FILE)
    assert radar.metadata['original_container'] == 'NEXRAD Level II'

    radar = pyart.io.read(pyart.testing.NEXRAD_ARCHIVE_FILE)
    assert radar.metadata['original_container'] == 'NEXRAD Level II'


def test_autoread_nexrad_cdm():
    tmpfile = tempfile.mkstemp(suffix='.nc', dir='.')[1]
    f = open(tmpfile, 'wb')
    f.write(bz2.BZ2File(pyart.testing.NEXRAD_CDM_FILE).read())
    f.close()
    radar = pyart.io.read(tmpfile)
    os.remove(tmpfile)
    assert radar.metadata['original_container'] == 'NEXRAD Level II'


def test_autoread_raises():
    f = StringIO('0000000000000000000')
    assert_raises(TypeError, pyart.io.read, f)


def test_determine_filetype():
    headers = [
        ('\x00\x00\x03\xf8\x00\x007>\x00\x00\x00\x01', 'MDV'),
        ('\x89HDF\r\n\x1a\n\x02\x08\x08\x00', 'NETCDF4'),
        ('CDF', 'NETCDF3'),
        ('AR2V0006.501', 'WSR88D'),
        ('\x1b\x00\x08\x00\x00\x08\xb7\x07\x00\x00\x00\x00', 'SIGMET'),
        ('BZh91AY&SY\xbd\x12', 'BZ2'),
        ('UF', 'UF'),                   # not from a real file
        ('SSWB', 'DORADE'),             # not from a real file
        ('RSL', 'RSL'),                 # not from a real file
        ('\x0e\x03\x13\x01', 'HDF4'),   # not from a real file
        ('000000000000', 'UNKNOWN'),
    ]
    for i in headers:
        string, filetype = i
        check_filetype.description = 'determine filetype: %s' % (filetype)
        yield check_filetype, string, filetype


def check_filetype(string, filetype):
    f = StringIO(string)
    assert pyart.io.auto_read.determine_filetype(f) == filetype
