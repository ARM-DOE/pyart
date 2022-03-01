""" Unit Tests for Py-ART's io/mdv.py module. """

import bz2
from io import BytesIO

import pytest

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


@pytest.mark.skipif(not pyart.io.rsl._RSL_AVAILABLE,
                    reason="TRMM RSL is not installed.")
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
    radar = pyart.io.read(pyart.testing.NEXRAD_ARCHIVE_MSG31_COMPRESSED_FILE,
                          delay_field_loading=True)
    assert radar.metadata['original_container'] == 'NEXRAD Level II'

    radar = pyart.io.read(pyart.testing.NEXRAD_ARCHIVE_MSG31_FILE,
                          delay_field_loading=True)
    assert radar.metadata['original_container'] == 'NEXRAD Level II'


def test_autoread_nexrad_cdm():
    with pyart.testing.InTemporaryDirectory():
        tmpfile = 'tmp_nexrad.nc'
        with open(tmpfile, 'wb') as f:
            f.write(bz2.BZ2File(pyart.testing.NEXRAD_CDM_FILE).read())
        radar = pyart.io.read(tmpfile)
        assert radar.metadata['original_container'] == 'NEXRAD Level II'


def test_autoread_nexrad_level3():
    radar = pyart.io.read(pyart.testing.NEXRAD_LEVEL3_MSG19)
    assert radar.metadata['original_container'] == 'NEXRAD Level 3'


def test_autoread_raises():
    f = BytesIO(b'0000000000000000000')
    pytest.raises(TypeError, pyart.io.read, f)


headers = [
    (b'\x00\x00\x03\xf8\x00\x007>\x00\x00\x00\x01', 'MDV'),
    (b'\x89HDF\r\n\x1a\n\x02\x08\x08\x00', 'NETCDF4'),
    (b'CDF', 'NETCDF3'),
    (b'AR2V0006.501', 'WSR88D'),
    (b'SDUS54 KBMX ', 'NEXRADL3'),
    (b'\x1b\x00\x08\x00\x00\x08\xb7\x07\x00\x00\x00\x00', 'SIGMET'),
    (b'BZh91AY&SY\xbd\x12', 'BZ2'),
    (b'UF', 'UF'),                   # not from a real file
    (b'SSWB', 'DORADE'),             # not from a real file
    (b'RSL', 'RSL'),                 # not from a real file
    (b'\x0e\x03\x13\x01', 'HDF4'),   # not from a real file
    (b'000000000000', 'UNKNOWN'),
    ]
@pytest.mark.parametrize("i", headers)
def test_determine_filetype(i):
    string, filetype = i
    check_filetype.description = 'determine filetype: %s' % (filetype)
    check_filetype(string, filetype)


def check_filetype(string, filetype):
    f = BytesIO(string)
    assert pyart.io.auto_read.determine_filetype(f) == filetype
