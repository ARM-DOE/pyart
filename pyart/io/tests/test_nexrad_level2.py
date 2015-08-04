""" Unit Tests for Py-ART's io/nexrad_level2.py module. """

import datetime
import tempfile
import os

import numpy as np
from numpy.testing import assert_array_equal, assert_raises
from numpy.testing import assert_almost_equal

import pyart
import pyart.io.nexrad_level2 as nexrad_level2

# create a NEXRADLevel2File from a UNCOMPRESSED dummy file
# pyart/testing/data/example_nexrad_archive.bz2
UNCOMPRESSED_FILE = pyart.testing.NEXRAD_ARCHIVE_FILE
nfile = nexrad_level2.NEXRADLevel2File(UNCOMPRESSED_FILE, bzip=True)


# attributes
def test_msg31s():
    assert len(nfile.msg31s) == 7200


def test_nscans():
    assert nfile.nscans == 16


def test_scan_msgs():
    assert len(nfile.scan_msgs) == 16
    assert len(nfile.scan_msgs[0]) == 720
    assert len(nfile.scan_msgs[-1]) == 360
    assert nfile.scan_msgs[0][10] == 10
    assert nfile.scan_msgs[10][10] == 5050


def test_volume_header():
    assert nfile.volume_header['date'] == 15904
    assert nfile.volume_header['extension'] == b'501'
    assert nfile.volume_header['icao'] == b'KATX'
    assert nfile.volume_header['tape'] == b'AR2V0006.'
    assert nfile.volume_header['time'] == 71424000


# methods
def test_get_azimuth_angles():
    angles = nfile.get_azimuth_angles([0, 1])
    assert angles.shape == (1440, )
    assert round(angles[0]) == 350.0
    assert round(angles[10]) == 355.0


def test_get_data_ref():
    data = nfile.get_data('REF', 1832, [0, 1])
    assert data.shape == (1440, 1832)
    assert round(data[0, 0]) == -32.0
    assert np.all(data == -32)
    assert np.ma.is_masked(data[-1, -1])


def test_get_data_vel():
    data = nfile.get_data('VEL', 1192, [1, 3])
    assert data.shape == (1440, 1192)
    assert round(data[0, 0], 2) == -63.5
    assert np.all(data == -63.5)


def test_get_data_sw():
    data = nfile.get_data('SW', 1192, [1, 3])
    assert data.shape == (1440, 1192)
    assert round(data[0, 0], 2) == -63.5
    assert np.all(data == -63.5)


def test_get_data_zdr():
    data = nfile.get_data('ZDR', 1192, [0, 2])
    assert data.shape == (1440, 1192)
    assert round(data[0, 0], 3) == -7.875
    assert np.all(data == -7.875)


def test_get_data_phi():
    data = nfile.get_data('PHI', 1192, [0, 2])
    assert data.shape == (1440, 1192)
    assert round(data[0, 0], 1) == 180.5


def test_get_data_rho():
    data = nfile.get_data('RHO', 1192, [0, 2])
    assert data.shape == (1440, 1192)
    assert_almost_equal(data[0, 0], 0.208, 3)


def test_get_data_ref_raw():
    data = nfile.get_data('REF', 1832, [0, 1], True)
    assert data.shape == (1440, 1832)
    assert round(data[0, 0]) == 2
    assert np.all(data[:720, :] == 2)
    assert np.all(data[720:, :1192] == 2)
    assert data[-1, -1] == 1        # folded


def test_get_data_vel_raw():
    data = nfile.get_data('VEL', 1192, [1, 3], True)
    assert data.shape == (1440, 1192)
    assert round(data[0, 0], 2) == 2
    assert np.all(data == 2)


def test_get_data_sw_raw():
    data = nfile.get_data('SW', 1192, [1, 3], True)
    assert data.shape == (1440, 1192)
    assert round(data[0, 0]) == 2
    assert np.all(data == 2)


def test_get_data_zdr_raw():
    data = nfile.get_data('ZDR', 1192, [0, 2], True)
    assert data.shape == (1440, 1192)
    assert round(data[0, 0]) == 2
    assert np.all(data == 2)


def test_get_data_phi_raw():
    data = nfile.get_data('PHI', 1192, [0, 2], True)
    assert data.shape == (1440, 1192)
    assert round(data[0, 0]) == 514     # big-endiann packing of 0,0
    assert np.all(data == 514)


def test_get_data_rho_raw():
    data = nfile.get_data('RHO', 1192, [0, 2], True)
    assert data.shape == (1440, 1192)
    assert round(data[0, 0]) == 2
    assert np.all(data == 2)


def test_get_elevation_angles():
    angles = nfile.get_elevation_angles([0, 1])
    assert angles.shape == (1440, )
    assert round(angles[0], 3) == 0.747
    assert round(angles[10], 3) == 0.670


def test_get_nrays():
    assert nfile.get_nrays(0) == 720
    assert nfile.get_nrays(4) == 360


def test_get_range():
    _range = nfile.get_range(0, 'REF')
    assert_array_equal(_range, 2125 + 250 * np.arange(1832))
    _range = nfile.get_range(4, 'REF')
    assert_array_equal(_range, 2125 + 250 * np.arange(1352))


def test_get_times():
    time_start, times = nfile.get_times([0, 1])
    assert time_start == datetime.datetime(2013, 7, 17, 19, 50, 21)
    assert times.shape == (1440, )
    assert round(times[0], 2) == 0.65
    assert round(times[10], 2) == 0.92


def test_location():
    lat, lon, alt = nfile.location()
    assert round(lat, 1) == 48.2
    assert round(lon, 1) == -122.5
    assert round(alt) == 195


def test_scan_info():
    scan_info = nfile.scan_info()

    assert len(scan_info) == 16

    nrays = [s['nrays'] for s in scan_info]
    ref_nrays = [720, 720, 720, 720, 360, 360, 360, 360, 360, 360, 360, 360,
                 360, 360, 360, 360]
    assert nrays == ref_nrays

    assert scan_info[0]['moments'] == ['REF', 'ZDR', 'PHI', 'RHO']
    assert scan_info[1]['moments'] == ['REF', 'VEL', 'SW']
    assert scan_info[2]['moments'] == ['REF', 'ZDR', 'PHI', 'RHO']
    assert scan_info[3]['moments'] == ['REF', 'VEL', 'SW']
    all_moments = ['REF', 'VEL', 'SW', 'ZDR', 'PHI', 'RHO']
    assert scan_info[4]['moments'] == all_moments
    assert scan_info[5]['moments'] == all_moments
    assert scan_info[6]['moments'] == all_moments
    assert scan_info[7]['moments'] == all_moments
    assert scan_info[8]['moments'] == all_moments
    assert scan_info[9]['moments'] == all_moments
    assert scan_info[10]['moments'] == all_moments
    assert scan_info[11]['moments'] == all_moments
    assert scan_info[12]['moments'] == all_moments
    assert scan_info[13]['moments'] == all_moments
    assert scan_info[14]['moments'] == all_moments
    assert scan_info[15]['moments'] == all_moments

    assert scan_info[0]['ngates'] == [1832, 1192, 1192, 1192]
    assert scan_info[1]['ngates'] == [1192, 1192, 1192]
    assert scan_info[2]['ngates'] == [1676, 1192, 1192, 1192]
    assert scan_info[3]['ngates'] == [1192, 1192, 1192]
    assert scan_info[4]['ngates'] == [1352, 1192, 1192, 1192, 1192, 1192]
    assert scan_info[5]['ngates'] == [1112, 1112, 1112, 1112, 1112, 1112]
    assert scan_info[6]['ngates'] == [940, 940, 940, 940, 940, 940]
    assert scan_info[7]['ngates'] == [800, 800, 800, 800, 800, 800]
    assert scan_info[8]['ngates'] == [704, 704, 704, 704, 704, 704]
    assert scan_info[9]['ngates'] == [540, 540, 540, 540, 540, 540]
    assert scan_info[10]['ngates'] == [500, 500, 500, 500, 500, 500]
    assert scan_info[11]['ngates'] == [460, 460, 460, 460, 460, 460]
    assert scan_info[12]['ngates'] == [388, 388, 388, 388, 388, 388]
    assert scan_info[13]['ngates'] == [332, 332, 332, 332, 332, 332]
    assert scan_info[14]['ngates'] == [280, 280, 280, 280, 280, 280]
    assert scan_info[15]['ngates'] == [240, 240, 240, 240, 240, 240]


# create a NEXRADLevel2File from a COMPRESSED file
# pyart/testing/data/example_nexrad_archive_compressed.ar2v
COMPRESSED_FILE = pyart.testing.NEXRAD_ARCHIVE_COMPRESSED_FILE
cfile = nexrad_level2.NEXRADLevel2File(COMPRESSED_FILE)


def test_compressed_attributes():
    # the compressed archive only contains the first 120 radials
    assert len(cfile.msg31s) == 120
    assert cfile.nscans == 1
    assert len(cfile.scan_msgs) == 1
    assert len(cfile.scan_msgs[0]) == 120
    assert cfile.volume_header['date'] == 15904
    assert cfile.volume_header['extension'] == b'501'
    assert cfile.volume_header['icao'] == b'KATX'
    assert cfile.volume_header['tape'] == b'AR2V0006.'
    assert cfile.volume_header['time'] == 71424000


# methods
def test_compressed_get_azimuth_angles():
    angles = cfile.get_azimuth_angles([0])
    assert angles.shape == (120, )
    assert round(angles[0]) == 350.0
    assert round(angles[10]) == 355.0


def test_compressed_get_data():
    data = cfile.get_data('REF', 1832, [0])
    assert data.shape == (120, 1832)
    assert round(data[0, 0], 1) == 10.5

    data = cfile.get_data('PHI', 1192, [0])
    assert data.shape == (120, 1192)
    assert round(data[0, 0]) == 254

    data = cfile.get_data('RHO', 1192, [0])
    assert data.shape == (120, 1192)
    assert_almost_equal(data[0, 0], 0.688, 3)

    data = cfile.get_data('ZDR', 1192, [0])
    assert data.shape == (120, 1192)
    assert round(data[0, 0], 3) == -7.875


def test_compressed_get_data_raw():
    data = cfile.get_data('REF', 1832, [0], True)
    assert data.shape == (120, 1832)
    assert round(data[0, 0]) == 87

    data = cfile.get_data('PHI', 1192, [0], True)
    assert data.shape == (120, 1192)
    assert round(data[0, 0]) == 722

    data = cfile.get_data('RHO', 1192, [0], True)
    assert data.shape == (120, 1192)
    assert round(data[0, 0]) == 146

    data = cfile.get_data('ZDR', 1192, [0], True)
    assert data.shape == (120, 1192)
    assert round(data[0, 0]) == 2


def test_compressed_get_elevation_angles():
    angles = cfile.get_elevation_angles([0])
    assert angles.shape == (120, )
    assert round(angles[0], 3) == 0.747
    assert round(angles[10], 3) == 0.670


def test_compressed_get_nrays():
    assert cfile.get_nrays(0) == 120


def test_compressed_get_range():
    _range = cfile.get_range(0, 'REF')
    assert_array_equal(_range, 2125 + 250 * np.arange(1832))


def test_compressed_get_times():
    time_start, times = cfile.get_times([0])
    assert time_start == datetime.datetime(2013, 7, 17, 19, 50, 21)
    assert times.shape == (120, )
    assert round(times[0], 2) == 0.65
    assert round(times[10], 2) == 0.92


def test_compressed_location():
    lat, lon, alt = cfile.location()
    assert round(lat, 1) == 48.2
    assert round(lon, 1) == -122.5
    assert round(alt) == 195


def test_compressed_scan_info():
    scan_info = cfile.scan_info()

    assert len(scan_info) == 1
    assert scan_info[0]['moments'] == ['REF', 'ZDR', 'PHI', 'RHO']
    assert scan_info[0]['ngates'] == [1832, 1192, 1192, 1192]
    assert scan_info[0]['nrays'] == 120


def test_bad_compression_header():

    # read the beginning of the compressed file
    f = open(COMPRESSED_FILE, 'rb')
    head = f.read(36)
    f.close()

    # corrupt the compression header and write to disk
    head = head[:28] + b'XX' + head[30:]
    tmpfile = tempfile.mkstemp(dir='.')[1]
    corrupt_file = open(tmpfile, 'wb')
    corrupt_file.write(head)
    corrupt_file.close()

    # should raise IOError
    assert_raises(IOError, nexrad_level2.NEXRADLevel2File, tmpfile)
    os.remove(tmpfile)
