""" Unit Tests for Py-ART's io/mdv_common.py module. """

import warnings
from datetime import datetime
from io import BytesIO

import numpy as np
from numpy.testing import assert_almost_equal
import pytest

import pyart
from pyart.io.mdv_common import MdvFile

#################
# MdvFile tests #
#################


mdvfile = MdvFile(pyart.testing.MDV_PPI_FILE)


def test_master_header():
    # test the master header
    ref_master_header = {
        'chunk_hdr_offset': 2464,
        'data_collection_type': 0,
        'data_dimension': 0,
        'data_ordering': 0,
        'data_set_info': 'MDV radar volume file created by Dsr2Vol.',
        'data_set_name': 'C-SAPR',
        'data_set_source': 'ARM SGP C-SAPR',
        'field_grids_differ': 0,
        'field_hdr_offset': 1024,
        'grid_orientation': 1,
        'index_number': 611,
        'native_vlevel_type': 9,
        'nchunks': 3,
        'nfields': 1,
        'max_nx': 110,
        'max_ny': 360,
        'max_nz': 1,
        'num_data_times': 1,
        'record_len1': 1016,
        'record_len2': 1016,
        'revision_number': 1,
        'sensor_alt': 0.32760000228881836,
        'sensor_lat': 36.79615783691406,
        'sensor_lon': -97.45054626464844,
        'struct_id': 14142,
        'time_begin': 1305889260,
        'time_centroid': 1305889595,
        'time_end': 1305889595,
        'time_expire': 1305890265,
        'time_gen': 1305889595,
        'time_written': 1305889668,
        'unused_fl3212': (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0),
        'unused_si325': (0, 0, 0, 0, 0),
        'user_data': 0,
        'user_data_fl326': (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        'user_data_si328': (0, 0, 0, 0, 0, 0, 0, 0),
        'user_time': 0,
        'vlevel_hdr_offset': 1440,
        'vlevel_included': 1,
        'vlevel_type': 9}

    for k, v in ref_master_header.items():
        print("checking key:", k)
        assert mdvfile.master_header[k] == v


def test_field_header():
    # test a few of the field_header
    assert len(mdvfile.field_headers) == 1
    assert mdvfile.field_headers[0]['bad_data_value'] == 0.0
    assert mdvfile.field_headers[0]['dz_constant'] == 0
    assert mdvfile.field_headers[0]['encoding_type'] == 2
    assert mdvfile.field_headers[0]['field_name'] == "DBZ_F"
    assert mdvfile.field_headers[0]['nz'] == 1


def test_vlevel_headers():
    # test the vlevel_headers
    assert len(mdvfile.vlevel_headers) == 1
    assert mdvfile.vlevel_headers[0]['record_len1'] == 1016
    assert mdvfile.vlevel_headers[0]['record_len2'] == 1016
    assert mdvfile.vlevel_headers[0]['unused_fl32'] == (0.0, 0.0, 0.0,
                                                        0.0, 0.0)
    assert mdvfile.vlevel_headers[0]['unused_si32'] == (0, 0, 0, 0)
    assert mdvfile.vlevel_headers[0]['struct_id'] == 14144
    assert len(mdvfile.vlevel_headers[0]['level']) == 122
    assert mdvfile.vlevel_headers[0]['level'][0] == 0.75
    assert mdvfile.vlevel_headers[0]['level'][10] == 11.699999809265137
    assert mdvfile.vlevel_headers[0]['level'][17] == 0.0
    assert len(mdvfile.vlevel_headers[0]['type']) == 122
    assert mdvfile.vlevel_headers[0]['type'][0:17] == (9,) * 17
    assert mdvfile.vlevel_headers[0]['type'][18] == 0


def test_chunk_headers():
    # test chunk_headers
    assert len(mdvfile.chunk_headers) == 3
    ref_chunk_header_0 = {
        'chunk_data_offset': 68580,
        'chunk_id': 3,
        'info': 'DsRadar params',
        'record_len1': 504,
        'record_len2': 504,
        'size': 240,
        'struct_id': 14145,
        'unused_si32': (0, 0)}
    for k, v in ref_chunk_header_0.items():
        print("checking key:", k)
        assert mdvfile.chunk_headers[0][k] == v


def test_calib_info():
    # test calib_info
    ref_calib_info = {
        'I0_h_co_dbm': 0.0,
        'I0_h_cx_dbm': 0.0,
        'I0_v_co_dbm': 0.0,
        'I0_v_cx_dbm': 0.0,
        'antenna_gain_h_db': 44.5,
        'antenna_gain_v_db': 44.5,
        'beamwidth_h_deg': 1.0,
        'beamwidth_v_deg': 1.0,
        'coupler_fwd_loss_h_db': 50.0,
        'coupler_fwd_loss_v_db': 50.0,
        'day': 15,
        'filter_loss_db': 1.600000023841858,
        'hour': 18,
        'ldr_h_bias_db': 0.0,
        'ldr_v_bias_db': 0.0,
        'minute': 49,
        'month': 4,
        'noise_h_co_dbm': -73.5886001586914,
        'noise_h_cx_dbm': -73.5886001586914,
        'noise_source_h_dbm': -9999.0,
        'noise_source_v_dbm': -9999.0,
        'noise_v_co_dbm': -73.13050079345703,
        'noise_v_cx_dbm': -73.13050079345703,
        'power_meas_loss_h_db': 4.96999979019165,
        'power_meas_loss_v_db': 5.840000152587891,
        'pulse_width_us': 0.800000011920929,
        'radar_constant_h_db': 73.25920104980469,
        'radar_constant_v_db': 73.25920104980469,
        'radar_name': 'C-SAPR',
        'rx_gain_h_co_dbm': 33.31079864501953,
        'rx_gain_h_cx_dbm': 33.31079864501953,
        'rx_gain_v_co_dbm': 34.24689865112305,
        'rx_gain_v_cx_dbm': 34.24689865112305,
        'rx_slope_h_co_db': 1.0030399560928345,
        'rx_slope_h_cx_db': 1.0030399560928345,
        'rx_slope_v_co_db': 0.9963300228118896,
        'rx_slope_v_cx_db': 0.9963300228118896,
        'second': 36,
        'spare': (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0),
        'sun_h_co_dbm': -9999.0,
        'sun_h_cx_dbm': -9999.0,
        'sun_v_co_dbm': -9999.0,
        'sun_v_cx_dbm': -9999.0,
        'system_phidp_deg': -60.0,
        'test_pulse_h_dbm': -9999.0,
        'test_pulse_v_dbm': -9999.0,
        'twoway_radome_loss_h_db': 0.5,
        'twoway_radome_loss_v_db': 0.5,
        'twoway_waveguide_loss_h_db': 4.0,
        'twoway_waveguide_loss_v_db': 4.0,
        'wavelength_cm': 5.329639911651611,
        'xmit_power_h_dbm': 83.6500015258789,
        'xmit_power_v_dbm': 83.6500015258789,
        'year': 2011,
        'zdr_bias_db': -1.2799999713897705,
        'zh1km_co_dbz': -33.63949966430664,
        'zh1km_cx_dbz': -33.63949966430664,
        'zv1km_co_dbz': -34.11470031738281,
        'zv1km_cx_dbz': -34.11470031738281}

    for k, v in ref_calib_info.items():
        print("checking key:", k)
        assert mdvfile.calib_info[k] == v


def test_radar_info():

    # check radar_info
    ref_radar_info = {
        'altitude_km': 0.32760000228881836,
        'antenna_gain_db': 44.5,
        'field_flag': 0,
        'follow_mode': 1,
        'gate_spacing_km': 0.11991698294878006,
        'horiz_beam_width_deg': 1.0,
        'latitude_deg': 36.79615783691406,
        'longitude_deg': -97.45054626464844,
        'measXmitPowerDbmH_dbm': 82.7537612915039,
        'measXmitPowerDbmV_dbm': 82.65692138671875,
        'nfields': 1,
        'nfields_current': 0,
        'ngates': 110,
        'polarization': 5,
        'prf_hz': 1239.9993896484375,
        'prf_mode': 1,
        'prt2_s': 0.0,
        'prt_s': 0.0008064520079642534,
        'pulse_width_us': 0.800000011920929,
        'radar_constant': 73.25920104980469,
        'radar_id': 0,
        'radar_name': 'C-SAPR',
        'radar_type': 0,
        'receiver_gain_db': 33.31079864501953,
        'receiver_mds_dbm': -106.89939880371094,
        'samples_per_beam': 128,
        'scan_mode': 8,
        'scan_type': 0,
        'scan_type_name': 'Default_s',
        'spare_floats': (0.0, 0.0, 0.0, 0.0),
        'spare_ints': (0, 0),
        'start_range_km': 0.11787839233875275,
        'system_gain_db': 77.81079864501953,
        'unambig_range_km': 117.99627685546875,
        'unambig_vel_mps': 16.52467918395996,
        'vert_beam_width_deg': 1.0,
        'wavelength_cm': 5.330544471740723,
        'xmit_peak_pwr_watts': 231739552.0}

    for k, v in ref_radar_info.items():
        print("checking key:", k)
        assert mdvfile.radar_info[k] == v


def test_geometry():
    # test geometry attributes
    az_deg, range_km, el_deg = mdvfile._calc_geometry()
    assert np.all(az_deg == np.arange(360))
    assert len(range_km) == 110
    assert_almost_equal(range_km[10], 1.32, 2)
    assert len(el_deg) == 1
    assert el_deg[0] == 0.75


def test_geometry_raises():
    mdvfile = MdvFile(pyart.testing.MDV_GRID_FILE)
    pytest.raises(NotImplementedError, mdvfile._calc_geometry)


def test_mdv_time():
    # test times attribute
    assert mdvfile.times['time_begin'] == datetime(2011, 5, 20, 11, 1)
    assert mdvfile.times['time_centroid'] == datetime(2011, 5, 20, 11, 6, 35)
    assert mdvfile.times['time_end'] == datetime(2011, 5, 20, 11, 6, 35)


def test_cart():
    # test cartography information
    assert mdvfile.projection == 'ppi'
    carts = mdvfile._make_carts_dict()
    assert carts['x'].shape == (1, 360, 110)
    assert carts['y'].shape == (1, 360, 110)
    assert carts['z'].shape == (1, 360, 110)
    assert_almost_equal(carts['x'][0, 1, 2], 6.24, 2)
    assert_almost_equal(carts['y'][0, 1, 2], 357.63, 2)
    assert_almost_equal(carts['z'][0, 1, 2], 4.69, 2)


def test_fields():
    # test fields
    ref_fields = ['DBZ_F', ]
    for a, b in zip(ref_fields, mdvfile.fields):
        assert a == b


def test_fileptr():
    # test fileptr
    assert hasattr(mdvfile.fileptr, 'read')


def test_read_one_field():

    mdvfile = MdvFile(pyart.testing.MDV_PPI_FILE)
    # extract a field
    assert mdvfile.fields_data[0] is None
    sweeps = mdvfile.read_a_field(0)
    assert sweeps.shape == (1, 360, 110)
    assert_almost_equal(sweeps[0, 1, 2], 13.19, 2)
    assert mdvfile.fields_data[0] is not None

    # reread, should be faster
    sweeps = mdvfile.read_a_field(0)
    assert sweeps.shape == (1, 360, 110)
    assert_almost_equal(sweeps[0, 1, 2], 13.19, 2)
    assert mdvfile.fields_data[0] is not None


def test_read_all_fields():

    mdvfile = MdvFile(pyart.testing.MDV_PPI_FILE)
    # read all fields
    assert mdvfile.fields_data[0] is None
    mdvfile.read_all_fields()
    assert mdvfile.fields_data[0] is not None
    assert_almost_equal(mdvfile.fields_data[0][0, 1, 2], 13.19, 2)


def test_read_all_fields_on_creation():
    mdvfile2 = MdvFile(pyart.testing.MDV_PPI_FILE, read_fields=True)
    assert mdvfile2.fields_data[0] is not None
    assert_almost_equal(mdvfile2.fields_data[0][0, 1, 2], 13.19, 2)
    mdvfile2.close()


def test_mdvfile_radar_stubs():
    mdvfile = pyart.io.mdv_common.MdvFile(None)
    # These methods are included in MdvFile as a stub for a future
    # write_radar_mdv function, test that the dictionaries they return have the
    # correct length.
    dic = mdvfile._get_chunk_header()
    assert len(dic) == 8
    assert isinstance(dic, dict)

    dic = mdvfile._get_radar_info()
    assert len(dic) == 38
    assert isinstance(dic, dict)

    dic = mdvfile._get_calib()
    assert len(dic) == 59
    assert isinstance(dic, dict)

    dic = mdvfile._get_compression_info()
    assert len(dic) == 5
    assert isinstance(dic, dict)

    dic = mdvfile._get_levels_info(10)
    assert len(dic) == 2
    assert isinstance(dic, dict)


def test_rhi_cart():
    mdvfile = pyart.io.mdv_common.MdvFile(pyart.testing.MDV_RHI_FILE)
    assert mdvfile.projection == 'rhi'
    carts = mdvfile._make_carts_dict()
    assert carts['x'].shape == (1, 283, 125)
    assert carts['y'].shape == (1, 283, 125)
    assert carts['z'].shape == (1, 283, 125)
    assert_almost_equal(carts['x'][0, 1, 2], -52.63, 2)
    assert_almost_equal(carts['y'][0, 1, 2], -332.31, 2)
    assert_almost_equal(carts['z'][0, 1, 2], 121.47, 2)


class Mdv_common_Tests(object):
    """
    Class for declaring unit tests for the io/mdv_common.py module which
    require a temporary file which will be removed at the end of the test.
    """

    def test_mdv_file_read_write_radar(self):
        with warnings.catch_warnings(record=True) as w:
            mdvfile_orig = pyart.io.mdv_common.MdvFile(
                pyart.testing.MDV_PPI_FILE)
            mdvfile_orig.read_all_fields()

            inmemfile = BytesIO()
            mdvfile_orig.write(inmemfile)
            inmemfile.seek(0)
            # check that a UserWarning was issued since zlib compression used
            assert len(w) == 1
            assert issubclass(w[-1].category, UserWarning)

            mdvfile = pyart.io.mdv_common.MdvFile(inmemfile)
            self.check_mdvfile_ppi(mdvfile)

    def test_read_write_file_objects(self):
        # read from and write two a file object
        f = open(pyart.testing.MDV_PPI_FILE, 'rb')
        mdvfile = pyart.io.mdv_common.MdvFile(f)
        mdvfile.read_all_fields()
        self.check_mdvfile_ppi(mdvfile)

        # write out the file using an file handler
        f2 = BytesIO()
        with warnings.catch_warnings():
            # ignore the UserWarning about the non-implement compression_type
            warnings.simplefilter("ignore", category=UserWarning)
            mdvfile.write(f2)
        f.close()
        f2.seek(0)

        # re-read and check object
        mdvfile = pyart.io.mdv_common.MdvFile(f2)
        self.check_mdvfile_ppi(mdvfile)

    @staticmethod
    def check_mdvfile_ppi(mdvfile):
        # check some parameters
        assert mdvfile.master_header['time_end'] == 1305889595
        assert mdvfile.master_header['chunk_hdr_offset'] == 2464
        assert len(mdvfile.field_headers) == 1
        assert mdvfile.field_headers[0]['bad_data_value'] == 0.0
        assert len(mdvfile.vlevel_headers) == 1
        assert mdvfile.vlevel_headers[0]['record_len1'] == 1016
        assert mdvfile.chunk_headers[0]['chunk_id'] == 3
        assert mdvfile.calib_info['day'] == 15
        assert mdvfile.radar_info['scan_mode'] == 8
        az_deg, range_km, el_deg = mdvfile._calc_geometry()
        assert np.all(az_deg == np.arange(360))
        assert len(range_km) == 110
        assert_almost_equal(range_km[10], 1.32, 2)
        assert len(el_deg) == 1
        assert el_deg[0] == 0.75
        assert mdvfile.times['time_begin'] == datetime(2011, 5, 20, 11, 1)
        assert len(mdvfile.elevations) == 0

    def test_mdv_file_read_write_radar_rhi(self):
        # write and read the RHI file which contains a elevation chunk
        mdvfile = pyart.io.mdv_common.MdvFile(pyart.testing.MDV_RHI_FILE)
        mdvfile.read_all_fields()
        inmemfile = BytesIO()
        with warnings.catch_warnings():
            # ignore the UserWarning about the non-implement compression_type
            warnings.simplefilter("ignore", category=UserWarning)
            mdvfile.write(inmemfile)
        inmemfile.seek(0)
        mdvfile2 = pyart.io.mdv_common.MdvFile(inmemfile)

        # verify that the elevations are similar
        assert len(mdvfile2.elevations) == len(mdvfile.elevations)
        assert mdvfile2.elevations[0] == mdvfile.elevations[0]
