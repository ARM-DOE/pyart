""" tests for py_mdv module of pyart.io """

from pyart.io import py_mdv
from datetime import datetime
import numpy as np
from os.path import join, dirname

mdv = py_mdv.read_mdv(join(dirname(__file__), '110635.mdv'))


def test_master_header():
    # test the master header
    ref_master_header = {
        'chunk_hdr_offset': 29824,
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
        'nfields': 20,
        'ngates': 983,
        'nrays': 360,
        'nsweeps': 17,
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
        'vlevel_hdr_offset': 9344,
        'vlevel_included': 1,
        'vlevel_type': 9}

    for k, v in ref_master_header.iteritems():
        assert mdv.master_header[k] == v


def test_field_header():
    # test a few of the field_header
    assert len(mdv.field_headers) == 20
    assert mdv.field_headers[0]['bad_data_value'] == 0.0
    assert mdv.field_headers[1]['dz_constant'] == 0
    assert mdv.field_headers[2]['encoding_type'] == 2
    assert mdv.field_headers[3]['field_name'] == "DBZ_F"
    assert mdv.field_headers[4]['nsweeps'] == 17


def test_vlevel_headers():
    # test the vlevel_headers
    assert len(mdv.vlevel_headers) == 20
    assert mdv.vlevel_headers[0]['record_len1'] == 1016
    assert mdv.vlevel_headers[0]['record_len2'] == 1016
    assert mdv.vlevel_headers[0]['unused_fl32'] == (0.0, 0.0, 0.0, 0.0, 0.0)
    assert mdv.vlevel_headers[0]['unused_si32'] == (0, 0, 0, 0)
    assert mdv.vlevel_headers[0]['struct_id'] == 14144
    assert len(mdv.vlevel_headers[0]['level']) == 122
    assert mdv.vlevel_headers[0]['level'][0] == 0.75
    assert mdv.vlevel_headers[0]['level'][10] == 11.699999809265137
    assert mdv.vlevel_headers[0]['level'][17] == 0.0
    assert len(mdv.vlevel_headers[0]['type']) == 122
    assert mdv.vlevel_headers[0]['type'][0:17] == (9,) * 17
    assert mdv.vlevel_headers[0]['type'][18] == 0


def test_chunk_headers():
    # test chunk_headers
    assert len(mdv.chunk_headers) == 3
    ref_chunk_header_0 = {
        'chunk_data_offset': 190457943,
        'chunk_id': 3,
        'info': 'DsRadar params',
        'record_len1': 504,
        'record_len2': 504,
        'size': 240,
        'struct_id': 14145,
        'unused_si32': (0, 0)}
    for k, v in ref_chunk_header_0.iteritems():
        assert mdv.chunk_headers[0][k] == v


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

    for k, v in ref_calib_info.iteritems():
        assert mdv.calib_info[k] == v


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
        'nfields': 62,
        'nfields_current': 0,
        'ngates': 983,
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

    for k, v in ref_radar_info.iteritems():
        assert mdv.radar_info[k] == v


def test_geometry():
    # test geometry attributes
    assert np.all(mdv.az_deg == np.arange(360))

    assert len(mdv.range_km) == 983
    assert mdv.range_km[10] == 1.3170482218265533

    assert len(mdv.el_deg) == 17
    assert mdv.el_deg[12] == 17.5


def test_time():
    # test times attribute
    assert mdv.times['time_begin'] == datetime(2011, 5, 20, 11, 1)
    assert mdv.times['time_centroid'] == datetime(2011, 5, 20, 11, 6, 35)
    assert mdv.times['time_end'] == datetime(2011, 5, 20, 11, 6, 35)


def test_cart():
    # test cartography information
    assert mdv.scan_type == 'ppi'
    assert mdv.carts['x'].shape == (17, 360, 983)
    assert mdv.carts['y'].shape == (17, 360, 983)
    assert mdv.carts['z'].shape == (17, 360, 983)
    assert round(mdv.carts['x'][1, 2, 3], 2) == 16.67
    assert round(mdv.carts['y'][1, 2, 3], 2) == 477.23
    assert round(mdv.carts['z'][1, 2, 3], 2) == 10.02


def test_fields():
    # test fields
    ref_fields = ['DBMHC', 'DBMVC', 'DBZ', 'DBZ_F', 'DBZVC', 'DBZVC_F',
                  'VEL', 'VEL_F', 'WIDTH', 'WIDTH_F', 'ZDR', 'ZDR_F',
                  'RHOHV', 'RHOHV_F', 'PHIDP', 'PHIDP_F', 'KDP', 'KDP_F',
                  'NCP', 'NCP_F']

    for a, b in zip(ref_fields, mdv.fields):
        assert a == b


def test_fileptr():
    # test fileptr
    assert type(mdv.fileptr) == file


def test_read_one_field():
    # extract a field
    assert hasattr(mdv, 'DBMHC') is False
    sweeps = mdv.read_a_field(0)
    assert sweeps.shape == (17, 360, 983)
    assert round(sweeps[1, 2, 3], 2) == -37.95
    assert hasattr(mdv, 'DBMHC') is True

    # compression level for last read field
    ref_current_compression_info = {
        'magic_cookie': 4160223223,
        'nbytes_coded': 502227,
        'nbytes_compressed': 502251,
        'nbytes_uncompressed': 707760,
        'spare': (0, 0)}

    for k, v in ref_current_compression_info.iteritems():
        assert mdv.current_compression_info[k] == v


def test_read_all_fields():

    # read all fields
    assert hasattr(mdv, 'PHIDP') is False
    mdv.get_all_fields()
    assert hasattr(mdv, 'DBMHC') is True
    assert hasattr(mdv, 'DBMVC') is True
    assert hasattr(mdv, 'DBZ') is True
    assert hasattr(mdv, 'DBZ_F') is True
    assert hasattr(mdv, 'DBZVC') is True
    assert hasattr(mdv, 'DBZVC_F') is True
    assert hasattr(mdv, 'VEL') is True
    assert hasattr(mdv, 'VEL_F') is True
    assert hasattr(mdv, 'WIDTH') is True
    assert hasattr(mdv, 'WIDTH_F') is True
    assert hasattr(mdv, 'ZDR') is True
    assert hasattr(mdv, 'ZDR_F') is True
    assert hasattr(mdv, 'RHOHV') is True
    assert hasattr(mdv, 'RHOHV_F') is True
    assert hasattr(mdv, 'PHIDP') is True
    assert hasattr(mdv, 'PHIDP_F') is True
    assert hasattr(mdv, 'KDP') is True
    assert hasattr(mdv, 'KDP_F') is True
    assert hasattr(mdv, 'NCP') is True
    assert hasattr(mdv, 'NCP_F') is True

    assert round(mdv.DBMHC[1, 2, 3], 2) == -37.95
    assert round(mdv.DBMVC[1, 2, 3], 2) == -40.63
    assert round(mdv.DBZ[1, 2, 3], 2) == 35.31
    assert round(mdv.DBZ_F[1, 2, 3], 2) == 17.6
    assert round(mdv.DBZVC[1, 2, 3], 2) == 32.64
    assert round(mdv.DBZVC_F[1, 2, 3], 2) == 15.85
    assert round(mdv.VEL[1, 2, 3], 2) == -0.93
    assert round(mdv.VEL_F[1, 2, 3], 2) == -0.18
    assert round(mdv.WIDTH[1, 2, 3], 2) == 0.63
    assert round(mdv.WIDTH_F[1, 2, 3], 2) == 7.26
    assert round(mdv.ZDR[1, 2, 3], 2) == 1.4
    assert round(mdv.ZDR_F[1, 2, 3], 2) == 0.47
    assert round(mdv.RHOHV[1, 2, 3], 2) == 0.98
    assert round(mdv.RHOHV_F[1, 2, 3], 2) == 0.72
    assert round(mdv.PHIDP[1, 2, 3], 2) == -155.72
    assert round(mdv.PHIDP_F[1, 2, 3], 2) == -155.72
    assert round(mdv.KDP[1, 2, 3], 2) == 0.1
    assert round(mdv.KDP_F[1, 2, 3], 2) == 0.1
    assert round(mdv.NCP[1, 2, 3], 2) == 0.98
    assert round(mdv.NCP_F[1, 2, 3], 2) == 0.39


test_read_all_fields.slow = True
