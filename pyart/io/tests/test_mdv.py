""" Unit Tests for Py-ART's io/mdv.py module. """

from datetime import datetime

import numpy as np
from numpy.ma.core import MaskedArray

import pyart

############################################
# read_mdv tests (verify radar attributes) #
############################################

# read in the sample file and create a a Radar object
radar = pyart.io.read_mdv(pyart.testing.MDV_PPI_FILE)


# time attribute
def test_time():
    assert 'comment' in radar.time.keys()
    assert 'long_name' in radar.time.keys()
    assert 'standard_name' in radar.time.keys()
    assert 'units' in radar.time.keys()
    assert 'calendar' in radar.time.keys()
    assert 'data' in radar.time.keys()
    assert radar.time['units'] == 'seconds since 2011-05-20T11:01:00Z'
    assert radar.time['data'].shape == (360, )
    assert round(radar.time['data'][200]) == 187.


# range attribute
def test_range():
    assert 'long_name' in radar.range
    assert 'standard_name' in radar.range
    assert 'meters_to_center_of_first_gate' in radar.range
    assert 'meters_between_gates' in radar.range
    assert 'units' in radar.range
    assert 'data' in radar.range
    assert 'spacing_is_constant' in radar.range
    assert radar.range['data'].shape == (110, )
    assert round(radar.range['data'][0]) == 118.0


# fields attribute is tested later


# metadata attribute
def test_metadata():
    assert 'instrument_name' in radar.metadata
    assert 'source' in radar.metadata


# scan_type attribute
def test_scan_type():
    assert radar.scan_type == 'ppi'


# latitude attribute
def test_latitude():
    assert 'data' in radar.latitude
    assert 'standard_name' in radar.latitude
    assert 'units' in radar.latitude
    assert radar.latitude['data'].shape == (1, )
    assert round(radar.latitude['data']) == 37.0


# longitude attribute
def test_longitude():
    assert 'data' in radar.longitude
    assert 'standard_name' in radar.longitude
    assert 'units' in radar.longitude
    assert radar.longitude['data'].shape == (1, )
    assert round(radar.longitude['data']) == -97.0


# altitude attribute
def test_altitude():
    assert 'data' in radar.altitude
    assert 'standard_name' in radar.altitude
    assert 'units' in radar.altitude
    assert 'positive' in radar.altitude
    assert radar.altitude['data'].shape == (1, )
    assert round(radar.altitude['data']) == 328.0


# altitude_agl attribute
def test_altitude_agl():
    assert radar.altitude_agl is None


# sweep_number attribute
def test_sweep_number():
    assert 'standard_name' in radar.sweep_number
    assert np.all(radar.sweep_number['data'] == range(1))


# sweep_mode attribute
def test_sweep_mode():
    assert 'standard_name' in radar.sweep_mode
    assert radar.sweep_mode['data'].shape == (1, )
    assert np.all(radar.sweep_mode['data'] == ['azimuth_surveillance'])


# fixed_angle attribute
def test_fixed_angle():
    assert 'standard_name' in radar.fixed_angle
    assert 'units' in radar.fixed_angle
    assert radar.fixed_angle['data'].shape == (1, )
    assert round(radar.fixed_angle['data'][0], 2) == 0.75


# sweep_start_ray_index attribute
def test_sweep_start_ray_index():
    assert 'long_name' in radar.sweep_start_ray_index
    assert radar.sweep_start_ray_index['data'].shape == (1, )
    assert round(radar.sweep_start_ray_index['data'][0]) == 0.0


# sweep_end_ray_index attribute
def test_sweep_end_ray_index():
    assert 'long_name' in radar.sweep_end_ray_index
    assert radar.sweep_end_ray_index['data'].shape == (1, )
    assert round(radar.sweep_end_ray_index['data'][0]) == 359.0


# target_scan_rate attribute
def test_target_scan_rate():
    assert radar.target_scan_rate is None


# azimuth attribute
def test_azimuth():
    assert 'standard_name' in radar.azimuth
    assert 'long_name' in radar.azimuth
    assert 'units' in radar.azimuth
    assert 'axis' in radar.azimuth
    assert round(radar.azimuth['data'][0]) == 0.0
    assert round(radar.azimuth['data'][10]) == 10.0


# elevation attribute
def test_elevation():
    assert 'standard_name' in radar.elevation
    assert 'long_name' in radar.azimuth
    assert 'units' in radar.elevation
    assert 'axis' in radar.elevation
    assert radar.elevation['data'].shape == (360, )
    assert round(radar.elevation['data'][0], 2) == 0.75


# scan_rate attribute
def test_scan_rate():
    assert radar.scan_rate is None


# antenna_transition attribute
def test_antenna_transition():
    assert radar.antenna_transition is None


# instrument_parameters attribute
def test_instument_parameters():
    # instrument_parameter sub-convention
    keys = ['prt', 'unambiguous_range', 'prt_mode', 'nyquist_velocity']
    for k in keys:
        description = 'instrument_parameters: %s' % k
        check_instrument_parameter.description = description
        yield check_instrument_parameter, k


def check_instrument_parameter(param):
    assert param in radar.instrument_parameters
    param_dic = radar.instrument_parameters[param]
    assert param_dic['meta_group'] == 'instrument_parameters'


# radar_calibration attribute
def test_radar_calibration():
    assert radar.radar_calibration is None


# ngates attribute
def test_ngates():
    assert radar.ngates == 110


# nrays attribute
def test_nrays():
    assert radar.nrays == 360


# nsweeps attribute
def test_nsweeps():
    assert radar.nsweeps == 1


####################
# fields attribute #
####################


def test_field_dics():
    fields = ['reflectivity_horizontal', ]
    for field in fields:
        description = "field : %s, dictionary" % field
        check_field_dic.description = description
        yield check_field_dic, field


def check_field_dic(field):
    """ Check that the required keys are present in a field dictionary. """
    assert 'standard_name' in radar.fields[field]
    assert 'units' in radar.fields[field]
    assert '_FillValue' in radar.fields[field]
    assert 'coordinates' in radar.fields[field]


def test_field_shapes():
    fields = ['reflectivity_horizontal', ]
    for field in fields:
        description = "field : %s, shape" % field
        check_field_shape.description = description
        yield check_field_shape, field


def check_field_shape(field):
    assert radar.fields[field]['data'].shape == (360, 110)


def test_field_types():
    fields = {'reflectivity_horizontal': MaskedArray, }
    for field, field_type in fields.iteritems():
        description = "field : %s, type" % field
        check_field_type.description = description
        yield check_field_type, field, field_type


def check_field_type(field, field_type):
    assert type(radar.fields[field]['data']) is field_type


def test_field_first_points():
    # these values can be found using:
    # [round(radar.fields[f]['data'][0,0]) for f in radar.fields]
    fields = {'reflectivity_horizontal': 24.0}
    for field, field_value in fields.iteritems():
        description = "field : %s, first point" % field
        check_field_first_point.description = description
        yield check_field_first_point, field, field_value


def check_field_first_point(field, value):
    assert round(radar.fields[field]['data'][0, 0]) == value


#################
# MdvFile tests #
#################


mdvfile = pyart.io.mdv.MdvFile(pyart.testing.MDV_PPI_FILE)


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
        'ngates': 110,
        'nrays': 360,
        'nsweeps': 1,
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

    for k, v in ref_master_header.iteritems():
        print "checking key:", k
        assert mdvfile.master_header[k] == v


def test_field_header():
    # test a few of the field_header
    assert len(mdvfile.field_headers) == 1
    assert mdvfile.field_headers[0]['bad_data_value'] == 0.0
    assert mdvfile.field_headers[0]['dz_constant'] == 0
    assert mdvfile.field_headers[0]['encoding_type'] == 2
    assert mdvfile.field_headers[0]['field_name'] == "DBZ_F"
    assert mdvfile.field_headers[0]['nsweeps'] == 1


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
    for k, v in ref_chunk_header_0.iteritems():
        print "checking key:", k
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

    for k, v in ref_calib_info.iteritems():
        print "checking key:", k
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

    for k, v in ref_radar_info.iteritems():
        print "checking key:", k
        assert mdvfile.radar_info[k] == v


def test_geometry():
    # test geometry attributes
    assert np.all(mdvfile.az_deg == np.arange(360))
    assert len(mdvfile.range_km) == 110
    assert round(mdvfile.range_km[10], 2) == 1.32
    assert len(mdvfile.el_deg) == 1
    assert mdvfile.el_deg[0] == 0.75


def test_mdv_time():
    # test times attribute
    assert mdvfile.times['time_begin'] == datetime(2011, 5, 20, 11, 1)
    assert mdvfile.times['time_centroid'] == datetime(2011, 5, 20, 11, 6, 35)
    assert mdvfile.times['time_end'] == datetime(2011, 5, 20, 11, 6, 35)


def test_cart():
    # test cartography information
    assert mdvfile.scan_type == 'ppi'
    assert mdvfile.carts['x'].shape == (1, 360, 110)
    assert mdvfile.carts['y'].shape == (1, 360, 110)
    assert mdvfile.carts['z'].shape == (1, 360, 110)
    assert round(mdvfile.carts['x'][0, 1, 2], 2) == 6.24
    assert round(mdvfile.carts['y'][0, 1, 2], 2) == 357.63
    assert round(mdvfile.carts['z'][0, 1, 2], 2) == 4.69


def test_fields():
    # test fields
    ref_fields = ['DBZ_F', ]
    for a, b in zip(ref_fields, mdvfile.fields):
        assert a == b


def test_fileptr():
    # test fileptr
    assert type(mdvfile.fileptr) == file


def test_read_one_field():
    # extract a field
    assert hasattr(mdvfile, 'DBZ_F') is False
    sweeps = mdvfile.read_a_field(0)
    assert sweeps.shape == (1, 360, 110)
    assert round(sweeps[0, 1, 2], 2) == 13.19
    assert hasattr(mdvfile, 'DBZ_F') is True

    # reread, should be faster
    sweeps = mdvfile.read_a_field(0)
    assert sweeps.shape == (1, 360, 110)
    assert round(sweeps[0, 1, 2], 2) == 13.19
    assert hasattr(mdvfile, 'DBZ_F') is True

    delattr(mdvfile, 'DBZ_F')


def test_read_all_fields():

    # read all fields
    assert hasattr(mdvfile, 'DBZ_F') is False
    mdvfile.read_all_fields()
    assert hasattr(mdvfile, 'DBZ_F') is True
    assert round(mdvfile.DBZ_F[0, 1, 2], 2) == 13.19


def test_read_all_fields_on_creation():
    mdvfile2 = pyart.io.mdv.MdvFile(pyart.testing.MDV_PPI_FILE,
                                    read_fields=True)
    assert hasattr(mdvfile2, 'DBZ_F') is True
    assert round(mdvfile2.DBZ_F[0, 1, 2], 2) == 13.19
    mdvfile2.close()


#############
# RHI tests #
#############

RADAR_RHI = pyart.io.read_mdv(pyart.testing.MDV_RHI_FILE)


# nsweeps attribute
def test_rhi_nsweeps():
    assert RADAR_RHI.nsweeps == 1


# sweep_number attribute
def test_rhi_sweep_number():
    assert 'standard_name' in RADAR_RHI.sweep_number
    assert np.all(RADAR_RHI.sweep_number['data'] == range(1))


# sweep_mode attribute
def test_rhi_sweep_mode():
    assert 'standard_name' in RADAR_RHI.sweep_mode
    assert RADAR_RHI.sweep_mode['data'].shape == (1, )
    assert np.all(RADAR_RHI.sweep_mode['data'] == ['rhi'])


# fixed_angle attribute
def test_rhi_fixed_angle():
    assert 'standard_name' in RADAR_RHI.fixed_angle
    assert 'units' in RADAR_RHI.fixed_angle
    assert RADAR_RHI.fixed_angle['data'].shape == (1, )
    assert round(RADAR_RHI.fixed_angle['data'][0], 2) == 189.


# sweep_start_ray_index attribute
def test_rhi_sweep_start_ray_index():
    assert 'long_name' in RADAR_RHI.sweep_start_ray_index
    assert RADAR_RHI.sweep_start_ray_index['data'].shape == (1, )
    assert round(RADAR_RHI.sweep_start_ray_index['data'][0]) == 0.0


# sweep_end_ray_index attribute
def test_rhi_sweep_end_ray_index():
    assert 'long_name' in RADAR_RHI.sweep_end_ray_index
    assert RADAR_RHI.sweep_end_ray_index['data'].shape == (1, )
    assert round(RADAR_RHI.sweep_end_ray_index['data'][0]) == 282.0


# azimuth attribute
def test_rhi_azimuth():
    assert 'standard_name' in RADAR_RHI.azimuth
    assert 'long_name' in RADAR_RHI.azimuth
    assert 'units' in RADAR_RHI.azimuth
    assert 'axis' in RADAR_RHI.azimuth
    assert round(RADAR_RHI.azimuth['data'][0]) == 189.0
    assert round(RADAR_RHI.azimuth['data'][10]) == 189.0


# elevation attribute
def test_rhi_elevation():
    assert 'standard_name' in RADAR_RHI.elevation
    assert 'long_name' in RADAR_RHI.azimuth
    assert 'units' in RADAR_RHI.elevation
    assert 'axis' in RADAR_RHI.elevation
    assert RADAR_RHI.elevation['data'].shape == (283, )
    assert round(RADAR_RHI.elevation['data'][0], 2) == 19.6
