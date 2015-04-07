"""
pyart.io.mdv
============

Utilities for reading of MDV files.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    MdvFile

.. autosummary::
    :toctree: generated/

    read_mdv

"""
# Code is adapted from Nitin Bharadwaj's Matlab code

import struct
import gzip
import zlib
import StringIO
import datetime

import numpy as np
from netCDF4 import date2num

from ..config import FileMetadata, get_fillvalue
from ..core.radar import Radar
from .common import make_time_unit_str
from .common import radar_coords_to_cart


def read_mdv(filename, field_names=None, additional_metadata=None,
             file_field_names=False, exclude_fields=None):
    """
    Read a MDV file.

    Parameters
    ----------
    filename : str
        Name of MDV file to read or file-like object pointing to the
        beginning of such a file.
    field_names : dict, optional
        Dictionary mapping MDV data type names to radar field names. If a
        data type found in the file does not appear in this dictionary or has
        a value of None it will not be placed in the radar.fields dictionary.
        A value of None, the default, will use the mapping defined in the
        Py-ART configuration file.
    additional_metadata : dict of dicts, optional
        Dictionary of dictionaries to retrieve metadata from during this read.
        This metadata is not used during any successive file reads unless
        explicitly included.  A value of None, the default, will not
        introduct any addition metadata and the file specific or default
        metadata as specified by the Py-ART configuration file will be used.
    file_field_names : bool, optional
        True to use the MDV data type names for the field names. If this
        case the field_names parameter is ignored. The field dictionary will
        likely only have a 'data' key, unless the fields are defined in
        `additional_metadata`.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters.

    Returns
    -------
    radar : Radar
        Radar object containing data from MDV file.

    Notes
    -----
    Currently this function can only read polar MDV files with fields
    compressed with gzip or zlib.

    """
    # create metadata retrieval object
    filemetadata = FileMetadata('mdv', field_names, additional_metadata,
                                file_field_names, exclude_fields)

    mdvfile = MdvFile(filename)

    # value attributes
    az_deg, range_km, el_deg = mdvfile._calc_geometry()
    naz = len(az_deg)
    nele = len(el_deg)
    scan_type = mdvfile.projection

    if scan_type not in ['ppi', 'rhi']:
        raise NotImplementedError('No support for scan_type %s.' % scan_type)

    # time
    time = filemetadata('time')
    units = make_time_unit_str(mdvfile.times['time_begin'])
    time['units'] = units
    time_start = date2num(mdvfile.times['time_begin'], units)
    time_end = date2num(mdvfile.times['time_end'], units)
    time['data'] = np.linspace(time_start, time_end, naz * nele)

    # range
    _range = filemetadata('range')
    _range['data'] = np.array(range_km * 1000.0, dtype='float32')
    _range['meters_to_center_of_first_gate'] = _range['data'][0]
    _range['meters_between_gates'] = (_range['data'][1] - _range['data'][0])

    # fields
    fields = {}
    mdv_fields = mdvfile._make_fields_list()
    for mdv_field in set(mdv_fields):
        field_name = filemetadata.get_field_name(mdv_field)
        if field_name is None:
            continue

        # grab data from MDV object, mask and reshape
        data = mdvfile.read_a_field(mdv_fields.index(mdv_field))
        data[np.where(np.isnan(data))] = get_fillvalue()
        data[np.where(data == 131072)] = get_fillvalue()
        data = np.ma.masked_equal(data, get_fillvalue())
        data.shape = (data.shape[0] * data.shape[1], data.shape[2])

        # create and store the field dictionary
        field_dic = filemetadata(field_name)
        field_dic['data'] = data
        field_dic['_FillValue'] = get_fillvalue()
        fields[field_name] = field_dic

    # metadata
    metadata = filemetadata('metadata')
    for meta_key, mdv_key in MDV_METADATA_MAP.iteritems():
        metadata[meta_key] = mdvfile.master_header[mdv_key]

    # latitude
    latitude = filemetadata('latitude')
    latitude['data'] = np.array([mdvfile.radar_info['latitude_deg']],
                                dtype='float64')
    # longitude
    longitude = filemetadata('longitude')
    longitude['data'] = np.array([mdvfile.radar_info['longitude_deg']],
                                 dtype='float64')
    # altitude
    altitude = filemetadata('altitude')
    altitude['data'] = np.array([mdvfile.radar_info['altitude_km'] * 1000.0],
                                dtype='float64')

    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
    # sweep_end_ray_index
    sweep_number = filemetadata('sweep_number')
    sweep_mode = filemetadata('sweep_mode')
    fixed_angle = filemetadata('fixed_angle')
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')
    len_time = len(time['data'])

    if mdvfile.scan_type == 'ppi':
        nsweeps = nele
        sweep_number['data'] = np.arange(nsweeps, dtype='int32')
        sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])
        fixed_angle['data'] = np.array(el_deg, dtype='float32')
        sweep_start_ray_index['data'] = np.arange(0, len_time, naz,
                                                  dtype='int32')
        sweep_end_ray_index['data'] = np.arange(naz-1, len_time, naz,
                                                dtype='int32')

    elif mdvfile.scan_type == 'rhi':
        nsweeps = naz
        sweep_number['data'] = np.arange(nsweeps, dtype='int32')
        sweep_mode['data'] = np.array(nsweeps * ['rhi'])
        fixed_angle['data'] = np.array(az_deg, dtype='float32')
        sweep_start_ray_index['data'] = np.arange(0, len_time, nele,
                                                  dtype='int32')
        sweep_end_ray_index['data'] = np.arange(nele - 1, len_time, nele,
                                                dtype='int32')

    # azimuth, elevation
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')

    if scan_type == 'ppi':
        azimuth['data'] = np.tile(az_deg, nele)
        elevation['data'] = np.array(el_deg).repeat(naz)

    elif scan_type == 'rhi':
        azimuth['data'] = np.array(az_deg).repeat(nele)
        elevation['data'] = np.tile(el_deg, naz)

    # instrument parameters
    # we will set 4 parameters in the instrument_parameters dict
    # prt, prt_mode, unambiguous_range, and nyquist_velocity

    # TODO prt mode: Need to fix this.. assumes dual if two prts
    if mdvfile.radar_info['prt2_s'] == 0.0:
        prt_mode_str = 'fixed'
    else:
        prt_mode_str = 'dual'

    prt_mode = filemetadata('prt_mode')
    prt = filemetadata('prt')
    unambiguous_range = filemetadata('unambiguous_range')
    nyquist_velocity = filemetadata('nyquist_velocity')
    beam_width_h = filemetadata('radar_beam_width_h')
    beam_width_v = filemetadata('radar_beam_width_v')

    prt_mode['data'] = np.array([prt_mode_str] * nsweeps)
    prt['data'] = np.array([mdvfile.radar_info['prt_s']] * nele * naz,
                           dtype='float32')

    urange_m = mdvfile.radar_info['unambig_range_km'] * 1000.0
    unambiguous_range['data'] = np.array([urange_m] * naz * nele,
                                         dtype='float32')

    uvel_mps = mdvfile.radar_info['unambig_vel_mps']
    nyquist_velocity['data'] = np.array([uvel_mps] * naz * nele,
                                        dtype='float32')
    beam_width_h['data'] = np.array(
        [mdvfile.radar_info['horiz_beam_width_deg']], dtype='float32')
    beam_width_v['data'] = np.array(
        [mdvfile.radar_info['vert_beam_width_deg']], dtype='float32')

    instrument_parameters = {'prt_mode': prt_mode, 'prt': prt,
                             'unambiguous_range': unambiguous_range,
                             'nyquist_velocity': nyquist_velocity,
                             'radar_beam_width_h': beam_width_h,
                             'radar_beam_width_v': beam_width_v}

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)


# mapping from MDV name space to CF-Radial name space
MDV_METADATA_MAP = {'instrument_name': 'data_set_source',
                    'source': 'data_set_info'}

# Information about the MDV file structure
MDV_CHUNK_INFO_LEN = 480
MDV_INFO_LEN = 512
MDV_LONG_FIELD_LEN = 64
MDV_MAX_PROJ_PARAMS = 8
MDV_MAX_VLEVELS = 122
MDV_NAME_LEN = 128
MDV_SHORT_FIELD_LEN = 16
MDV_TRANSFORM_LEN = 16
MDV_UNITS_LEN = 16
MDV_N_COORD_LABELS = 3
MDV_COORD_UNITS_LEN = 32

# (x,y) in degrees. Simple latitude-longitude grid.
# Also known as the Simple Cylindrical or Platte Carree projection.
PROJ_LATLON = 0
# (x,y) in km. Lambert Conformal Conic projection.
PROJ_LAMBERT_CONF = 3
# (x,y) in km. Polar Stereographic projection.
PROJ_POLAR_STEREO = 5
# Cartesian, (x,y) in km. This is a simple line-of-sight
# projection used for single radar sites. The formal name is
# Oblique Lambert Azimuthal projection.
PROJ_FLAT = 8
# radar data in native Plan Position Indicator (PPI) coordinates of
# range, azimuth angle and elevation angle. x is radial range (km),
# y is azimuth angle (deg), z is elev angle (deg).
PROJ_POLAR_RADAR = 9
# (x,y) in km. Oblique Stereographic projection.
PROJ_OBLIQUE_STEREO = 12
# radar data in native Range Height Indicator (RHI) coordinates.
# x is radial range (km), y is elev angle (deg), z is az angle (deg).
PROJ_RHI_RADAR = 13
#  ***************** COMPRESSION *******************
COMPRESSION_NONE = 0  # no compression
COMPRESSION_ZLIB = 3  # Lempel-Ziv
COMPRESSION_BZIP = 4  # bzip2
COMPRESSION_GZIP = 5  # Lempel-Ziv in gzip format

#  ***************** COMPRESSION CODE *******************
TA_NOT_COMPRESSED = 791621423
GZIP_COMPRESSED = 4160223223

#  ***************** TRANSFORM *******************
DATA_TRANSFORM_NONE = 0  # None
DATA_TRANSFORM_LOG = 1  # Natural log

#  ***************** BIT ENCODING *******************
ENCODING_INT8 = 1  # unsigned 8 bit integer
ENCODING_INT16 = 2  # unsigned 16 bit integer
ENCODING_FLOAT32 = 5  # 32 bit IEEE floating point

#  ***************** CHUNK HEADER and DATA *******************
CHUNK_DSRADAR_PARAMS = 3
CHUNK_DSRADAR_ELEVATIONS = 10
CHUNK_DSRADAR_CALIB = 7
DS_LABEL_LEN = 40
NCHAR_DS_RADAR_PARAMS = 2 * DS_LABEL_LEN
DS_RADAR_CALIB_NAME_LEN = 16
DS_RADAR_CALIB_MISSING = -9999.0


class MdvFile:
    """
    A file object for MDV data.

    A `MdvFile` object stores metadata and data from a MDV file.  Metadata is
    stored in dictionaries as attributes of the object, field data is
    stored as NumPy ndarrays as attributes with the field name. By default
    only metadata is read initially and field data must be read using the
    `read_a_field` or `read_all_fields` methods.  This behavior can be changed
    by setting the `read_fields` parameter to True.

    Parameters
    ----------
    filename : str or file-like
        Name of MDV file to read or file-like object pointing to the
        beginning of such a file.
    debug : bool
        True to print out debugging information, False to supress
    read_fields : bool
        True to read all field during initalization, False (default) only
        reads metadata.

    Notes
    -----
    This Object is not stable enough to be considered a general MDV lib, nor is
    that our intention, but with careful use it shall provide full read/write capacity.

    """
    # ftm for use in the struct lib
    # mapper are used to convert vector to dics, they are of the following type (var_name,inicial pos, final pos, is_string)
    master_header_fmt = '>28i 8i i 5i 6f 3f 12f 512c 128c 128c i'
    master_header_mapper = [ ("record_len1",0,1 , False), ("struct_id",1,2 , False), ("revision_number",2,3 , False), ("time_gen",3,4 , False), ("user_time",4,5 , False), ("time_begin",5,6 , False), ("time_end",6,7 , False), ("time_centroid",7,8 , False), ("time_expire",8,9 , False), ("num_data_times",9,10 , False), ("index_number",10,11 , False), ("data_dimension",11,12 , False), ("data_collection_type",12,13 , False), ("user_data",13,14 , False), ("native_vlevel_type",14,15 , False), ("vlevel_type",15,16 , False), ("vlevel_included",16,17 , False), ("grid_orientation",17,18 , False), ("data_ordering",18,19 , False), ("nfields",19,20 , False), ("max_nx",20,21 , False), ("max_ny",21,22 , False), ("max_nz",22,23 , False), ("nchunks",23,24 , False), ("field_hdr_offset",24,25 , False), ("vlevel_hdr_offset",25,26 , False), ("chunk_hdr_offset",26,27 , False), ("field_grids_differ",27,28 , False), ("user_data_si328",28,36 , False), ("time_written",36,37 , False), ("unused_si325",37,42 , False), ("user_data_fl326",42,48 , False), ("sensor_lon",48,49 , False), ("sensor_lat",49,50 , False), ("sensor_alt",50,51 , False), ("unused_fl3212",51,63 , False), ("data_set_info",63,575 , True), ("data_set_name",575,703 , True), ("data_set_source",703,831 , True), ("record_len2",831,832 , False), ]
    
    field_header_fmt = '>17i 10i 9i 4i f f 8f 12f 4f 5f 64c 16c 16c 16c 16c i'
    field_header_mapper = [ ("record_len1",0,1,False), ("struct_id",1,2,False), ("field_code",2,3,False), ("user_time1",3,4,False), ("forecast_delta",4,5,False), ("user_time2",5,6,False), ("user_time3",6,7,False), ("forecast_time",7,8,False), ("user_time4",8,9,False), ("nx",9,10,False), ("ny",10,11,False), ("nz",11,12,False), ("proj_type",12,13,False), ("encoding_type",13,14,False), ("data_element_nbytes",14,15,False), ("field_data_offset",15,16,False), ("volume_size",16,17,False), ("user_data_si32",17,27,False), ("compression_type",27,28,False), ("transform_type",28,29,False), ("scaling_type",29,30,False), ("native_vlevel_type",30,31,False), ("vlevel_type",31,32,False), ("dz_constant",32,33,False), ("data_dimension",33,34,False), ("zoom_clipped",34,35,False), ("zoom_no_overlap",35,36,False), ("unused_si32",36,40,False), ("proj_origin_lat",40,41,False), ("proj_origin_lon",41,42,False), ("proj_param",42,50,False), ("vert_reference",50,51,False), ("grid_dx",51,52,False), ("grid_dy",52,53,False), ("grid_dz",53,54,False), ("grid_minx",54,55,False), ("grid_miny",55,56,False), ("grid_minz",56,57,False), ("scale",57,58,False), ("bias",58,59,False), ("bad_data_value",59,60,False), ("missing_data_value",60,61,False), ("proj_rotation",61,62,False), ("user_data_fl32",62,66,False), ("min_value",66,67,False), ("max_value",67,68,False), ("min_value_orig_vol",68,69,False), ("max_value_orig_vol",69,70,False), ("unused_fl32",70,71,False), ("field_name_long",71,135,True), ("field_name",135,151,True), ("units",151,167,True), ("transform",167,183,True), ("unused_char",183,199,True), ("record_len2",199,200,False),]
    
    vlevel_header_fmt = '>i i 122i 4i 122f 5f i'
    vlevel_header_mapper = [("record_len1",0,1,False), ("struct_id",1,2,False), ("type",2,124,False), ("unused_si32",124,128,False), ("level",128,250,False), ("unused_fl32",250,255,False), ("record_len2",255,256,False),] 
    
    chunk_header_fmt = '>5i 2i 480c i'
    chunk_header_mapper = [("record_len1",0,1,False), ("struct_id",1,2,False), ("chunk_id",2,3,False), ("chunk_data_offset",3,4,False), ("size",4,5,False), ("unused_si32",5,7,False), ("info",7,487,True), ("record_len2",487,488,False),]
    
    compression_info_fmt = '>I I I I 2I'
    compression_info_mapper = [("magic_cookie",0,1,False), ("nbytes_uncompressed",1,2,False), ("nbytes_compressed",2,3,False), ("nbytes_coded",3,4,False), ("spare",4,6,False),] 
    
    radar_info_fmt = '>12i 2i 22f 4f 40c 40c'
    radar_info_mapper = [("radar_id",0,1,False), ("radar_type",1,2,False), ("nfields",2,3,False), ("ngates",3,4,False), ("samples_per_beam",4,5,False), ("scan_type",5,6,False), ("scan_mode",6,7,False), ("nfields_current",7,8,False), ("field_flag",8,9,False), ("polarization",9,10,False), ("follow_mode",10,11,False), ("prf_mode",11,12,False), ("spare_ints",12,14,False), ("radar_constant",14,15,False), ("altitude_km",15,16,False), ("latitude_deg",16,17,False), ("longitude_deg",17,18,False), ("gate_spacing_km",18,19,False), ("start_range_km",19,20,False), ("horiz_beam_width_deg",20,21,False), ("vert_beam_width_deg",21,22,False), ("pulse_width_us",22,23,False), ("prf_hz",23,24,False), ("wavelength_cm",24,25,False), ("xmit_peak_pwr_watts",25,26,False), ("receiver_mds_dbm",26,27,False), ("receiver_gain_db",27,28,False), ("antenna_gain_db",28,29,False), ("system_gain_db",29,30,False), ("unambig_vel_mps",30,31,False), ("unambig_range_km",31,32,False), ("measXmitPowerDbmH_dbm",32,33,False), ("measXmitPowerDbmV_dbm",33,34,False), ("prt_s",34,35,False), ("prt2_s",35,36,False), ("spare_floats",36,40,False), ("radar_name",40,80,True), ("scan_type_name",80,120,True),]
    
    calib_fmt = '>16c 6i 51f 14f'
    calib_mapper = [("radar_name",0,16,True), ("year",16,17,False), ("month",17,18,False), ("day",18,19,False), ("hour",19,20,False), ("minute",20,21,False), ("second",21,22,False), ("wavelength_cm",22,23,False), ("beamwidth_h_deg",23,24,False), ("beamwidth_v_deg",24,25,False), ("antenna_gain_h_db",25,26,False), ("antenna_gain_v_db",26,27,False), ("pulse_width_us",27,28,False), ("xmit_power_h_dbm",28,29,False), ("xmit_power_v_dbm",29,30,False), ("twoway_waveguide_loss_h_db",30,31,False), ("twoway_waveguide_loss_v_db",31,32,False), ("twoway_radome_loss_h_db",32,33,False), ("twoway_radome_loss_v_db",33,34,False), ("filter_loss_db",34,35,False), ("radar_constant_h_db",35,36,False), ("radar_constant_v_db",36,37,False), ("noise_h_co_dbm",37,38,False), ("noise_h_cx_dbm",38,39,False), ("noise_v_co_dbm",39,40,False), ("noise_v_cx_dbm",40,41,False), ("rx_gain_h_co_dbm",41,42,False), ("rx_gain_h_cx_dbm",42,43,False), ("rx_gain_v_co_dbm",43,44,False), ("rx_gain_v_cx_dbm",44,45,False), ("zh1km_co_dbz",45,46,False), ("zh1km_cx_dbz",46,47,False), ("zv1km_co_dbz",47,48,False), ("zv1km_cx_dbz",48,49,False), ("sun_h_co_dbm",49,50,False), ("sun_h_cx_dbm",50,51,False), ("sun_v_co_dbm",51,52,False), ("sun_v_cx_dbm",52,53,False), ("noise_source_h_dbm",53,54,False), ("noise_source_v_dbm",54,55,False), ("power_meas_loss_h_db",55,56,False), ("power_meas_loss_v_db",56,57,False), ("coupler_fwd_loss_h_db",57,58,False), ("coupler_fwd_loss_v_db",58,59,False), ("zdr_bias_db",59,60,False), ("ldr_h_bias_db",60,61,False), ("ldr_v_bias_db",61,62,False), ("system_phidp_deg",62,63,False), ("test_pulse_h_dbm",63,64,False), ("test_pulse_v_dbm",64,65,False), ("rx_slope_h_co_db",65,66,False), ("rx_slope_h_cx_db",66,67,False), ("rx_slope_v_co_db",67,68,False), ("rx_slope_v_cx_db",68,69,False), ("I0_h_co_dbm",69,70,False), ("I0_h_cx_dbm",70,71,False), ("I0_v_co_dbm",71,72,False), ("I0_v_cx_dbm",72,73,False), ("spare",73,87,False), ]

    
    def __init__(self, filename, debug=False, read_fields=False):
        """ initalize MdvFile from filename (str). If filename=None create empty object"""
        if debug:
            print "Opening file for reading: ", filename
        if filename==None:
            self.fileptr = None #will creat empty struct, for filling and writing after
        elif hasattr(filename, 'read'):
            self.fileptr = filename
        else:
            self.fileptr = open(filename, 'rb')

        if debug:
            print "Getting master header"
        self.master_header = self._get_master_header()

        if debug:
            print "getting field headers"
        nfields = self.master_header['nfields']
        self.field_headers = self._get_field_headers(nfields)

        if debug:
            print "getting vlevel headers"
        self.vlevel_headers = self._get_vlevel_headers(nfields)

        if debug:
            print "getting chunk headers"
        nchunks = self.master_header['nchunks']
        self.chunk_headers = self._get_chunk_headers(nchunks)

        if debug:
            print "Getting Chunk Data"
        self.chunk_data = [None]*self.master_header['nchunks'] #will store raw chunk data, use for unkown chunk information
        self.radar_info, self.elevations, self.calib_info = self._get_chunks(debug)

        if self.master_header['nfields']>0:
            if self.field_headers[0]['proj_type'] == PROJ_LATLON:
                self.projection = 'latlon'
            elif self.field_headers[0]['proj_type'] == PROJ_LAMBERT_CONF:
                self.projection = 'lambert_conform'
            elif self.field_headers[0]['proj_type'] == PROJ_POLAR_STEREO:
                self.projection = 'polar_stereographic'
            elif self.field_headers[0]['proj_type'] == PROJ_FLAT:
                self.projection = 'flat'
            elif self.field_headers[0]['proj_type'] == PROJ_POLAR_RADAR:
                self.projection = 'ppi'
            elif self.field_headers[0]['proj_type'] == PROJ_OBLIQUE_STEREO:
                self.projection = 'oblique_stereographic'
            elif self.field_headers[0]['proj_type'] == PROJ_RHI_RADAR:
                self.projection = 'rhi'

#        if debug:
#            print "Calculating Radar coordinates"
#        az_deg, range_km, el_deg = self._calc_geometry()
#        self.az_deg = np.array(az_deg, dtype='float32')
#        self.range_km = np.array(range_km, dtype='float32')
#        self.el_deg = np.array(el_deg, dtype='float32')

        if debug:
            print "Making usable time objects"
        self.times = self._make_time_dict()

#        if debug:
#            print "Calculating cartesian coordinates"
#        self.carts = self._make_carts_dict()

#        if debug:
#            print "indexing fields"
#        self.fields = self._make_fields_list()
        
        self.fields_data = [None]*self.master_header["nfields"]
        self.compr_data=[[None]*head['nz'] for head in self.field_headers]
        if read_fields:
            if debug:
                print "Reading all fields"
            self.read_all_fields()
        return

    ##################
    # public methods #
    ##################
    def write(self, filename, debug=False):
        """ Write a MdvFile to filename (stg) """
        if debug:
            print "Opening file for writing:", filename
        if hasattr(filename, 'write'):
            self.fileptr = filename
        else:
            self.fileptr = open(filename, 'wb')
        file_start=self.fileptr.tell()
        
        #first write fields so one can calculate the offsets
        #put zero in headers

        headers_size=1024+(416+1024)*self.master_header["nfields"]+512*self.master_header["nchunks"]
        self.fileptr.write("\x00"*headers_size)
        
        if debug:
            print "Writing Fields Data"
        for ifield in range(self.master_header["nfields"]):
            self.write_a_field(ifield)
        
        #write chunks
        if debug:
            print "Writing Chunk Data"
        self._write_chunks(debug)
        #calculate offsets
        self._calc_file_offsets()
        self.fileptr.seek(file_start)
        #write headers
        if debug:
            print "Writing master header"
        self._write_master_header()

        if debug:
            print "Writing field headers"
        self._write_field_headers(self.master_header["nfields"])

        if debug:
            print "Writing vlevel headers"
        self._write_vlevel_headers(self.master_header["nfields"])

        if debug:
            print "Writing chunk headers"
        self._write_chunk_headers(self.master_header["nchunks"])
        # close file
        #XXX should I really do that? what if it's a file-like struct?
        if debug:
            print "Closing file"
        self.fileptr.close()
    
    
    def read_a_field(self, fnum, debug=False):
        """
        Read a field from the MDV file.

        Parameters
        ----------
        fnum : int
            Field number to read.
        debug : bool
            True to print debugging information, False to supress.

        Returns
        -------
        field_data : array
            Field data.  This data is also stored as a object attribute under
            the field name.

        See Also
        --------
        read_all_fields : Read all fields in the MDV file.

        """

        field_header = self.field_headers[fnum]
        # if the field has already been read, return it
        if self.fields_data[fnum]!=None:
            if debug:
                print "Getting data from the object."
            return self.fields_data[fnum]

        # field has not yet been read, populate the object and return
        if debug:
            print "No data found in object, populating"

        nz = field_header['nz']
        ny = field_header['ny']
        nx = field_header['nx']

        # read the header
        field_data = np.zeros([nz, ny, nx], dtype='float32')
        self.fileptr.seek(field_header['field_data_offset'])
        self._get_levels_info(nz)  # dict not used, but need to seek.

        for sw in xrange(nz):
            if debug:
                print "doing levels ", sw

            # get the compressed level data
            compr_info = self._get_compression_info()
            compr_data = self.fileptr.read(compr_info['nbytes_coded'])
            encoding_type = field_header['encoding_type']
            if encoding_type == ENCODING_INT8:
                fmt = '>%iB' % (nx * ny)
                np_form = '>B'
            elif encoding_type == ENCODING_INT16:
                fmt = '>%iH' % (nx * ny)
                np_form = '>H'
            elif encoding_type == ENCODING_FLOAT32:
                fmt = '>%if' % (nx * ny)
                np_form = '>f'
            else:
                raise ValueError('unknown encoding: ', encoding_type)

            # decompress the level data
            if compr_info['magic_cookie'] == 0xf7f7f7f7:
                cd_fobj = StringIO.StringIO(compr_data)
                gzip_file_handle = gzip.GzipFile(fileobj=cd_fobj)
                decompr_data = gzip_file_handle.read(struct.calcsize(fmt))
                gzip_file_handle.close()
            elif compr_info['magic_cookie'] == 0xf5f5f5f5:
                decompr_data = zlib.decompress(compr_data)
            else:
                raise NotImplementedError('unsupported compression mode')
                # With sample data it should be possible to write
                # decompressor for other modes, the compression magic
                # cookies for these modes are:
                # 0x2f2f2f2f : TA_NOT_COMPRESSED
                # 0xf8f8f8f8 : GZIP_NOT_COMPRSSED
                # 0xf3f3f3f3 : BZIP_COMPRESSED
                # 0xf4f4f4f4 : BZIP_NOT_COMPRESSED
                # 0xf6f6f6f6 : ZLIB_NOT_COMPRESSED
            self.read_data[fnum][sw]=decompr_data
            # read the decompressed data, reshape and mask
            sw_data = np.fromstring(decompr_data, np_form).astype('float32')
            sw_data.shape = (ny, nx)
            mask =  (sw_data == field_header['bad_data_value']) | (sw_data == field_header['missing_data_value'])
            np.putmask(sw_data, mask, [np.NaN])
            
            # scale and offset the data, store in field_data
            scale = field_header['scale']
            bias = field_header['bias']
            field_data[sw, :, :] = sw_data * scale + bias
            
        # store data as object attribute and return
        self.fields_data[fnum]= field_data
        return field_data

    def write_a_field(self, fnum, debug=False):
        """ write field number 'fnum' to mdv file """
        # the file pointer must be set at the correct location prior to call
        field_header = self.field_headers[fnum]
        if field_header['compression_type']!=3:
            import warnings
            warnings.warn("compression_type not implemented, converting to zlib")
            field_header['compression_type']=3

        field_data = self.fields_data[fnum]
        nz = field_header['nz']
        #save file posicion
        field_start=self.fileptr.tell()
        # write zeros to vlevel_offsets and vlevel_nbytes
        self.fileptr.write("\x00"*4*2*nz)
        field_size=0 
        vlevel_offsets=[0]*nz
        vlevel_nbytes=[0]*nz
        for sw in xrange(nz):
            vlevel_offsets[sw]=field_size
            scale = field_header['scale']
            bias = field_header['bias']
            sw_data = np.round((field_data[sw, :, :]-bias)/scale)
            if hasattr(sw_data, 'mask'):
                sw_data = np.where( sw_data.mask , field_header['bad_data_value'], sw_data)
            
            encoding_type = field_header['encoding_type']
            if encoding_type == ENCODING_INT8:
                np_form = '>B'
            elif encoding_type == ENCODING_INT16:
                np_form = '>H'
            elif encoding_type == ENCODING_FLOAT32:
                np_form = '>f'
            else:
                raise ValueError('unknown encoding: ', encoding_type)
            uncompr_data = np.array(sw_data,dtype=np_form).tostring()
            compr_data = zlib.compress(uncompr_data)
            if len(compr_data)>len(uncompr_data):
                magic = 0xf6f6f6f6
                compr_data = uncompr_data
            else:
                magic = 0xf5f5f5f5
            compr_info = {
                'magic_cookie': magic,
                'nbytes_uncompressed': len(uncompr_data),
                'nbytes_compressed': len(compr_data) + 24,
                'nbytes_coded': len(compr_data),
                'spare': [0,0],
            }

            self._write_compression_info(compr_info)
            self.fileptr.write(compr_data)
            field_size=field_size+len(compr_data) + 24
            vlevel_nbytes[sw]=len(compr_data) + 24
        #go back and rewrite vlevel_offsets and vlevel_nbytes
        field_end=self.fileptr.tell()
        self.fileptr.seek(field_start)
        fmt = '>%iI %iI' % (nz, nz)
        string = struct.pack(fmt, *(vlevel_offsets + vlevel_nbytes))
        self.fileptr.write(string)
        self.fileptr.seek(field_end)
        field_header["volume_size"]=field_size+2*4*nz


    def read_all_fields(self):
        """ Read all fields, storing data to field name attributes. """
        for i in xrange(self.master_header['nfields']):
            self.read_a_field(i)

    def close(self):
        """ Close the MDV file. """
        self.fileptr.close()

    ###################
    # private methods #
    ###################

    # get_ methods for reading headers

    def _get_master_header(self):
        """ Read the MDV master header, return a dict. """
        # the file pointer must be set at the correct location prior to call
        if self.fileptr==None:
            l = [0]*self.master_header_mapper[-1][2]
            l[0] = 1016
            l[1] = 14142
            l[2] = 1
            l[9] = 1
            l[16] = 1
            l[17] = 1
            l[831] = 1016
            for item in self.master_header_mapper:#if string convert to char
                if item[3]:
                    l[item[1]:item[2]] = [chr(a) for a in l[item[1]:item[2]]]
        else:
            l = struct.unpack(self.master_header_fmt, self.fileptr.read(struct.calcsize(self.master_header_fmt)))
        d = {}
        for item in self.master_header_mapper:
            if item[3]:
                d[item[0]] = "".join(l[item[1]:item[2]]).strip('\x00')
            else:
                if item[2]==item[1]+1:
                    d[item[0]] = l[item[1]]
                else:
                    d[item[0]] = l[item[1]:item[2]]
        return d

    def _write_master_header(self):
        """ Write the MDV master header. """
        # the file pointer must be set at the correct location prior to call
        d = self.master_header
        l=[0]*self.master_header_mapper[-1][2]
        for item in self.master_header_mapper:
            if item[3]:
                l[item[1]:item[2]] = (list(d[item[0]].encode("ASCII"))+['\x00']*(item[2]-item[1]))[0:item[2]-item[1]]#convert str to list of char and complet with zero
            else:
                if item[2]==item[1]+1:
                    l[item[1]] = d[item[0]]
                else:
                    l[item[1]:item[2]] = d[item[0]]
        string = struct.pack(self.master_header_fmt, *l)
        self.fileptr.write(string)

    def _get_field_headers(self, nfields):
        """ Read nfields field headers, return a list of dicts. """
        # the file pointer must be set at the correct location prior to call
        return [self._get_field_header() for i in range(nfields)]

    def _write_field_headers(self, nfields):
        """ Write nfields field headers. """
        # the file pointer must be set at the correct location prior to call
        for i in range(nfields):
            self._write_field_header(self.field_headers[i])

    def _get_field_header(self):
        """ Read a single field header, return a dict. """
        # the file pointer must be set at the correct location prior to call
        if self.fileptr==None:
            l = [0]*self.field_header_mapper[-1][2]
            l[0] = 408
            l[1] = 14143
            l[57] = 1 #scale
            l[199] = 408
            for item in self.field_header_mapper:#if string convert to char
                if item[3]:
                    l[item[1]:item[2]] = [chr(a) for a in l[item[1]:item[2]]]
        else:
            l = struct.unpack(self.field_header_fmt, self.fileptr.read(struct.calcsize(self.field_header_fmt)))
        d = {}
        for item in self.field_header_mapper:
            if item[3]:
                d[item[0]] = "".join(l[item[1]:item[2]]).strip('\x00')
            else:
                if item[2]==item[1]+1:
                    d[item[0]] = l[item[1]]
                else:
                    d[item[0]] = l[item[1]:item[2]]
        return d

    def _write_field_header(self,d):
        """ Write the a single field header. """
        # the file pointer must be set at the correct location prior to call
        l=[0]*self.field_header_mapper[-1][2]
        for item in self.field_header_mapper:
            if item[3]:
                l[item[1]:item[2]] = (list(d[item[0]].encode("ASCII"))+['\x00']*(item[2]-item[1]))[0:item[2]-item[1]]#convert str to list of char and complet with zero
            else:
                if item[2]==item[1]+1:
                    l[item[1]] = d[item[0]]
                else:
                    l[item[1]:item[2]] = d[item[0]]
        string = struct.pack(self.field_header_fmt, *l)
        self.fileptr.write(string)

    def _get_vlevel_headers(self, nfields):
        """ Read nfields vlevel headers, return a list of dicts. """
        # the file pointer must be set at the correct location prior to call
        return [self._get_vlevel_header() for i in range(nfields)]

    def _write_vlevel_headers(self, nfields):
        """ Write nfields vlevel headers"""
        # the file pointer must be set at the correct location prior to call
        for i in range(nfields):
            self._write_vlevel_header(self.vlevel_headers[i])

    def _get_vlevel_header(self):
        """ Read a single vlevel header, return a dict. """
        # the file pointer must be set at the correct location prior to call
        if self.fileptr==None:
            l = [0]*self.vlevel_header_mapper[-1][2]
            l[0] = 1016
            l[1] = 14144
            l[255] = 1016
            for item in self.vlevel_header_mapper:#if string convert to char
                if item[3]:
                    l[item[1]:item[2]] = [chr(a) for a in l[item[1]:item[2]]]
        else:
            l = struct.unpack(self.vlevel_header_fmt, self.fileptr.read(struct.calcsize(self.vlevel_header_fmt)))
        d = {}
        for item in self.vlevel_header_mapper:
            if item[3]:
                d[item[0]] = "".join(l[item[1]:item[2]]).strip('\x00')
            else:
                if item[2]==item[1]+1:
                    d[item[0]] = l[item[1]]
                else:
                    d[item[0]] = l[item[1]:item[2]]
        return d

    def _write_vlevel_header(self,d):
        """  Write the a single vfield header. """
        # the file pointer must be set at the correct location prior to call
        l=[0]*self.vlevel_header_mapper[-1][2]
        for item in self.vlevel_header_mapper:
            if item[3]:
                l[item[1]:item[2]] = (list(d[item[0]].encode("ASCII"))+['\x00']*(item[2]-item[1]))[0:item[2]-item[1]]#convert str to list of char and complet with zero
            else:
                if item[2]==item[1]+1:
                    l[item[1]] = d[item[0]]
                else:
                    l[item[1]:item[2]] = d[item[0]]
        string = struct.pack(self.vlevel_header_fmt, *l)
        self.fileptr.write(string)

    def _get_chunk_headers(self, nchunks):
        """ Get nchunk chunk headers, return a list of dicts. """
        # the file pointer must be set at the correct location prior to call
        return [self._get_chunk_header() for i in range(nchunks)]

    def _write_chunk_headers(self, nchunks):
        """ Write nchunk chunk headers. """
        # the file pointer must be set at the correct location prior to call
        for i in range(nchunks):
            self._write_chunk_header(self.chunk_headers[i])

    def _get_chunk_header(self):
        """ Get a single chunk header, return a dict. """
        # the file pointer must be set at the correct location prior to call
        if self.fileptr==None:
            l = [0]*self.chunk_header_mapper[-1][2]
            l[0] = 504
            l[1] = 14145
            l[487] = 504
            for item in self.chunk_header_mapper:#if string convert to char
                if item[3]:
                    l[item[1]:item[2]] = [chr(a) for a in l[item[1]:item[2]]]
        else:
            l = struct.unpack(self.chunk_header_fmt, self.fileptr.read(struct.calcsize(self.chunk_header_fmt)))
        d = {}
        for item in self.chunk_header_mapper:
            if item[3]:
                d[item[0]] = "".join(l[item[1]:item[2]]).strip('\x00')
            else:
                if item[2]==item[1]+1:
                    d[item[0]] = l[item[1]]
                else:
                    d[item[0]] = l[item[1]:item[2]]
        return d

    def _write_chunk_header(self,d):
        """  Write the a single chunk header. """
        # the file pointer must be set at the correct location prior to call
        l=[0]*self.chunk_header_mapper[-1][2]
        for item in self.chunk_header_mapper:
            if item[3]:
                l[item[1]:item[2]] = (list(d[item[0]].encode("ASCII"))+['\x00']*(item[2]-item[1]))[0:item[2]-item[1]]#convert str to list of char and complet with zero
            else:
                if item[2]==item[1]+1:
                    l[item[1]] = d[item[0]]
                else:
                    l[item[1]:item[2]] = d[item[0]]
        string = struct.pack(self.chunk_header_fmt, *l)
        self.fileptr.write(string)

    def _get_chunks(self, debug=False):
        """ Get data in chunks, return radar_info, elevations, calib_info. """
        # the file pointer must be set at the correct location prior to call
        radar_info, elevations, calib_info = None, [], None
        for cnum,curr_chunk_header in enumerate(self.chunk_headers):

            chunk_id = curr_chunk_header['chunk_id']
            self.fileptr.seek(curr_chunk_header['chunk_data_offset'])

            if chunk_id == CHUNK_DSRADAR_PARAMS:
                if debug:
                    print 'Getting radar info'
                radar_info = self._get_radar_info()

            elif chunk_id == CHUNK_DSRADAR_CALIB:
                if debug:
                    print 'getting elevations'
                elevations = self._get_elevs(curr_chunk_header['size'])

            elif chunk_id == CHUNK_DSRADAR_ELEVATIONS:
                if debug:
                    print 'getting cal'
                calib_info = self._get_calib()

            else:
                if debug:
                    print 'getting unknown chunk %i'%chunk_id
                self.chunk_data[cnum] = self._get_unknown_chunk(cnum)

        return radar_info, elevations, calib_info

    def _write_chunks(self, debug=False):
        """ write chunks data """
        # the file pointer must be set at the correct location prior to call
        for curr_chunk_header in self.chunk_headers:
            chunk_id = curr_chunk_header['chunk_id']
            
            if chunk_id == CHUNK_DSRADAR_PARAMS:
                if debug:
                    print 'writing radar info'
                radar_info = self._write_radar_info(self.radar_info)
            
            elif chunk_id == CHUNK_DSRADAR_ELEVATIONS:
                if debug:
                    print 'writing elevations'
                elevations = self._write_elevs(self.elevations)
            
            elif chunk_id == CHUNK_DSRADAR_CALIB:
                if debug:
                    print 'writing cal'
                calib_info = self._write_calib(self.calib_info)
            
            else:
                if debug:
                    print 'writing unknown chunk %i'%chunk_id
                self._write_unknown_chunk(self,self.chunk_data[cnum])
            
    def _get_radar_info(self):
        """ Get the radar information, return dict. """
        # the file pointer must be set at the correct location prior to call
        if self.fileptr==None:
            l = [0]*self.radar_info_mapper[-1][2]
            for item in self.radar_info_mapper:#if string convert to char
                if item[3]:
                    l[item[1]:item[2]] = [chr(a) for a in l[item[1]:item[2]]]
        else:
            l = struct.unpack(self.radar_info_fmt, self.fileptr.read(struct.calcsize(self.radar_info_fmt)))
        d = {}
        for item in self.radar_info_mapper:
            if item[3]:
                d[item[0]] = "".join(l[item[1]:item[2]]).strip('\x00')
            else:
                if item[2]==item[1]+1:
                    d[item[0]] = l[item[1]]
                else:
                    d[item[0]] = l[item[1]:item[2]]
        return d

    def _write_radar_info(self,d):
        """  Write radar information. """
        # the file pointer must be set at the correct location prior to call
        l=[0]*self.radar_info_mapper[-1][2]
        for item in self.radar_info_mapper:
            if item[3]:
                l[item[1]:item[2]] = (list(d[item[0]].encode("ASCII"))+['\x00']*(item[2]-item[1]))[0:item[2]-item[1]]#convert str to list of char and complet with zero
            else:
                if item[2]==item[1]+1:
                    l[item[1]] = d[item[0]]
                else:
                    l[item[1]:item[2]] = d[item[0]]
        string = struct.pack(self.radar_info_fmt, *l)
        self.fileptr.write(string)

    def _get_elevs(self, nbytes):
        """ Return an array of elevation read from current file position. """
        # the file pointer must be set at the correct location prior to call
        SIZE_FLOAT = 4.0
        nelevations = np.floor(nbytes / SIZE_FLOAT)
        fmt = '>%df' % (nelevations)
        l = struct.unpack(fmt, self.fileptr.read(struct.calcsize(fmt)))
        return np.array(l)

    def _write_elevs(self, l):
        """ Write an array of elevation. """
        # the file pointer must be set at the correct location prior to call
        fmt = '>%df' % (len(l))
        string = struct.pack(fmt, *l)
        self.fileptr.write(string)

    def _get_calib(self):
        """ Get the calibration information, return a dict. """
        # the file pointer must be set at the correct location prior to call
        if self.fileptr==None:
            l = [0]*self.calib_mapper[-1][2]
            for item in self.calib_mapper:#if string convert to char
                if item[3]:
                    l[item[1]:item[2]] = [chr(a) for a in l[item[1]:item[2]]]
        else:
            l = struct.unpack(self.calib_fmt, self.fileptr.read(struct.calcsize(self.calib_fmt)))
        d = {}
        for item in self.calib_mapper:
            if item[3]:
                d[item[0]] = "".join(l[item[1]:item[2]]).strip('\x00')
            else:
                if item[2]==item[1]+1:
                    d[item[0]] = l[item[1]]
                else:
                    d[item[0]] = l[item[1]:item[2]]
        return d

    def _write_calib(self,d):
        """  Write calibration information. """
        # the file pointer must be set at the correct location prior to call
        l=[0]*self.calib_mapper[-1][2]
        for item in self.calib_mapper:
            if item[3]:
                l[item[1]:item[2]] = (list(d[item[0]].encode("ASCII"))+['\x00']*(item[2]-item[1]))[0:item[2]-item[1]]#convert str to list of char and complet with zero
            else:
                if item[2]==item[1]+1:
                    l[item[1]] = d[item[0]]
                else:
                    l[item[1]:item[2]] = d[item[0]]
        string = struct.pack(self.calib_fmt, *l)
        self.fileptr.write(string)

    def _get_compression_info(self):
        """ Get compression infomation, return a dict. """
        # the file pointer must be set at the correct location prior to call
        if self.fileptr==None:
            l = [0]*self.compression_info_mapper[-1][2]
        else:
            l = struct.unpack(self.compression_info_fmt, self.fileptr.read(struct.calcsize(self.compression_info_fmt)))
        d = {}
        for item in self.compression_info_mapper:
            if item[3]:
                d[item[0]] = "".join(l[item[1]:item[2]]).strip('\x00')
            else:
                if item[2]==item[1]+1:
                    d[item[0]] = l[item[1]]
                else:
                    d[item[0]] = l[item[1]:item[2]]
        return d

    def _write_compression_info(self,d):
        """ Write compression infomation"""
        # the file pointer must be set at the correct location prior to call
        l=[0]*self.compression_info_mapper[-1][2]
        for item in self.compression_info_mapper:
            if item[3]:
                l[item[1]:item[2]] = (list(d[item[0]].encode("ASCII"))+['\x00']*(item[2]-item[1]))[0:item[2]-item[1]]#convert str to list of char and complet with zero
            else:
                if item[2]==item[1]+1:
                    l[item[1]] = d[item[0]]
                else:
                    l[item[1]:item[2]] = d[item[0]]
        string = struct.pack(self.compression_info_fmt, *l)
        self.fileptr.write(string)

    def _get_unknown_chunk(self,cnum):
        """ Get raw data from chunk """
        # the file pointer must be set at the correct location prior to call
        size = self.chunk_headers[cnum]['size']
        return self.fileptr.read(size)

    def _write_unknown_chunk(self,data):
        """ Write raw data from chunk """
        # the file pointer must be set at the correct location prior to call
        self.fileptr.write(data)

    def _get_levels_info(self, nlevels):
        """ Get nlevel information, return a dict. """
        # the file pointer must be set at the correct location prior to call
        fmt = '>%iI %iI' % (nlevels, nlevels)
        if self.fileptr:
            l = struct.unpack(fmt, self.fileptr.read(struct.calcsize(fmt)))
        else:
            l = [0]*2*nlevels
        d = {}
        d['vlevel_offsets'] = l[:nlevels]
        d['vlevel_nbytes'] = l[nlevels:2*nlevels]
        return d

    def _write_levels_info(self, nlevels,d):
        """ write levels information, return a dict. """
        # the file pointer must be set at the correct location prior to call
        fmt = '%iI %iI' % (nlevels, nlevels)
        l= d['vlevel_offsets']+d['vlevel_nbytes']
        string = struct.pack(fmt, *l)
        self.fileptr.write(string)

    def _calc_file_offsets(self):
        self.master_header["field_hdr_offset"]=1024
        self.master_header["vlevel_hdr_offset"]=1024+416*self.master_header["nfields"]
        self.master_header["chunk_hdr_offset"]=1024+(416+1024)*self.master_header["nfields"]
        
        file_pos=self.master_header["chunk_hdr_offset"]+512*self.master_header["nchunks"]
        for i in range(self.master_header["nfields"]):
            self.field_headers[i]["field_data_offset"]=file_pos
            file_pos=file_pos+self.field_headers[i]["volume_size"]
        
        for i in range(self.master_header["nchunks"]):
            self.chunk_headers[i]["chunk_data_offset"]=file_pos
            file_pos=file_pos+self.chunk_headers[i]["size"]

    def _make_time_dict(self):
        """ Return a time dictionary. """
        t_base = datetime.datetime(1970, 1, 1, 00, 00)
        tb = datetime.timedelta(seconds=self.master_header['time_begin'])
        te = datetime.timedelta(seconds=self.master_header['time_end'])
        tc = datetime.timedelta(seconds=self.master_header['time_centroid'])
        return {'time_begin': t_base + tb, 'time_end': t_base + te,
                'time_centroid': t_base + tc}

    def _time_dict_into_header(self):
        """ Complete time information in master_header from the time dict """
        t_base = datetime.datetime(1970, 1, 1, 00, 00)
        self.master_header['time_begin'] = (self.times['time_begin']- t_base).total_seconds()
        self.master_header['time_end'] = (self.times['time_end']- t_base).total_seconds()
        self.master_header['time_centroid'] = (self.times['time_centroid']- t_base).total_seconds()



    
    # misc. methods
    #XXX move some where else, there are not general mdv operations

    def _calc_geometry(self):
        """ Calculate geometry, return az_deg, range_km, el_deg. """
        nsweeps = self.master_header['max_nz']
        nrays = self.master_header['max_ny']
        ngates = self.master_header['max_nx']
        grid_minx = self.field_headers[0]['grid_minx']
        grid_miny = self.field_headers[0]['grid_miny']
        grid_dx = self.field_headers[0]['grid_dx']
        grid_dy = self.field_headers[0]['grid_dy']

        range_km = grid_minx + np.arange(ngates) * grid_dx

        if self.field_headers[0]['proj_type'] == PROJ_RHI_RADAR:
            self.scan_type = 'rhi'
            el_deg = grid_miny + np.arange(nrays) * grid_dy
            az_deg = self.vlevel_headers[0]['level'][0:nsweeps]

        if self.field_headers[0]['proj_type'] == PROJ_POLAR_RADAR:
            self.scan_type = 'ppi'
            az_deg = grid_miny + np.arange(nrays) * grid_dy
            el_deg = self.vlevel_headers[0]['level'][0:nsweeps]

        return az_deg, range_km, el_deg



    def _make_carts_dict(self):
        """ Return a carts dictionary, distances in meters. """

        # simple calculation involving 4/3 earth radius
        nsweeps = self.master_header['max_nz']
        nrays = self.master_header['max_ny']
        ngates = self.master_header['max_nx']
        xx = np.empty([nsweeps, nrays, ngates], dtype=np.float32)
        yy = np.empty([nsweeps, nrays, ngates], dtype=np.float32)
        zz = np.empty([nsweeps, nrays, ngates], dtype=np.float32)

        if self.scan_type == 'rhi':
            rg, ele = np.meshgrid(self.range_km, self.el_deg)
            rg = np.array(rg, dtype=np.float64)
            ele = np.array(ele, dtype=np.float64)
            for aznum in xrange(nsweeps):
                azg = np.ones(rg.shape, dtype=np.float64) * self.az_deg[aznum]
                x, y, z = radar_coords_to_cart(rg, azg, ele)
                zz[aznum, :, :] = z
                xx[aznum, :, :] = x
                yy[aznum, :, :] = y

        elif self.scan_type == 'ppi':
            rg, azg = np.meshgrid(self.range_km, self.az_deg)
            rg = np.array(rg, dtype=np.float64)
            azg = np.array(azg, dtype=np.float64)
            for elnum in xrange(nsweeps):
                ele = np.ones(rg.shape, dtype=np.float64) * self.el_deg[elnum]
                x, y, z = radar_coords_to_cart(rg, azg, ele)
                zz[elnum, :, :] = z
                xx[elnum, :, :] = x
                yy[elnum, :, :] = y

        return {'x': xx, 'y': yy, 'z': zz}

    def _make_fields_list(self):
        """ Return a list of fields. """
        fh = self.field_headers
        return [fh[i]['field_name'] for i in range(len(fh))]
