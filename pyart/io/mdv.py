"""
pyart.io.mdv
============

Utilities for reading of MDV files.

Code is adapted from Nitin Bharadwaj's Matlab code

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    MdvFile

.. autosummary::
    :toctree: generated/

    read_mdv

"""

import struct
import gzip
import StringIO
import datetime

import numpy as np
from netCDF4 import date2num

from .radar import Radar
from .common import COMMON2STANDARD, get_metadata, make_time_unit_str
from .common import radar_coords_to_cart


def read_mdv(filename):
    """
    Read a MDV file.

    Parameters
    ----------
    filename : str
        Name of MDV file to read data from.

    Returns
    -------
    radar : Radar
        Radar object containing data from MDV file.

    Notes
    -----
    Currently this function can only read polar MDV files which are gzipped.
    Support for cartesian and non-gzipped file are planned.

    """
    mdvfile = MdvFile(filename)

    # value attributes
    naz = len(mdvfile.az_deg)
    nele = len(mdvfile.el_deg)
    scan_type = mdvfile.scan_type

    if scan_type not in ['ppi', 'rhi']:
        raise NotImplementedError('No support for scan_type %s.' % scan_type)

    # time
    time = get_metadata('time')
    units = make_time_unit_str(mdvfile.times['time_begin'])
    time['units'] = units
    time_start = date2num(mdvfile.times['time_begin'], units)
    time_end = date2num(mdvfile.times['time_end'], units)
    time['data'] = np.linspace(time_start, time_end, naz * nele)

    # range
    _range = get_metadata('range')
    _range['data'] = np.array(mdvfile.range_km * 1000.0, dtype='float32')
    _range['meters_to_center_of_first_gate'] = _range['data'][0]
    gate_0 = _range['data'][0]
    gate_1 = _range['data'][1]
    _range['meters_between_gates'] = (gate_1 - gate_0)

    # fields
    fields = {}
    # Transfer only the fields that we have valid names for and are present.
    for field in set(mdvfile.fields) & set(COMMON2STANDARD.keys()):

        # grab data from MDV object, mask and reshape
        data = mdvfile.read_a_field(mdvfile.fields.index(field))
        data[np.where(np.isnan(data))] = -9999.0
        data[np.where(data == 131072)] = -9999.0
        data = np.ma.masked_equal(data, -9999.0)
        data.shape = (data.shape[0] * data.shape[1], data.shape[2])

        # create and store the field dictionary
        fielddict = get_metadata(field)
        fielddict['data'] = data
        fielddict['_FillValue'] = -9999.0
        fields[COMMON2STANDARD[field]] = fielddict

    # metadata
    metadata = {}
    for meta_key, mdv_key in MDV_METADATA_MAP.iteritems():
        metadata[meta_key] = mdvfile.master_header[mdv_key]

    # additional required CF/Radial metadata set to blank strings
    metadata['title'] = ''
    metadata['institution'] = ''
    metadata['references'] = ''
    metadata['source'] = ''
    metadata['history'] = ''
    metadata['comment'] = ''

    # latitude
    latitude = get_metadata('latitude')
    latitude['data'] = np.array([mdvfile.radar_info['latitude_deg']],
                                dtype='float64')
    # longitude
    longitude = get_metadata('longitude')
    longitude['data'] = np.array([mdvfile.radar_info['longitude_deg']],
                                 dtype='float64')
    # altitude
    altitude = get_metadata('altitude')
    altitude['data'] = np.array([mdvfile.radar_info['altitude_km'] * 1000.0],
                                dtype='float64')

    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
    # sweep_end_ray_index
    sweep_number = get_metadata('sweep_number')
    sweep_mode = get_metadata('sweep_mode')
    fixed_angle = get_metadata('fixed_angle')
    sweep_start_ray_index = get_metadata('sweep_start_ray_index')
    sweep_end_ray_index = get_metadata('sweep_end_ray_index')
    len_time = len(time['data'])

    if mdvfile.scan_type == 'ppi':
        nsweeps = nele
        sweep_number['data'] = np.arange(nsweeps, dtype='int32')
        sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])
        fixed_angle['data'] = np.array(mdvfile.el_deg, dtype='float32')
        sweep_start_ray_index['data'] = np.arange(0, len_time, naz,
                                                  dtype='int32')
        sweep_end_ray_index['data'] = np.arange(naz-1, len_time, naz,
                                                dtype='int32')

    elif mdvfile.scan_type == 'rhi':
        nsweeps = naz
        sweep_number['data'] = np.arange(nsweeps, dtype='int32')
        sweep_mode['data'] = np.array(nsweeps * ['rhi'])
        fixed_angle['data'] = np.array(mdvfile.az_deg, dtype='float32')
        sweep_start_ray_index['data'] = np.arange(0, len_time, nele,
                                                  dtype='int32')
        sweep_end_ray_index['data'] = np.arange(nele - 1, len_time, nele,
                                                dtype='int32')

    # azimuth, elevation
    azimuth = get_metadata('azimuth')
    elevation = get_metadata('elevation')

    if scan_type == 'ppi':
        azimuth['data'] = np.tile(mdvfile.az_deg, nele)
        elevation['data'] = np.array(mdvfile.el_deg).repeat(naz)

    elif scan_type == 'rhi':
        azimuth['data'] = np.array(mdvfile.az_deg).repeat(nele)
        elevation['data'] = np.tile(mdvfile.el_deg, naz)

    # instrument parameters
    # we will set 4 parameters in the instrument_parameters dict
    # prt, prt_mode, unambiguous_range, and nyquist_velocity

    # TODO prt mode: Need to fix this.. assumes dual if two prts
    if mdvfile.radar_info['prt2_s'] == 0.0:
        prt_mode_str = 'fixed'
    else:
        prt_mode_str = 'dual'

    prt_mode = get_metadata('prt_mode')
    prt = get_metadata('prt')
    unambiguous_range = get_metadata('unambiguous_range')
    nyquist_velocity = get_metadata('nyquist_velocity')

    prt_mode['data'] = np.array([prt_mode_str] * nsweeps)
    prt['data'] = np.array([mdvfile.radar_info['prt_s']] * nele * naz,
                           dtype='float32')

    urange_m = mdvfile.radar_info['unambig_range_km'] * 1000.0
    unambiguous_range['data'] = np.array([urange_m] * naz * nele,
                                         dtype='float32')

    uvel_mps = mdvfile.radar_info['unambig_vel_mps']
    nyquist_velocity['data'] = np.array([uvel_mps] * naz * nele,
                                        dtype='float32')

    instrument_parameters = {'prt_mode': prt_mode, 'prt': prt,
                             'unambiguous_range': unambiguous_range,
                             'nyquist_velocity': nyquist_velocity}

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
    filename : str
        Name of MDV file.
    debug : bool
        True to print out debugging information, False to supress
    read_fields : bool
        True to read all field during initalization, False (default) only
        reads metadata.

    """

    def __init__(self, filename, debug=False, read_fields=False):
        """ initalize MdvFile from filename (str). """

        if debug:
            print "Opening ", filename
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
        radar_info, elevations, calib_info = self._get_chunks()
        if radar_info is not None:
            self.radar_info = radar_info
        if elevations is not None:
            self.elevations = elevations
        if calib_info is not None:
            self.calib_info = calib_info

        if debug:
            print "Calculating Radar coordinates"
        az_deg, range_km, el_deg = self._calc_geometry()
        self.az_deg = np.array(az_deg, dtype='float32')
        self.range_km = np.array(range_km, dtype='float32')
        self.el_deg = np.array(el_deg, dtype='float32')

        if debug:
            print "Making usable time objects"
        self.times = self._make_time_dict()

        if debug:
            print "Calculating cartesian coordinates"
        self.carts = self._make_carts_dict()

        if debug:
            print "indexing fields"
        self.fields = self._make_fields_list()

        if read_fields:
            if debug:
                print "Reading all fields"
            self.read_all_fields()
        return

    ##################
    # public methods #
    ##################

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
        if hasattr(self, field_header['field_name']):
            if debug:
                print "Getting data from the object."
            return getattr(self, field_header['field_name'])

        # field has not yet been read, populate the object and return
        if debug:
            print "No data found in object, populating"

        nsweeps = self.master_header['nsweeps']
        nrays = self.master_header['nrays']
        ngates = self.master_header['ngates']

        # read the header
        field_data = np.zeros([nsweeps, nrays, ngates], dtype='float32')
        self.fileptr.seek(field_header['field_data_offset'])
        self._get_sweep_info(nsweeps)  # dict not used, but need to seek.

        for sw in xrange(nsweeps):
            if debug:
                print "doing sweep ", sw

            # get the compressed sweep data
            compr_info = self._get_compression_info()
            compr_data = self.fileptr.read(compr_info['nbytes_coded'])
            cd_fobj = StringIO.StringIO(compr_data)

            # decompress the sweep data
            gzip_file_handle = gzip.GzipFile(fileobj=cd_fobj)
            encoding_type = field_header['encoding_type']
            if encoding_type == ENCODING_INT8:
                fmt = '>%iB' % (ngates * nrays)
                np_form = '>B'
            elif encoding_type == ENCODING_INT16:
                fmt = '>%iH' % (ngates * nrays)
                np_form = '>H'
            elif encoding_type == ENCODING_FLOAT32:
                fmt = '>%if' % (ngates * nrays)
                np_form = '>f'
            else:
                raise ValueError('unknown encoding: ', encoding_type)
            decompr_data = gzip_file_handle.read(struct.calcsize(fmt))
            gzip_file_handle.close()

            # read the decompressed data, reshape and mask
            sw_data = np.fromstring(decompr_data, np_form).astype('float32')
            sw_data.shape = (nrays, ngates)
            mask = sw_data == field_header['bad_data_value']
            np.putmask(sw_data, mask, [np.NaN])

            # scale and offset the data, store in field_data
            scale = field_header['scale']
            bias = field_header['bias']
            field_data[sw, :, :] = sw_data * scale + bias

        # store data as object attribute and return
        setattr(self, field_header['field_name'], field_data)
        return field_data

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
        fmt = '>28i 8i i 5i 6f 3f 12f 512c 128c 128c i'
        l = struct.unpack(fmt, self.fileptr.read(struct.calcsize(fmt)))
        d = {}
        d["record_len1"] = l[0]             # 28i: 1
        d["struct_id"] = l[1]               # 28i: 2
        d["revision_number"] = l[2]         # 28i: 3
        d["time_gen"] = l[3]                # 28i: 4
        d["user_time"] = l[4]               # 28i: 5
        d["time_begin"] = l[5]              # 28i: 6
        d["time_end"] = l[6]                # 28i: 7
        d["time_centroid"] = l[7]           # 28i: 8
        d["time_expire"] = l[8]             # 28i: 9
        d["num_data_times"] = l[9]          # 28i: 10
        d["index_number"] = l[10]           # 28i: 11
        d["data_dimension"] = l[11]         # 28i: 12
        d["data_collection_type"] = l[12]   # 28i: 13
        d["user_data"] = l[13]              # 28i: 14
        d["native_vlevel_type"] = l[14]     # 28i: 15
        d["vlevel_type"] = l[15]            # 28i: 16
        d["vlevel_included"] = l[16]        # 28i: 17
        d["grid_orientation"] = l[17]       # 28i: 18
        d["data_ordering"] = l[18]          # 28i: 19
        d["nfields"] = l[19]                # 28i: 20
        d["ngates"] = l[20]                 # 28i: 21
        d["nrays"] = l[21]                  # 28i: 22
        d["nsweeps"] = l[22]                # 28i: 23
        d["nchunks"] = l[23]                # 28i: 24
        d["field_hdr_offset"] = l[24]       # 28i: 25
        d["vlevel_hdr_offset"] = l[25]      # 28i: 26
        d["chunk_hdr_offset"] = l[26]       # 28i: 27
        d["field_grids_differ"] = l[27]     # 28i: 28
        d["user_data_si328"] = l[28:36]     # 8i
        d["time_written"] = l[36]           # i
        d["unused_si325"] = l[37:42]        # 5i
        d["user_data_fl326"] = l[42:48]     # 6f
        d["sensor_lon"] = l[48]             # 3f : 1
        d["sensor_lat"] = l[49]             # 3f : 2
        d["sensor_alt"] = l[50]             # 3f : 3
        d["unused_fl3212"] = l[51:63]       # 12f
        d["data_set_info"] = ''.join(l[63:575]).strip('\x00')       # 512c
        d["data_set_name"] = ''.join(l[575:703]).strip('\x00')      # 128c
        d["data_set_source"] = ''.join(l[703:831]).strip('\x00')    # 128c
        d["record_len2"] = l[831]           # i
        return d

    def _get_field_headers(self, nfields):
        """ Read nfields field headers, return a list of dicts. """
        return [self._get_field_header() for i in range(nfields)]

    def _get_field_header(self):
        """ Read a single field header, return a dict. """

        fmt = '>17i 10i 9i 4i f f 8f 12f 4f 5f 64c 16c 16c 16c 16c i'
        l = struct.unpack(fmt, self.fileptr.read(struct.calcsize(fmt)))
        d = {}
        d["record_len1"] = l[0]             # 17i: 1
        d["struct_id"] = l[1]               # 17i: 2
        d["field_code"] = l[2]              # 17i: 3
        d["user_time1"] = l[3]              # 17i: 4
        d["forecast_delta"] = l[4]          # 17i: 5
        d["user_time2"] = l[5]              # 17i: 6
        d["user_time3"] = l[6]              # 17i: 7
        d["forecast_time"] = l[7]           # 17i: 8
        d["user_time4"] = l[8]              # 17i: 9
        d["ngates"] = l[9]                  # 17i: 10
        d["nrays"] = l[10]                  # 17i: 11
        d["nsweeps"] = l[11]                # 17i: 12
        d["proj_type"] = l[12]              # 17i: 13
        d["encoding_type"] = l[13]          # 17i: 14
        d["data_element_nbytes"] = l[14]    # 17i: 15
        d["field_data_offset"] = l[15]      # 17i: 16
        d["volume_size"] = l[16]            # 17i: 17
        d["user_data_si32"] = l[17:27]      # 10i
        d["compression_type"] = l[27]       # 9i: 1
        d["transform_type"] = l[28]         # 9i: 2
        d["scaling_type"] = l[29]           # 9i: 3
        d["native_vlevel_type"] = l[30]     # 9i: 4
        d["vlevel_type"] = l[31]            # 9i: 5
        d["dz_constant"] = l[32]            # 9i: 6
        d["data_dimension"] = l[33]         # 9i: 7
        d["zoom_clipped"] = l[34]           # 9i: 8
        d["zoom_no_overlap"] = l[35]        # 9i: 9
        d["unused_si32"] = l[36:40]         # 4i
        d["proj_origin_lat"] = l[40]        # f
        d["proj_origin_lon"] = l[41]        # f
        d["proj_param"] = l[42:50]          # 8f
        d["vert_reference"] = l[50]         # 12f: 1
        d["grid_dx"] = l[51]                # 12f: 2
        d["grid_dy"] = l[52]                # 12f: 3
        d["grid_dz"] = l[53]                # 12f: 4
        d["grid_minx"] = l[54]              # 12f: 5
        d["grid_miny"] = l[55]              # 12f: 6
        d["grid_minz"] = l[56]              # 12f: 7
        d["scale"] = l[57]                  # 12f: 8
        d["bias"] = l[58]                   # 12f: 9
        d["bad_data_value"] = l[59]         # 12f: 10
        d["missing_data_value"] = l[60]     # 12f: 11
        d["proj_rotation"] = l[61]          # 12f: 12
        d["user_data_fl32"] = l[62:66]      # 4f
        d["min_value"] = l[66]              # 5f: 1
        d["max_value"] = l[67]              # 5f: 2
        d["min_value_orig_vol"] = l[68]     # 5f: 3
        d["max_value_orig_vol"] = l[69]     # 5f: 4
        d["unused_fl32"] = l[70]            # 5f: 5
        d["field_name_long"] = ''.join(l[71:135]).strip('\x00')     # 64c
        d["field_name"] = ''.join(l[135:151]).strip('\x00')         # 16c
        d["units"] = ''.join(l[151:167]).strip('\x00')              # 16c
        d["transform"] = ''.join(l[167:183]).strip('\x00')          # 16c
        d["unused_char"] = ''.join(l[183:199]).strip('\x00')        # 16c
        d["record_len2"] = l[199]           # i
        return d

    def _get_vlevel_headers(self, nfields):
        """ Read nfields vlevel headers, return a list of dicts. """
        return [self._get_vlevel_header() for i in range(nfields)]

    def _get_vlevel_header(self):
        """ Read a single vlevel header, return a dict. """
        fmt = '>i i 122i 4i 122f 5f i'
        l = struct.unpack(fmt, self.fileptr.read(struct.calcsize(fmt)))
        d = {}
        d["record_len1"] = l[0]         # i
        d["struct_id"] = l[1]           # i
        d["type"] = l[2:124]            # 122i
        d["unused_si32"] = l[124:128]   # 4i
        d["level"] = l[128:250]         # 122f
        d["unused_fl32"] = l[250:255]   # 5f
        d["record_len2"] = l[255]       # i
        return d

    def _get_chunk_headers(self, nchunks):
        """ Get nchunk chunk headers, return a list of dicts. """
        return [self._get_chunk_header() for i in range(nchunks)]

    def _get_chunk_header(self):
        """ Get a single chunk header, return a dict. """
        fmt = '>5i 2i 480c i'
        l = struct.unpack(fmt, self.fileptr.read(struct.calcsize(fmt)))
        d = {}
        d["record_len1"] = l[0]         # 5i: 1
        d["struct_id"] = l[1]           # 5i: 2
        d["chunk_id"] = l[2]            # 5i: 3
        d["chunk_data_offset"] = l[3]   # 5i: 4
        d["size"] = l[4]                # 5i: 5
        d["unused_si32"] = l[5:7]       # 2i
        d["info"] = ''.join(l[7:487]).strip('\x00')     # 480c
        d["record_len2"] = l[487]       # i
        return d

    def _get_chunks(self, debug=False):
        """ Get data in chunks, return radar_info, elevations, calib_info. """
        radar_info, elevations, calib_info = None, None, None
        for curr_chunk_header in self.chunk_headers:

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

        return radar_info, elevations, calib_info

    def _get_radar_info(self):
        """ Get the radar information, return dict. """
        fmt = '>12i 2i 22f 4f 40c 40c'
        l = struct.unpack(fmt, self.fileptr.read(struct.calcsize(fmt)))
        d = {}
        d["radar_id"] = l[0]                # 12i: 1
        d["radar_type"] = l[1]              # 12i: 2
        d["nfields"] = l[2]                 # 12i: 3
        d["ngates"] = l[3]                  # 12i: 4
        d["samples_per_beam"] = l[4]        # 12i: 5
        d["scan_type"] = l[5]               # 12i: 6
        d["scan_mode"] = l[6]               # 12i: 7
        d["nfields_current"] = l[7]         # 12i: 8
        d["field_flag"] = l[8]              # 12i: 9
        d["polarization"] = l[9]            # 12i: 10
        d["follow_mode"] = l[10]            # 12i: 11
        d["prf_mode"] = l[11]               # 12i: 12
        d["spare_ints"] = l[12:14]          # 2i
        d["radar_constant"] = l[14]         # 22f: 1
        d["altitude_km"] = l[15]            # 22f: 2
        d["latitude_deg"] = l[16]           # 22f: 3
        d["longitude_deg"] = l[17]          # 22f: 4
        d["gate_spacing_km"] = l[18]        # 22f: 5
        d["start_range_km"] = l[19]         # 22f: 6
        d["horiz_beam_width_deg"] = l[20]   # 22f: 7
        d["vert_beam_width_deg"] = l[21]    # 22f: 8
        d["pulse_width_us"] = l[22]         # 22f: 9
        d["prf_hz"] = l[23]                 # 22f: 10
        d["wavelength_cm"] = l[24]          # 22f: 11
        d["xmit_peak_pwr_watts"] = l[25]    # 22f: 12
        d["receiver_mds_dbm"] = l[26]       # 22f: 13
        d["receiver_gain_db"] = l[27]       # 22f: 14
        d["antenna_gain_db"] = l[28]        # 22f: 15
        d["system_gain_db"] = l[29]         # 22f: 16
        d["unambig_vel_mps"] = l[30]        # 22f: 17
        d["unambig_range_km"] = l[31]       # 22f: 18
        d["measXmitPowerDbmH_dbm"] = l[32]  # 22f: 19
        d["measXmitPowerDbmV_dbm"] = l[33]  # 22f: 20
        d["prt_s"] = l[34]                  # 22f: 21
        d["prt2_s"] = l[35]                 # 22f: 22
        d["spare_floats"] = l[36:40]        # 4f
        d["radar_name"] = ''.join(l[40:80]).strip('\x00')       # 40c
        d["scan_type_name"] = ''.join(l[80:120]).strip('\x00')   # 40c
        return d

    def _get_elevs(self, nbytes):
        """ Return an array of elevation read from current file position. """
        SIZE_FLOAT = 4.0
        nelevations = np.floor(nbytes / SIZE_FLOAT)
        fmt = '%df' % (nelevations)
        l = struct.unpack(fmt, self.fileptr.read(struct.calcsize(fmt)))
        return np.array(l)

    def _get_calib(self):
        """ Get the calibration information, return a dict. """
        fmt = '>16c 6i 51f 14f'
        l = struct.unpack(fmt, self.fileptr.read(struct.calcsize(fmt)))
        d = {}
        d["radar_name"] = ''.join(l[0:16]).strip('\x00')    # c16
        d["year"] = l[16]                           # 6i: 1
        d["month"] = l[17]                          # 6i: 2
        d["day"] = l[18]                            # 6i: 3
        d["hour"] = l[19]                           # 6i: 4
        d["minute"] = l[20]                         # 6i: 5
        d["second"] = l[21]                         # 6i: 6
        d["wavelength_cm"] = l[22]                  # 51f: 1
        d["beamwidth_h_deg"] = l[23]                # 51f: 2
        d["beamwidth_v_deg"] = l[24]                # 51f: 3
        d["antenna_gain_h_db"] = l[25]              # 51f: 4
        d["antenna_gain_v_db"] = l[26]              # 51f: 5
        d["pulse_width_us"] = l[27]                 # 51f: 6
        d["xmit_power_h_dbm"] = l[28]               # 51f: 7
        d["xmit_power_v_dbm"] = l[29]               # 51f: 8
        d["twoway_waveguide_loss_h_db"] = l[30]     # 51f: 9
        d["twoway_waveguide_loss_v_db"] = l[31]     # 51f: 10
        d["twoway_radome_loss_h_db"] = l[32]        # 51f: 11
        d["twoway_radome_loss_v_db"] = l[33]        # 51f: 12
        d["filter_loss_db"] = l[34]                 # 51f: 13
        d["radar_constant_h_db"] = l[35]            # 51f: 14
        d["radar_constant_v_db"] = l[36]            # 51f: 15
        d["noise_h_co_dbm"] = l[37]                 # 51f: 16
        d["noise_h_cx_dbm"] = l[38]                 # 51f: 17
        d["noise_v_co_dbm"] = l[39]                 # 51f: 18
        d["noise_v_cx_dbm"] = l[40]                 # 51f: 19
        d["rx_gain_h_co_dbm"] = l[41]               # 51f: 20
        d["rx_gain_h_cx_dbm"] = l[42]               # 51f: 21
        d["rx_gain_v_co_dbm"] = l[43]               # 51f: 22
        d["rx_gain_v_cx_dbm"] = l[44]               # 51f: 23
        d["zh1km_co_dbz"] = l[45]                   # 51f: 24
        d["zh1km_cx_dbz"] = l[46]                   # 51f: 25
        d["zv1km_co_dbz"] = l[47]                   # 51f: 26
        d["zv1km_cx_dbz"] = l[48]                   # 51f: 27
        d["sun_h_co_dbm"] = l[49]                   # 51f: 28
        d["sun_h_cx_dbm"] = l[50]                   # 51f: 29
        d["sun_v_co_dbm"] = l[51]                   # 51f: 30
        d["sun_v_cx_dbm"] = l[52]                   # 51f: 31
        d["noise_source_h_dbm"] = l[53]             # 51f: 32
        d["noise_source_v_dbm"] = l[54]             # 51f: 33
        d["power_meas_loss_h_db"] = l[55]           # 51f: 34
        d["power_meas_loss_v_db"] = l[56]           # 51f: 35
        d["coupler_fwd_loss_h_db"] = l[57]          # 51f: 36
        d["coupler_fwd_loss_v_db"] = l[58]          # 51f: 37
        d["zdr_bias_db"] = l[59]                    # 51f: 38
        d["ldr_h_bias_db"] = l[60]                  # 51f: 39
        d["ldr_v_bias_db"] = l[61]                  # 51f: 40
        d["system_phidp_deg"] = l[62]               # 51f: 41
        d["test_pulse_h_dbm"] = l[63]               # 51f: 42
        d["test_pulse_v_dbm"] = l[64]               # 51f: 43
        d["rx_slope_h_co_db"] = l[65]               # 51f: 44
        d["rx_slope_h_cx_db"] = l[66]               # 51f: 45
        d["rx_slope_v_co_db"] = l[67]               # 51f: 46
        d["rx_slope_v_cx_db"] = l[68]               # 51f: 47
        d["I0_h_co_dbm"] = l[69]                    # 51f: 48
        d["I0_h_cx_dbm"] = l[70]                    # 51f: 49
        d["I0_v_co_dbm"] = l[71]                    # 51f: 50
        d["I0_v_cx_dbm"] = l[72]                    # 51f: 51
        d["spare"] = l[73:87]                       # 14f
        return d

    def _get_compression_info(self):
        """ Get compression infomation, return a dict. """
        # the file pointer must be set at the correct location prior to call
        fmt = '>I I I I 2I'
        l = struct.unpack(fmt, self.fileptr.read(struct.calcsize(fmt)))
        d = {}
        d['magic_cookie'] = l[0]
        d['nbytes_uncompressed'] = l[1]
        d['nbytes_compressed'] = l[2]
        d['nbytes_coded'] = l[3]
        d['spare'] = l[4:6]
        return d

    def _get_sweep_info(self, nsweeps):
        """ Get sweep information, return a dict. """
        # the file pointer must be set at the correct location prior to call
        fmt = '%iI %iI' % (nsweeps, nsweeps)
        l = struct.unpack(fmt, self.fileptr.read(struct.calcsize(fmt)))
        d = {}
        d['vlevel_offsets'] = l[:nsweeps]
        d['vlevel_nbytes'] = l[nsweeps:2*nsweeps]
        return d

    # misc. methods

    def _calc_geometry(self):
        """ Calculate geometry, return az_deg, range_km, el_deg. """
        nsweeps = self.master_header['nsweeps']
        nrays = self.master_header['nrays']
        ngates = self.master_header['ngates']
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

    def _make_time_dict(self):
        """ Return a time dictionary. """
        t_base = datetime.datetime(1970, 1, 1, 00, 00)
        tb = datetime.timedelta(seconds=self.master_header['time_begin'])
        te = datetime.timedelta(seconds=self.master_header['time_end'])
        tc = datetime.timedelta(seconds=self.master_header['time_centroid'])
        return {'time_begin': t_base + tb, 'time_end': t_base + te,
                'time_centroid': t_base + tc}

    def _make_carts_dict(self):
        """ Return a carts dictionary, distances in meters. """

        # simple calculation involving 4/3 earth radius
        nsweeps = self.master_header['nsweeps']
        nrays = self.master_header['nrays']
        ngates = self.master_header['ngates']
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
