"""
Functions and classes common between MDV grid and radar files.

"""

# Code is adapted from Nitin Bharadwaj's Matlab code

import struct
import bz2
import gzip
import zlib
from io import BytesIO
import datetime

import numpy as np

from ..core.transforms import antenna_to_cartesian

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
TA_NOT_COMPRESSED = 0x2f2f2f2f
GZIP_COMPRESSED = 0xf7f7f7f7
GZIP_NOT_COMPRESSED = 0xf8f8f8f8
BZIP_COMPRESSED = 0xf3f3f3f3
BZIP_NOT_COMPRESSED = 0xf4f4f4f4
ZLIB_COMPRESSED = 0xf5f5f5f5
ZLIB_NOT_COMPRESSED = 0xf6f6f6f6
RL8_COMPRESSION = 0xfe0103fd

#  ***************** TRANSFORM *******************
DATA_TRANSFORM_NONE = 0  # None
DATA_TRANSFORM_LOG = 1  # Natural log

#  ***************** BIT ENCODING *******************
ENCODING_INT8 = 1  # unsigned 8 bit integer
ENCODING_INT16 = 2  # unsigned 16 bit integer
ENCODING_FLOAT32 = 5  # 32 bit IEEE floating point

#  ***************** CHUNK HEADER and DATA *******************
CHUNK_DSRADAR_PARAMS = 3
CHUNK_DSRADAR_ELEVATIONS = 7
CHUNK_DSRADAR_CALIB = 10
DS_LABEL_LEN = 40
NCHAR_DS_RADAR_PARAMS = 2 * DS_LABEL_LEN
DS_RADAR_CALIB_NAME_LEN = 16
DS_RADAR_CALIB_MISSING = -9999.0


class MdvFile(object):
    """
    A file object for MDV data.

    A `MdvFile` object stores metadata and data from a MDV file. Metadata is
    stored in dictionaries as attributes of the object, field data is
    stored as NumPy ndarrays as attributes with the field name. By default
    only metadata is read initially and field data must be read using the
    `read_a_field` or `read_all_fields` methods. This behavior can be changed
    by setting the `read_fields` parameter to True.

    Parameters
    ----------
    filename : str, file-like or None.
        Name of MDV file to read or file-like object pointing to the
        beginning of such a file. None can be used to initalize an object
        which can be used for writing mdv files.
    debug : bool
        True to print out debugging information, False to supress.
    read_fields : bool
        True to read all field during initalization, False (default) only
        reads metadata.

    Notes
    -----
    This class is not stable enough for general purpose MDV reading/writing,
    nor is that the intention, but with care it can provide sufficient
    read/write capacity.

    """
    # formats for use in the struct lib
    # mapper are used to convert tuples to dics, the formats is as follows:
    # (var_name, initial_position_in_tuple, final_position_in_tuple)
    master_header_fmt = b'>28i 8i i 5i 6f 3f 12f 512s 128s 128s i'
    master_header_mapper = [
        ("record_len1", 0, 1),
        ("struct_id", 1, 2),
        ("revision_number", 2, 3),
        ("time_gen", 3, 4),
        ("user_time", 4, 5),
        ("time_begin", 5, 6),
        ("time_end", 6, 7),
        ("time_centroid", 7, 8),
        ("time_expire", 8, 9),
        ("num_data_times", 9, 10),
        ("index_number", 10, 11),
        ("data_dimension", 11, 12),
        ("data_collection_type", 12, 13),
        ("user_data", 13, 14),
        ("native_vlevel_type", 14, 15),
        ("vlevel_type", 15, 16),
        ("vlevel_included", 16, 17),
        ("grid_orientation", 17, 18),
        ("data_ordering", 18, 19),
        ("nfields", 19, 20),
        ("max_nx", 20, 21),
        ("max_ny", 21, 22),
        ("max_nz", 22, 23),
        ("nchunks", 23, 24),
        ("field_hdr_offset", 24, 25),
        ("vlevel_hdr_offset", 25, 26),
        ("chunk_hdr_offset", 26, 27),
        ("field_grids_differ", 27, 28),
        ("user_data_si328", 28, 36),
        ("time_written", 36, 37),
        ("unused_si325", 37, 42),
        ("user_data_fl326", 42, 48),
        ("sensor_lon", 48, 49),
        ("sensor_lat", 49, 50),
        ("sensor_alt", 50, 51),
        ("unused_fl3212", 51, 63),
        ("data_set_info", 63, 64),
        ("data_set_name", 64, 65),
        ("data_set_source", 65, 66),
        ("record_len2", 66, 67)
    ]

    field_header_fmt = '>17i 10i 9i 4i f f 8f 12f 4f 5f 64s 16s 16s 16s 16s i'
    field_header_mapper = [
        ("record_len1", 0, 1),
        ("struct_id", 1, 2),
        ("field_code", 2, 3),
        ("user_time1", 3, 4),
        ("forecast_delta", 4, 5),
        ("user_time2", 5, 6),
        ("user_time3", 6, 7),
        ("forecast_time", 7, 8),
        ("user_time4", 8, 9),
        ("nx", 9, 10),
        ("ny", 10, 11),
        ("nz", 11, 12),
        ("proj_type", 12, 13),
        ("encoding_type", 13, 14),
        ("data_element_nbytes", 14, 15),
        ("field_data_offset", 15, 16),
        ("volume_size", 16, 17),
        ("user_data_si32", 17, 27),
        ("compression_type", 27, 28),
        ("transform_type", 28, 29),
        ("scaling_type", 29, 30),
        ("native_vlevel_type", 30, 31),
        ("vlevel_type", 31, 32),
        ("dz_constant", 32, 33),
        ("data_dimension", 33, 34),
        ("zoom_clipped", 34, 35),
        ("zoom_no_overlap", 35, 36),
        ("unused_si32", 36, 40),
        ("proj_origin_lat", 40, 41),
        ("proj_origin_lon", 41, 42),
        ("proj_param", 42, 50),
        ("vert_reference", 50, 51),
        ("grid_dx", 51, 52),
        ("grid_dy", 52, 53),
        ("grid_dz", 53, 54),
        ("grid_minx", 54, 55),
        ("grid_miny", 55, 56),
        ("grid_minz", 56, 57),
        ("scale", 57, 58),
        ("bias", 58, 59),
        ("bad_data_value", 59, 60),
        ("missing_data_value", 60, 61),
        ("proj_rotation", 61, 62),
        ("user_data_fl32", 62, 66),
        ("min_value", 66, 67),
        ("max_value", 67, 68),
        ("min_value_orig_vol", 68, 69),
        ("max_value_orig_vol", 69, 70),
        ("unused_fl32", 70, 71),
        ("field_name_long", 71, 72),
        ("field_name", 72, 73),
        ("units", 73, 74),
        ("transform", 74, 75),
        ("unused_char", 75, 76),
        ("record_len2", 76, 77)
    ]

    vlevel_header_fmt = '>i i 122i 4i 122f 5f i'
    vlevel_header_mapper = [
        ("record_len1", 0, 1),
        ("struct_id", 1, 2),
        ("type", 2, 124),
        ("unused_si32", 124, 128),
        ("level", 128, 250),
        ("unused_fl32", 250, 255),
        ("record_len2", 255, 256)
    ]

    chunk_header_fmt = '>5i 2i 480s i'
    chunk_header_mapper = [
        ("record_len1", 0, 1),
        ("struct_id", 1, 2),
        ("chunk_id", 2, 3),
        ("chunk_data_offset", 3, 4),
        ("size", 4, 5),
        ("unused_si32", 5, 7),
        ("info", 7, 8),
        ("record_len2", 8, 9)
    ]

    compression_info_fmt = '>I I I I 2I'
    compression_info_mapper = [
        ("magic_cookie", 0, 1),
        ("nbytes_uncompressed", 1, 2),
        ("nbytes_compressed", 2, 3),
        ("nbytes_coded", 3, 4),
        ("spare", 4, 6)
    ]

    radar_info_fmt = '>12i 2i 22f 4f 40s 40s'
    radar_info_mapper = [
        ("radar_id", 0, 1),
        ("radar_type", 1, 2),
        ("nfields", 2, 3),
        ("ngates", 3, 4),
        ("samples_per_beam", 4, 5),
        ("scan_type", 5, 6),
        ("scan_mode", 6, 7),
        ("nfields_current", 7, 8),
        ("field_flag", 8, 9),
        ("polarization", 9, 10),
        ("follow_mode", 10, 11),
        ("prf_mode", 11, 12),
        ("spare_ints", 12, 14),
        ("radar_constant", 14, 15),
        ("altitude_km", 15, 16),
        ("latitude_deg", 16, 17),
        ("longitude_deg", 17, 18),
        ("gate_spacing_km", 18, 19),
        ("start_range_km", 19, 20),
        ("horiz_beam_width_deg", 20, 21),
        ("vert_beam_width_deg", 21, 22),
        ("pulse_width_us", 22, 23),
        ("prf_hz", 23, 24),
        ("wavelength_cm", 24, 25),
        ("xmit_peak_pwr_watts", 25, 26),
        ("receiver_mds_dbm", 26, 27),
        ("receiver_gain_db", 27, 28),
        ("antenna_gain_db", 28, 29),
        ("system_gain_db", 29, 30),
        ("unambig_vel_mps", 30, 31),
        ("unambig_range_km", 31, 32),
        ("measXmitPowerDbmH_dbm", 32, 33),
        ("measXmitPowerDbmV_dbm", 33, 34),
        ("prt_s", 34, 35),
        ("prt2_s", 35, 36),
        ("spare_floats", 36, 40),
        ("radar_name", 40, 41),
        ("scan_type_name", 41, 42)
    ]

    calib_fmt = '>16s 6i 51f 14f'
    calib_mapper = [
        ("radar_name", 0, 1),
        ("year", 1, 2),
        ("month", 2, 3),
        ("day", 3, 4),
        ("hour", 4, 5),
        ("minute", 5, 6),
        ("second", 6, 7),
        ("wavelength_cm", 7, 8),
        ("beamwidth_h_deg", 8, 9),
        ("beamwidth_v_deg", 9, 10),
        ("antenna_gain_h_db", 10, 11),
        ("antenna_gain_v_db", 11, 12),
        ("pulse_width_us", 12, 13),
        ("xmit_power_h_dbm", 13, 14),
        ("xmit_power_v_dbm", 14, 15),
        ("twoway_waveguide_loss_h_db", 15, 16),
        ("twoway_waveguide_loss_v_db", 16, 17),
        ("twoway_radome_loss_h_db", 17, 18),
        ("twoway_radome_loss_v_db", 18, 19),
        ("filter_loss_db", 19, 20),
        ("radar_constant_h_db", 20, 21),
        ("radar_constant_v_db", 21, 22),
        ("noise_h_co_dbm", 22, 23),
        ("noise_h_cx_dbm", 23, 24),
        ("noise_v_co_dbm", 24, 25),
        ("noise_v_cx_dbm", 25, 26),
        ("rx_gain_h_co_dbm", 26, 27),
        ("rx_gain_h_cx_dbm", 27, 28),
        ("rx_gain_v_co_dbm", 28, 29),
        ("rx_gain_v_cx_dbm", 29, 30),
        ("zh1km_co_dbz", 30, 31),
        ("zh1km_cx_dbz", 31, 32),
        ("zv1km_co_dbz", 32, 33),
        ("zv1km_cx_dbz", 33, 34),
        ("sun_h_co_dbm", 34, 35),
        ("sun_h_cx_dbm", 35, 36),
        ("sun_v_co_dbm", 36, 37),
        ("sun_v_cx_dbm", 37, 38),
        ("noise_source_h_dbm", 38, 39),
        ("noise_source_v_dbm", 39, 40),
        ("power_meas_loss_h_db", 40, 41),
        ("power_meas_loss_v_db", 41, 42),
        ("coupler_fwd_loss_h_db", 42, 43),
        ("coupler_fwd_loss_v_db", 43, 44),
        ("zdr_bias_db", 44, 45),
        ("ldr_h_bias_db", 45, 46),
        ("ldr_v_bias_db", 46, 47),
        ("system_phidp_deg", 47, 48),
        ("test_pulse_h_dbm", 48, 49),
        ("test_pulse_v_dbm", 49, 50),
        ("rx_slope_h_co_db", 50, 51),
        ("rx_slope_h_cx_db", 51, 52),
        ("rx_slope_v_co_db", 52, 53),
        ("rx_slope_v_cx_db", 53, 54),
        ("I0_h_co_dbm", 54, 55),
        ("I0_h_cx_dbm", 55, 56),
        ("I0_v_co_dbm", 56, 57),
        ("I0_v_cx_dbm", 57, 58),
        ("spare", 58, 72)
    ]

    def __init__(self, filename, debug=False, read_fields=False):
        """ initalize """
        if debug:
            print("Opening file for reading: ", filename)
        if filename is None:
            # will creat empty structures, to be filled and written later
            self.fileptr = None
        elif hasattr(filename, 'read'):
            self.fileptr = filename
        else:
            self.fileptr = open(filename, 'rb')

        if debug:
            print("Getting master header")
        self.master_header = self._get_master_header()

        if debug:
            print("getting field headers")
        nfields = self.master_header['nfields']
        self.field_headers = self._get_field_headers(nfields)

        if debug:
            print("getting vlevel headers")
        self.vlevel_headers = self._get_vlevel_headers(nfields)

        if debug:
            print("getting chunk headers")
        nchunks = self.master_header['nchunks']
        self.chunk_headers = self._get_chunk_headers(nchunks)

        if debug:
            print("Getting Chunk Data")
        # store raw chunk data, used for unknown chunk information
        self.chunk_data = [None] * self.master_header['nchunks']
        self.radar_info, self.elevations, self.calib_info = self._get_chunks(
            debug)

        if self.master_header['nfields'] > 0:
            projections = {
                PROJ_LATLON: 'latlon',
                PROJ_LAMBERT_CONF: 'lambert_conform',
                PROJ_POLAR_STEREO: 'polar_stereographic',
                PROJ_FLAT: 'flat',
                PROJ_POLAR_RADAR: 'ppi',
                PROJ_OBLIQUE_STEREO: 'oblique_stereographic',
                PROJ_RHI_RADAR: 'rhi',
            }
            self.projection = projections[self.field_headers[0]['proj_type']]

        if debug:
            print("Making usable time objects")
        self.times = self._make_time_dict()

        if debug:
            print("indexing fields")
        self.fields = self._make_fields_list()

        self.fields_data = [None] * self.master_header["nfields"]

        if read_fields:
            if debug:
                print("Reading all fields")
            self.read_all_fields()
        return

    ##################
    # public methods #
    ##################
    def write(self, filename, debug=False):
        """
        Write object data to a MDV file.

        Note that the file is not explicitly closes, use x.close() to
        close file object when complete.

        Parameters
        ----------
        filename : str or file-like
            Filename or open file object to which data will be written.
        debug : bool, options
            True to print out debugging information, False to supress.

        """
        if debug:
            print("Opening file for writing:", filename)
        if hasattr(filename, 'write'):
            self.fileptr = filename
        else:
            self.fileptr = open(filename, 'wb')
        file_start = self.fileptr.tell()

        # write fields and chunk so that offsets can be calculated
        # headers are initially zeros
        headers_size = (1024 + (416 + 1024) * self.master_header["nfields"] +
                        512 * self.master_header["nchunks"])
        self.fileptr.write(b"\x00" * headers_size)

        if debug:
            print("Writing Fields Data")
        for ifield in range(self.master_header["nfields"]):
            self._write_a_field(ifield)

        # write chunks
        if debug:
            print("Writing Chunk Data")
        self._write_chunks(debug)

        # calculate offsets
        self._calc_file_offsets()
        self.fileptr.seek(file_start)

        # write headers
        if debug:
            print("Writing master header")
        self._write_master_header()

        if debug:
            print("Writing field headers")
        self._write_field_headers(self.master_header["nfields"])

        if debug:
            print("Writing vlevel headers")
        self._write_vlevel_headers(self.master_header["nfields"])

        if debug:
            print("Writing chunk headers")
        self._write_chunk_headers(self.master_header["nchunks"])

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
        if self.fields_data[fnum] is not None:
            if debug:
                print("Getting data from the object.")
            return self.fields_data[fnum]

        # field has not yet been read, populate the object and return
        if debug:
            print("No data found in object, populating")

        nz = field_header['nz']
        ny = field_header['ny']
        nx = field_header['nx']

        # read the header
        field_data = np.zeros([nz, ny, nx], dtype='float32')
        self.fileptr.seek(field_header['field_data_offset'])
        self._get_levels_info(nz)  # dict not used, but need to seek.

        for sw in range(nz):
            if debug:
                print("doing levels ", sw)

            # get the compressed level data
            compr_info = self._get_compression_info()
            if compr_info['magic_cookie'] == 0xfe0103fd:
                # Run length encoding only has 20 bytes of compression
                # information with slightly different order, back up
                # 4 bytes to read all of the compressed data.
                self.fileptr.seek(-4, 1)
                compr_data = self.fileptr.read(compr_info['spare'][0])
            else:
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
                raise NotImplementedError('encoding: ', encoding_type)

            # decompress the level data
            if compr_info['magic_cookie'] == GZIP_COMPRESSED:
                cd_fobj = BytesIO(compr_data)
                gzip_file_handle = gzip.GzipFile(fileobj=cd_fobj)
                decompr_data = gzip_file_handle.read(struct.calcsize(fmt))
                gzip_file_handle.close()
            elif compr_info['magic_cookie'] == ZLIB_COMPRESSED:
                decompr_data = zlib.decompress(compr_data)
            elif compr_info['magic_cookie'] == BZIP_COMPRESSED:
                decompr_data = bz2.decompress(compr_data)
            elif compr_info['magic_cookie'] == TA_NOT_COMPRESSED:
                decompr_data = compr_data
            elif compr_info['magic_cookie'] == GZIP_NOT_COMPRESSED:
                decompr_data = compr_data
            elif compr_info['magic_cookie'] == ZLIB_NOT_COMPRESSED:
                decompr_data = compr_data
            elif compr_info['magic_cookie'] == BZIP_NOT_COMPRESSED:
                decompr_data = compr_data
            elif compr_info['magic_cookie'] == RL8_COMPRESSION:
                # Run length encoding of 8-bit data
                # Compression info is in a different order, namely
                # int32 : RL8_FLAG (0xfe0103fd)
                # int32 : key
                # int32 : nbytes_array (bytes of encoded data with header)
                # int32 : nbytes_full (bytes of unencoded data, no header)
                # int32 : nbytes_coded (bytes of encoded data, no header)
                key = compr_info['nbytes_uncompressed']
                decompr_size = compr_info['nbytes_coded']
                decompr_data = _decode_rle8(compr_data, key, decompr_size)
            else:
                raise NotImplementedError('unknown compression mode')

            # read the decompressed data, reshape and mask
            sw_data = np.frombuffer(decompr_data, np_form).astype('float32')
            sw_data.shape = (ny, nx)
            mask = sw_data == field_header['bad_data_value']
            np.putmask(sw_data, mask, [np.NaN])

            # scale and offset the data, store in field_data
            scale = field_header['scale']
            bias = field_header['bias']
            field_data[sw, :, :] = sw_data * scale + bias

        # store data as object attribute and return
        self.fields_data[fnum] = field_data
        return field_data

    def read_all_fields(self):
        """ Read all fields, storing data to field name attributes. """
        for i in range(self.master_header['nfields']):
            self.read_a_field(i)

    def close(self):
        """ Close the MDV file. """
        self.fileptr.close()

    ###################
    # private methods #
    ###################

    def _write_a_field(self, fnum):
        """ write field number 'fnum' to mdv file. """
        # the file pointer must be set at the correct location prior to call
        field_header = self.field_headers[fnum]
        if field_header['compression_type'] != 3:
            import warnings
            warnings.warn(
                "compression_type not implemented, converting to zlib")
            field_header['compression_type'] = 3

        field_data = self.fields_data[fnum]
        nz = field_header['nz']
        # save file position
        field_start = self.fileptr.tell()
        # write zeros to vlevel_offsets and vlevel_nbytes, these will
        # replaced by the correct data later
        self.fileptr.write(b"\x00" * 4 * 2 * nz)
        field_size = 0
        vlevel_offsets = [0] * nz
        vlevel_nbytes = [0] * nz
        for sw in range(nz):
            vlevel_offsets[sw] = field_size
            # apply scaling, offset and masking to field data
            scale = field_header['scale']
            bias = field_header['bias']
            sw_data = (field_data[sw, :, :] - bias) / scale
            if hasattr(sw_data, 'mask'):
                sw_data = np.where(
                    sw_data.mask, field_header['bad_data_value'], sw_data)

            # encode field data to the correct type, round when necessary
            encoding_type = field_header['encoding_type']
            if encoding_type == ENCODING_INT8:
                sw_data = np.round(sw_data).astype(np.uint8)
                np_form = '>B'
            elif encoding_type == ENCODING_INT16:
                sw_data = np.round(sw_data).astype(np.uint16)
                np_form = '>H'
            elif encoding_type == ENCODING_FLOAT32:
                sw_data = sw_data.astype(np.float32)
                np_form = '>f'
            else:
                raise NotImplementedError('encoding: ', encoding_type)
            uncompr_data = np.array(sw_data, dtype=np_form).tostring()
            compr_data = zlib.compress(uncompr_data)
            if len(compr_data) > len(uncompr_data):
                magic = 0xf6f6f6f6
                compr_data = uncompr_data
            else:
                magic = 0xf5f5f5f5
            compr_info = {
                'magic_cookie': magic,
                'nbytes_uncompressed': len(uncompr_data),
                'nbytes_compressed': len(compr_data) + 24,
                'nbytes_coded': len(compr_data),
                'spare': [0, 0],
            }

            self._write_compression_info(compr_info)
            self.fileptr.write(compr_data)
            field_size = field_size + len(compr_data) + 24
            vlevel_nbytes[sw] = len(compr_data) + 24
        # rewrite vlevel_offsets and vlevel_nbytes with corrected data
        field_end = self.fileptr.tell()
        self.fileptr.seek(field_start)
        vlevels_dic = {'vlevel_offsets': vlevel_offsets,
                       'vlevel_nbytes': vlevel_nbytes}
        self._write_levels_info(nz, vlevels_dic)
        self.fileptr.seek(field_end)
        field_header["volume_size"] = field_size + 2 * 4 * nz

    # get_ methods for reading headers

    def _unpack_mapped_tuple(self, l, mapper):
        """ Create a dictionary from a tuple using a mapper. """
        d = {}
        for item in mapper:
            if item[2] == item[1] + 1:
                d[item[0]] = l[item[1]]
            else:
                d[item[0]] = l[item[1]:item[2]]
            if isinstance(d[item[0]], bytes):
                d[item[0]] = d[item[0]].decode('ascii').split('\x00', 1)[0]
        return d

    def _pack_mapped(self, d, mapper, fmt):
        """ Create a packed string using a mapper and format. """
        l = [0] * mapper[-1][2]
        for item in mapper:
            if item[2] == item[1] + 1:
                l[item[1]] = d[item[0]]
                if hasattr(l[item[1]], 'encode'):   # encode str/unicode
                    l[item[1]] = l[item[1]].encode('ascii')
            else:
                l[item[1]:item[2]] = d[item[0]]
        # cast to string as Python < 2.7.7 pack does not except unicode
        fmt = str(fmt)
        return struct.pack(fmt, *l)

    def _get_master_header(self):
        """ Read the MDV master header, return a dict. """
        # the file pointer must be set at the correct location prior to call
        if self.fileptr is None:
            l = [0] * self.master_header_mapper[-1][2]
            l[0] = 1016
            l[1] = 14142
            l[2] = 1
            l[9] = 1
            l[16] = 1
            l[17] = 1
            l[63] = ""
            l[64] = ""
            l[65] = ""
            l[66] = 1016
        else:
            l = struct.unpack(
                self.master_header_fmt,
                self.fileptr.read(struct.calcsize(self.master_header_fmt)))
        return self._unpack_mapped_tuple(l, self.master_header_mapper)

    def _write_master_header(self):
        """ Write the MDV master header. """
        # the file pointer must be set at the correct location prior to call
        d = self.master_header
        l = [0] * self.master_header_mapper[-1][2]
        for item in self.master_header_mapper:
            if item[2] == item[1] + 1:
                l[item[1]] = d[item[0]]
                if hasattr(l[item[1]], 'encode'):   # encode str/unicode
                    l[item[1]] = l[item[1]].encode('ascii')
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
        if self.fileptr is None:
            l = [0] * self.field_header_mapper[-1][2]
            l[0] = 408
            l[1] = 14143
            l[57] = 1  # scale
            l[71] = ""
            l[72] = ""
            l[73] = ""
            l[74] = ""
            l[75] = ""
            l[76] = 408
        else:
            l = struct.unpack(
                self.field_header_fmt,
                self.fileptr.read(struct.calcsize(self.field_header_fmt)))
        return self._unpack_mapped_tuple(l, self.field_header_mapper)

    def _write_field_header(self, d):
        """ Write the a single field header. """
        # the file pointer must be set at the correct location prior to call
        string = self._pack_mapped(
            d, self.field_header_mapper, self.field_header_fmt)
        self.fileptr.write(string)

    def _get_vlevel_headers(self, nfields):
        """ Read nfields vlevel headers, return a list of dicts. """
        # the file pointer must be set at the correct location prior to call
        return [self._get_vlevel_header() for i in range(nfields)]

    def _write_vlevel_headers(self, nfields):
        """ Write nfields vlevel headers. """
        # the file pointer must be set at the correct location prior to call
        for i in range(nfields):
            self._write_vlevel_header(self.vlevel_headers[i])

    def _get_vlevel_header(self):
        """ Read a single vlevel header, return a dict. """
        # the file pointer must be set at the correct location prior to call
        if self.fileptr is None:
            l = [0] * self.vlevel_header_mapper[-1][2]
            l[0] = 1016
            l[1] = 14144
            l[255] = 1016
        else:
            l = struct.unpack(
                self.vlevel_header_fmt,
                self.fileptr.read(struct.calcsize(self.vlevel_header_fmt)))
        return self._unpack_mapped_tuple(l, self.vlevel_header_mapper)

    def _write_vlevel_header(self, d):
        """  Write the a single vfield header. """
        # the file pointer must be set at the correct location prior to call
        string = self._pack_mapped(
            d, self.vlevel_header_mapper, self.vlevel_header_fmt)
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
        if self.fileptr is None:
            l = [0] * self.chunk_header_mapper[-1][2]
            l[0] = 504
            l[1] = 14145
            l[7] = ""
            l[8] = 504
        else:
            l = struct.unpack(
                self.chunk_header_fmt,
                self.fileptr.read(struct.calcsize(self.chunk_header_fmt)))
        return self._unpack_mapped_tuple(l, self.chunk_header_mapper)

    def _write_chunk_header(self, d):
        """  Write the a single chunk header. """
        # the file pointer must be set at the correct location prior to call
        string = self._pack_mapped(
            d, self.chunk_header_mapper, self.chunk_header_fmt)
        self.fileptr.write(string)

    def _get_chunks(self, debug=False):
        """ Get data in chunks, return radar_info, elevations, calib_info. """
        # the file pointer must be set at the correct location prior to call
        radar_info, elevations, calib_info = None, [], None
        for cnum, curr_chunk_header in enumerate(self.chunk_headers):

            chunk_id = curr_chunk_header['chunk_id']
            self.fileptr.seek(curr_chunk_header['chunk_data_offset'])

            if chunk_id == CHUNK_DSRADAR_PARAMS:
                if debug:
                    print('Getting radar info')
                radar_info = self._get_radar_info()

            elif chunk_id == CHUNK_DSRADAR_ELEVATIONS:
                if debug:
                    print('getting elevations')
                elevations = self._get_elevs(curr_chunk_header['size'])

            elif chunk_id == CHUNK_DSRADAR_CALIB:
                if debug:
                    print('getting cal')
                calib_info = self._get_calib()

            else:
                if debug:
                    print('getting unknown chunk %i' % chunk_id)
                self.chunk_data[cnum] = self._get_unknown_chunk(cnum)

        return radar_info, elevations, calib_info

    def _write_chunks(self, debug=False):
        """ write chunks data """
        # the file pointer must be set at the correct location prior to call
        for cnum, curr_chunk_header in enumerate(self.chunk_headers):
            chunk_id = curr_chunk_header['chunk_id']

            if chunk_id == CHUNK_DSRADAR_PARAMS:
                if debug:
                    print('writing radar info')
                self._write_radar_info(self.radar_info)

            elif chunk_id == CHUNK_DSRADAR_ELEVATIONS:
                if debug:
                    print('writing elevations')
                self._write_elevs(self.elevations)

            elif chunk_id == CHUNK_DSRADAR_CALIB:
                if debug:
                    print('writing cal')
                self._write_calib(self.calib_info)

            else:
                if debug:
                    print('writing unknown chunk %i' % chunk_id)
                self._write_unknown_chunk(self.chunk_data[cnum])

    def _get_radar_info(self):
        """ Get the radar information, return dict. """
        # the file pointer must be set at the correct location prior to call
        if self.fileptr is None:
            l = [0] * self.radar_info_mapper[-1][2]
            l[40] = ""
            l[41] = ""
        else:
            l = struct.unpack(
                self.radar_info_fmt,
                self.fileptr.read(struct.calcsize(self.radar_info_fmt)))
        return self._unpack_mapped_tuple(l, self.radar_info_mapper)

    def _write_radar_info(self, d):
        """  Write radar information. """
        # the file pointer must be set at the correct location prior to call
        string = self._pack_mapped(
            d, self.radar_info_mapper, self.radar_info_fmt)
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
        # cast to string as Python < 2.7.7 pack does not except unicode
        fmt = str(fmt)
        string = struct.pack(fmt, *l)
        self.fileptr.write(string)

    def _get_calib(self):
        """ Get the calibration information, return a dict. """
        # the file pointer must be set at the correct location prior to call
        if self.fileptr is None:
            l = [0] * self.calib_mapper[-1][2]
            l[0] = ""
        else:
            l = struct.unpack(
                self.calib_fmt,
                self.fileptr.read(struct.calcsize(self.calib_fmt)))
        return self._unpack_mapped_tuple(l, self.calib_mapper)

    def _write_calib(self, d):
        """  Write calibration information. """
        # the file pointer must be set at the correct location prior to call
        string = self._pack_mapped(
            d, self.calib_mapper, self.calib_fmt)
        self.fileptr.write(string)

    def _get_compression_info(self):
        """ Get compression infomation, return a dict. """
        # the file pointer must be set at the correct location prior to call
        if self.fileptr is None:
            l = [0] * self.compression_info_mapper[-1][2]
        else:
            l = struct.unpack(
                self.compression_info_fmt,
                self.fileptr.read(struct.calcsize(self.compression_info_fmt)))
        return self._unpack_mapped_tuple(l, self.compression_info_mapper)

    def _write_compression_info(self, d):
        """ Write compression infomation. """
        # the file pointer must be set at the correct location prior to call
        string = self._pack_mapped(
            d, self.compression_info_mapper, self.compression_info_fmt)
        self.fileptr.write(string)

    def _get_unknown_chunk(self, cnum):
        """ Get raw data from chunk. """
        # the file pointer must be set at the correct location prior to call
        size = self.chunk_headers[cnum]['size']
        return self.fileptr.read(size)

    def _write_unknown_chunk(self, data):
        """ Write raw data from chunk. """
        # the file pointer must be set at the correct location prior to call
        self.fileptr.write(data)

    def _get_levels_info(self, nlevels):
        """ Get nlevel information, return a dict. """
        # the file pointer must be set at the correct location prior to call
        fmt = '>%iI %iI' % (nlevels, nlevels)
        if self.fileptr:
            l = struct.unpack(fmt, self.fileptr.read(struct.calcsize(fmt)))
        else:
            l = [0] * 2 * nlevels
        d = {}
        d['vlevel_offsets'] = l[:nlevels]
        d['vlevel_nbytes'] = l[nlevels: 2 * nlevels]
        return d

    def _write_levels_info(self, nlevels, d):
        """ write levels information, return a dict. """
        # the file pointer must be set at the correct location prior to call
        fmt = '%iI %iI' % (nlevels, nlevels)
        l = d['vlevel_offsets'] + d['vlevel_nbytes']
        # cast to string as Python < 2.7.7 pack does not except unicode
        fmt = str(fmt)
        string = struct.pack(fmt, *l)
        self.fileptr.write(string)

    def _calc_file_offsets(self):
        """ Calculate file offsets. """
        self.master_header["field_hdr_offset"] = 1024
        self.master_header["vlevel_hdr_offset"] = (
            1024 + 416 * self.master_header["nfields"])
        self.master_header["chunk_hdr_offset"] = (
            1024 + (416 + 1024) * self.master_header["nfields"])

        file_pos = (self.master_header["chunk_hdr_offset"] +
                    512 * self.master_header["nchunks"])
        for i in range(self.master_header["nfields"]):
            self.field_headers[i]["field_data_offset"] = file_pos
            file_pos = file_pos + self.field_headers[i]["volume_size"]

        for i in range(self.master_header["nchunks"]):
            self.chunk_headers[i]["chunk_data_offset"] = file_pos
            file_pos = file_pos + self.chunk_headers[i]["size"]

    def _make_time_dict(self):
        """ Return a time dictionary. """
        t_base = datetime.datetime(1970, 1, 1, 00, 00)
        tb = datetime.timedelta(seconds=self.master_header['time_begin'])
        te = datetime.timedelta(seconds=self.master_header['time_end'])
        tc = datetime.timedelta(seconds=self.master_header['time_centroid'])
        return {'time_begin': t_base + tb, 'time_end': t_base + te,
                'time_centroid': t_base + tc}

    def _time_dict_into_header(self):
        """ Complete time information in master_header from the time dict. """
        self.master_header['time_begin'] = self._secs_since_epoch(
            self.times['time_begin'])
        self.master_header['time_end'] = self._secs_since_epoch(
            self.times['time_end'])
        self.master_header['time_centroid'] = self._secs_since_epoch(
            self.times['time_centroid'])

    def _secs_since_epoch(self, dt):
        """ Return the number of seconds since the epoch for a datetime. """
        epoch = datetime.datetime(1970, 1, 1, 0, 0)
        td = dt - epoch
        # use td.total_seconds() in Python 2.7+
        return int((td.microseconds + (td.seconds + td.days * 24 * 3600) *
                    10**6) / 10**6)

    # misc. methods
    # XXX move some where else, there are not general mdv operations
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
            el_deg = grid_miny + np.arange(nrays) * grid_dy
            az_deg = self.vlevel_headers[0]['level'][0:nsweeps]
        elif self.field_headers[0]['proj_type'] == PROJ_POLAR_RADAR:
            az_deg = grid_miny + np.arange(nrays) * grid_dy
            el_deg = self.vlevel_headers[0]['level'][0:nsweeps]
        else:
            proj_type = self.field_headers[0]['proj_type']
            message = ("Unsupported projection type: %i, " % (proj_type) +
                       "is MDV file in antenna coordinates?")
            raise NotImplementedError(message)

        return az_deg, range_km, el_deg

    def _make_carts_dict(self):
        """ Return a carts dictionary, distances in meters. """
        az_deg, range_km, el_deg = self._calc_geometry()
        # simple calculation involving 4/3 earth radius
        nsweeps = self.master_header['max_nz']
        nrays = self.master_header['max_ny']
        ngates = self.master_header['max_nx']
        xx = np.empty([nsweeps, nrays, ngates], dtype=np.float32)
        yy = np.empty([nsweeps, nrays, ngates], dtype=np.float32)
        zz = np.empty([nsweeps, nrays, ngates], dtype=np.float32)

        if self.projection == 'rhi':
            rg, ele = np.meshgrid(range_km, el_deg)
            rg = np.array(rg, dtype=np.float64)
            ele = np.array(ele, dtype=np.float64)
            for aznum in range(nsweeps):
                azg = np.ones(rg.shape, dtype=np.float64) * az_deg[aznum]
                x, y, z = antenna_to_cartesian(rg, azg, ele)
                zz[aznum, :, :] = z
                xx[aznum, :, :] = x
                yy[aznum, :, :] = y

        elif self.projection == 'ppi':
            rg, azg = np.meshgrid(range_km, az_deg)
            rg = np.array(rg, dtype=np.float64)
            azg = np.array(azg, dtype=np.float64)
            for elnum in range(nsweeps):
                ele = np.ones(rg.shape, dtype=np.float64) * el_deg[elnum]
                x, y, z = antenna_to_cartesian(rg, azg, ele)
                zz[elnum, :, :] = z
                xx[elnum, :, :] = x
                yy[elnum, :, :] = y

        return {'x': xx, 'y': yy, 'z': zz}

    def _make_fields_list(self):
        """ Return a list of fields. """
        fh = self.field_headers
        return [fh[i]['field_name'] for i in range(len(fh))]


def _decode_rle8(compr_data, key, decompr_size):
    """ Decode 8-bit MDV run length encoding. """
    # Encoding is described in section 7 of:
    # http://rap.ucar.edu/projects/IHL/RalHtml/protocols/mdv_file/
    # This function would benefit greate by being rewritten in Cython
    data = np.frombuffer(compr_data, dtype='>B')
    out = np.empty((decompr_size, ), dtype='uint8')
    data_ptr = 0
    out_ptr = 0
    while data_ptr != len(data):
        v = data[data_ptr]
        if v != key:  # not encoded
            out[out_ptr] = v
            data_ptr += 1
            out_ptr += 1
        else:   # run length encoded
            count = data[data_ptr+1]
            value = data[data_ptr+2]
            out[out_ptr:out_ptr+count] = value
            data_ptr += 3
            out_ptr += count
    return out.tostring()


class _MdvVolumeDataExtractor(object):
    """
    Class facilitating on demand extraction of data from a MDV file.

    Parameters
    ----------
    mdvfile : MdvFile
        Open MdvFile object to extract data from.
    field_num : int
        Field number of data to be extracted.
    fillvalue : int
        Value used to fill masked values in the returned array.
    two_dims : bool.
        True to combine the first and second dimension of the array when
        returning the data, False will return a three dimensional array.

    """

    def __init__(self, mdvfile, field_num, fillvalue, two_dims=True):
        """ initialize the object. """
        self.mdvfile = mdvfile
        self.field_num = field_num
        self.fillvalue = fillvalue
        self.two_dims = two_dims

    def __call__(self):
        """ Return an array containing data from the referenced volume. """
        # grab data from MDV object, mask and reshape
        data = self.mdvfile.read_a_field(self.field_num)
        data[np.where(np.isnan(data))] = self.fillvalue
        data[np.where(data == 131072)] = self.fillvalue
        data = np.ma.masked_equal(data, self.fillvalue)
        if self.two_dims:
            data.shape = (data.shape[0] * data.shape[1], data.shape[2])
        return data
