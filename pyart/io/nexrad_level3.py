"""
pyart.io.nexrad_level3
======================

Class for reading data from NEXRAD Level 3 files.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    NEXRADLevel3File

.. autosummary::
    :toctree: generated/

    nexrad_level3_message_code
    _datetime_from_mdate_mtime
    _structure_size
    _unpack_from_buf
    _unpack_structure
    _int16_to_float16


"""

from __future__ import division

# This file is part of the Py-ART, the Python ARM Radar Toolkit
# https://github.com/ARM-DOE/pyart

# Care has been taken to keep this file free from extraneous dependancies
# so that it can be used by other projects with no/minimal modification.

# Please feel free to use this file in other project provided the license
# below is followed.  Keeping the above comment lines would also be helpful
# to direct other back to the Py-ART project and the source of this file.


LICENSE = """
Copyright (c) 2013, UChicago Argonne, LLC
All rights reserved.

Copyright 2013 UChicago Argonne, LLC. This software was produced under U.S.
Government contract DE-AC02-06CH11357 for Argonne National Laboratory (ANL),
which is operated by UChicago Argonne, LLC for the U.S. Department of Energy.
The U.S. Government has rights to use, reproduce, and distribute this
software. NEITHER THE GOVERNMENT NOR UCHICAGO ARGONNE, LLC MAKES ANY
WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS
SOFTWARE. If software is modified to produce derivative works, such modified
software should be clearly marked, so as not to confuse it with the version
available from ANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions
are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    * Neither the name of UChicago Argonne, LLC, Argonne National
      Laboratory, ANL, the U.S. Government, nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY UCHICAGO ARGONNE, LLC AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL UCHICAGO ARGONNE, LLC OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import bz2
from datetime import datetime, timedelta
import struct

import numpy as np


class NEXRADLevel3File(object):
    """
    A Class for accessing data in NEXRAD Level III (3) files.

    Attributes
    ----------
    text_header : dic
        File textual header.
    msg_header : dic
        Message header.
    prod_descr : dic
        Product description.
    symbology_header : dict
        Symbology header.
    packet_header : dict
        Radial data array packet header.
    radial_headers : list of dicts
        List of radials headers.
    raw_data : array
        Raw unscaled, unmasked data.
    data : array
        Scaled, masked radial data.
    _fh : file-like
        File like object from which data is read.

    """

    def __init__(self, filename):
        """ initalize the object. """
        # read the entire file into memory
        if hasattr(filename, 'read'):
            fhandle = filename
        else:
            fhandle = open(filename, 'rb')
        buf = fhandle.read()    # string buffer containing file data
        self._fh = fhandle

        # Text header
        # Format of Text header is SDUSXX KYYYY DDHHMM\r\r\nAAABBB\r\r\n
        # Sometime additional padding is present before the Text header
        record_padding = buf.find(b'SDUS')
        if record_padding == -1:
            raise ValueError('Not a valid NEXRAD Level 3 file.')
        self.text_header = buf[:30 + record_padding]
        bpos = 30 + record_padding      # current reading position in buffer

        # Read and decode 18 byte Message Header Block
        self.msg_header = _unpack_from_buf(buf, bpos, MESSAGE_HEADER)
        if self.msg_header['code'] not in SUPPORTED_PRODUCTS:
            code = self.msg_header['code']
            raise NotImplementedError(
                'Level3 product with code %i is not supported' % (code))
        bpos += 18

        # Read and decode 102 byte Product Description Block
        self.prod_descr = _unpack_from_buf(buf, bpos, PRODUCT_DESCRIPTION)
        bpos += 102

        # uncompressed symbology block if necessary
        if buf[bpos:bpos+2] == b'BZ':
            buf2 = bz2.decompress(buf[bpos:])
        else:
            buf2 = buf[bpos:]

        self._read_symbology_block(buf2)

    def close(self):
        """ Close the file. """
        self._fh.close()

    def _read_symbology_block(self, buf2):
        """ Read symbology block. """
        # Read and decode symbology header
        self.symbology_header = _unpack_from_buf(buf2, 0, SYMBOLOGY_HEADER)

        # Read radial packets
        packet_code = struct.unpack('>h', buf2[16:18])[0]
        assert packet_code in SUPPORTED_PACKET_CODES
        self.packet_header = _unpack_from_buf(buf2, 16, RADIAL_PACKET_HEADER)
        self.radial_headers = []
        nbins = self.packet_header['nbins']
        nradials = self.packet_header['nradials']
        nbytes = _unpack_from_buf(buf2, 30, RADIAL_HEADER)['nbytes']
        if packet_code == 16 and nbytes != nbins:
            nbins = nbytes  # sometimes these do not match, use nbytes
        self.raw_data = np.empty((nradials, nbins), dtype='uint8')
        pos = 30

        for radial in self.raw_data:
            radial_header = _unpack_from_buf(buf2, pos, RADIAL_HEADER)
            pos += 6
            if packet_code == 16:
                radial[:] = np.frombuffer(buf2[pos:pos+nbins], '>u1')
                pos += radial_header['nbytes']
            else:
                assert packet_code == AF1F
                # decode run length encoding
                rle_size = radial_header['nbytes'] * 2
                rle = np.frombuffer(buf2[pos:pos+rle_size], dtype='>u1')
                colors = np.bitwise_and(rle, 0b00001111)
                runs = np.bitwise_and(rle, 0b11110000) // 16
                radial[:] = np.repeat(colors, runs)
                pos += rle_size
            self.radial_headers.append(radial_header)

    def get_location(self):
        """ Return the latitude, longitude and height of the radar. """
        latitude = self.prod_descr['latitude'] * 0.001
        longitude = self.prod_descr['longitude'] * 0.001
        height = self.prod_descr['height']
        return latitude, longitude, height

    def get_azimuth(self):
        """ Return an array of starting azimuth angles in degrees. """
        azimuths = [d['angle_start'] for d in self.radial_headers]
        return np.array(azimuths, dtype='float32') * 0.1

    def get_range(self):
        """ Return an array of gate range spacing in meters. """
        nbins = self.raw_data.shape[1]
        first_bin = self.packet_header['first_bin']
        range_scale = (self.packet_header['range_scale'] *
                       PRODUCT_RANGE_RESOLUTION[self.msg_header['code']])
        return np.arange(nbins, dtype='float32') * range_scale + first_bin

    def get_elevation(self):
        """ Return the sweep elevation angle in degrees. """
        hw30 = self.prod_descr['halfwords_30']
        elevation = struct.unpack('>h', hw30)[0] * 0.1
        return elevation

    def get_volume_start_datetime(self):
        """ Return a datetime of the start of the radar volume. """
        return _datetime_from_mdate_mtime(self.prod_descr['vol_scan_date'],
                                          self.prod_descr['vol_scan_time'])

    def get_data(self):
        """ Return an masked array containing the field data. """
        msg_code = self.msg_header['code']
        threshold_data = self.prod_descr['threshold_data']

        if msg_code in _8_OR_16_LEVELS:
            mdata = self._get_data_8_or_16_levels()

        elif msg_code in [134]:
            mdata = self._get_data_msg_134()

        elif msg_code in [94, 99, 182, 186]:
            hw31, hw32 = np.frombuffer(threshold_data[:4], '>i2')
            data = (self.raw_data - 2) * (hw32/10.) + hw31/10.
            mdata = np.ma.array(data, mask=self.raw_data < 2)

        elif msg_code in [32]:
            hw31, hw32 = np.frombuffer(threshold_data[:4], '>i2')
            data = (self.raw_data) * (hw32/10.) + hw31/10.
            mdata = np.ma.array(data, mask=self.raw_data < 2)

        elif msg_code in [138]:
            hw31, hw32 = np.frombuffer(threshold_data[:4], '>i2')
            data = self.raw_data * (hw32/100.) + hw31/100.
            mdata = np.ma.array(data)

        elif msg_code in [159, 161, 163]:
            scale, offset = np.frombuffer(threshold_data[:8], '>f4')
            data = (self.raw_data - offset) / (scale)
            mdata = np.ma.array(data, mask=self.raw_data < 2)

        elif msg_code in [170, 172, 173, 174, 175]:
            # units are 0.01 inches
            scale, offset = np.frombuffer(threshold_data[:8], '>f4')
            data = (self.raw_data - offset) / (scale) * 0.01
            mdata = np.ma.array(data, mask=self.raw_data < 1)

        elif msg_code in [165, 177]:
            # Corresponds to classifications in table on page 3-37
            mdata = np.ma.masked_equal(self.raw_data, 0)

        elif msg_code in [135]:
            mdata = np.ma.array(self.raw_data - 2, mask=self.raw_data <= 1)
            mdata[self.raw_data >= 128] -= np.uint8(128)

        else:
            assert msg_code in [34]
            # There does not seem to be any discussion on what this product
            # contains.
            mdata = np.ma.array(self.raw_data.copy())

        return mdata.astype('float32')

    def _get_data_8_or_16_levels(self):
        """ Return a masked array for products with 8 or 16 data levels. """
        thresh = np.frombuffer(self.prod_descr['threshold_data'], '>B')
        flags = thresh[::2]
        values = thresh[1::2]

        sign = np.choose(np.bitwise_and(flags, 0x01), [1, -1])
        bad = np.bitwise_and(flags, 0x80) == 128
        scale = 1.
        if flags[0] & 2**5:
            scale = 1/20.
        if flags[0] & 2**4:
            scale = 1/10.

        data_levels = values * sign * scale
        data_levels[bad] = -999 # sentinal for bad data points

        data = np.choose(self.raw_data, data_levels)
        mdata = np.ma.masked_equal(data, -999)
        return mdata

    def _get_data_msg_134(self):
        """ Return a masked array for product with message code 134. """
        hw31, hw32, hw33, hw34, hw35 = np.frombuffer(
            self.prod_descr['threshold_data'][:10], '>i2')
        linear_scale = _int16_to_float16(hw31)
        linear_offset = _int16_to_float16(hw32)
        log_start = hw33
        log_scale = _int16_to_float16(hw34)
        log_offset = _int16_to_float16(hw35)
        # linear scale data
        data = np.zeros(self.raw_data.shape, dtype=np.float32)
        lin = self.raw_data < log_start
        data[lin] = ((self.raw_data[lin] - linear_offset) / (linear_scale))
        # log scale data
        log = self.raw_data >= log_start
        data[log] = np.exp((self.raw_data[log] - log_offset) / (log_scale))
        mdata = np.ma.masked_array(data, mask=self.raw_data < 2)
        return mdata


def _datetime_from_mdate_mtime(mdate, mtime):
    """ Returns a datetime for a given message date and time. """
    epoch = datetime.utcfromtimestamp(0)
    return epoch + timedelta(days=mdate - 1, seconds=mtime)


def _structure_size(structure):
    """ Find the size of a structure in bytes. """
    return struct.calcsize('>' + ''.join([i[1] for i in structure]))


def _unpack_from_buf(buf, pos, structure):
    """ Unpack a structure from a buffer. """
    size = _structure_size(structure)
    return _unpack_structure(buf[pos:pos + size], structure)


def _unpack_structure(string, structure):
    """ Unpack a structure from a string """
    fmt = '>' + ''.join([i[1] for i in structure]) # NEXRAD is big-endian
    lst = struct.unpack(fmt, string)
    return dict(zip([i[0] for i in structure], lst))


def nexrad_level3_message_code(filename):
    """ Return the message (product) code for a NEXRAD Level 3 file. """
    fhl = open(filename, 'r')
    buf = fhl.read(48)
    fhl.close()
    msg_header = _unpack_from_buf(buf, 30, MESSAGE_HEADER)
    return msg_header['code']


# NEXRAD Level III file structures, sizes, and static data
# The deails on these structures are documented in:
# "INTERFACE CONTROL DOCUMENT FOR THE RPG TO CLASS 1 USER" RPG Build 13.0
# Document Number 2620001T
# Tables and page number refer to those in this document.


def _int16_to_float16(val):
    """ Convert a 16 bit interger into a 16 bit float. """
    # NEXRAD Level III float16 format defined on page 3-33.
    # Differs from IEEE 768-2008 format so np.float16 cannot be used.
    sign = (val & 0b1000000000000000) / 0b1000000000000000
    exponent = (val & 0b0111110000000000) / 0b0000010000000000
    fraction = (val & 0b0000001111111111)
    if exponent == 0:
        return (-1)**sign * 2 * (0 + (fraction/2**10.))
    else:
        return (-1)**sign * 2**(exponent-16) * (1 + fraction/2**10.)


_8_OR_16_LEVELS = [19, 20, 25, 27, 28, 30, 56, 78, 79, 80, 169, 171, 181]

PRODUCT_RANGE_RESOLUTION = {
    19: 1.,     # 124 nm
    20: 2.,     # 248 nm
    25: 0.25,   # 32 nm
    27: 1.,
    28: 0.25,
    30: 1.,
    32: 1.,
    34: 1.,
    56: 1.,
    78: 1.,
    79: 1.,
    80: 1.,
    94: 1.,
    99: 0.25,
    134: 1000.,
    135: 1000.,
    138: 1.,
    159: 0.25,
    161: 0.25,
    163: 0.25,
    165: 0.25,
    169: 1.,
    170: 1.,
    171: 1.,
    172: 1.,
    173: 1.,
    174: 1.,
    175: 1.,
    177: 0.25,
    181: 150.,
    182: 150.,
    186: 300.,
}


# format of structure elements
# Figure E-1, page E-1
BYTE = 'B'      # not in table but used in Product Description
INT2 = 'h'
INT4 = 'i'
UINT4 = 'I'
REAL4 = 'f'

# 3.3.1 Graphic Product Messages

# Graphic Product Message: Message Header Block
# 18 bytes, 9 halfwords
# Figure 3-3, page 3-7.
MESSAGE_HEADER = (
    ('code', INT2),     # message code
    ('date', INT2),     # date of message, days since 1 Jan, 1970
    ('time', INT4),     # time of message, seconds since midnight
    ('length', INT4),   # length of message in bytes
    ('source', INT2),   # Source ID
    ('dest', INT2),     # Destination ID
    ('nblocks', INT2),  # Number of blocks in the message (inclusive)
)

# Graphic Product Message: Product Description Block
# Description: section 3.3.1.1, page 3-3
# 102 bytes, 51 halfwords (halfwords 10-60)
# Figure 3-6, pages 3-31 and 3-32
PRODUCT_DESCRIPTION = (
    ('divider', INT2),          # Delineate blocks, -1
    ('latitude', INT4),         # Latitude of radar, degrees, + for north
    ('longitude', INT4),        # Longitude of radar, degrees, + for east
    ('height', INT2),           # Height of radar, feet abouve mean sea level
    ('product_code', INT2),     # NEXRAD product code
    ('operational_mode', INT2),  # 0 = Maintenance, 1 = Clean Air, 2 = Precip
    ('vcp', INT2),              # Volume Coverage Pattern of scan strategy
    ('sequence_num', INT2),     # Sequence Number of the request.
    ('vol_scan_num', INT2),     # Volume Scan number, 1 to 80.
    ('vol_scan_date', INT2),    # Volume Scan start date, days since 1/1/1970
    ('vol_scan_time', INT4),    # Volume Scan start time, sec since midnight
    ('product_date', INT2),     # Product Generation Date, days since 1/1/1970
    ('product_time', INT4),     # Product Generation Time, sec since midnight
    ('halfwords_27_28', '4s'),  # Product dependent parameters 1 and 2
    ('elevation_num', INT2),    # Elevation number within volume scan
    ('halfwords_30', '2s'),     # Product dependent parameter 3
    ('threshold_data', '32s'),  # Data to determine threshold level values
    ('halfwords_47_53', '14s'),  # Product dependent parameters 4-10
    ('version', BYTE),          # Version, 0
    ('spot_blank', BYTE),       # 1 = Spot blank ON, 0 = Blanking OFF
    ('offet_symbology', INT4),  # halfword offset to Symbology block
    ('offset_graphic', INT4),   # halfword offset to Graphic block
    ('offset_tabular', INT4)    # halfword offset to Tabular block
)

# Graphic Product Message: Product Symbology Block
# Description
# 16 byte header
# Figure 3-6 (Sheet 8), pages 3-40

SYMBOLOGY_HEADER = (
    ('divider', INT2),          # Delineate blocks, -1
    ('id', INT2),               # Block ID, 1
    ('block_length', INT4),     # Length of block in bytes
    ('layers', INT2),           # Number of data layers
    ('layer_divider', INT2),    # Delineate data layers, -1
    ('layer_length', INT4)      # Length of data layer in bytes
    # Display data packets
)

# Digital Radial Data Array Packet - Packet Code 16 (Sheet 2)
# Figure 3-11c (Sheet 1 and 2), page 3-120.
# and
# Radial Data Packet - Packet Code AF1F
# Figure 3-10 (Sheet 1 and 2), page 3-113.
AF1F = -20705       # struct.unpack('>h', 'AF1F'.decode('hex'))
SUPPORTED_PACKET_CODES = [16, AF1F]       # elsewhere
RADIAL_PACKET_HEADER = (
    ('packet_code', INT2),      # Packet Code, Type 16
    ('first_bin', INT2),        # Location of first range bin.
    ('nbins', INT2),            # Number of range bins.
    ('i_sweep_center', INT2),   # I coordinate of center of sweep.
    ('j_sweep_center', INT2),   # J coordinate of center of sweep.
    ('range_scale', INT2),      # Range Scale factor
    ('nradials', INT2)          # Total number of radials in the product
)

RADIAL_HEADER = (
    ('nbytes', INT2),           # Number of bytes in the radial.
    ('angle_start', INT2),      # Starting angle at which data was collected.
    ('angle_delta', INT2)       # Delta angle from previous radial.
)

# A list of the NEXRAD Level 3 Product supported by this module taken
# from the "Message Code for Products" Table III pages 3-15 to 3-22
# All the supported products have a Radial Image Message format.
#   Code    # Product Name
#   -----   -----------------------
SUPPORTED_PRODUCTS = [
    19,     # Base Reflectivity
    20,     # Base Reflectivity
    25,     # Base Velocity
    27,     # Base Velocity
    28,     # Base Spectrum Width
    30,     # Base Spectrum Width
    32,     # Digital Hybrid Scan
    34,     # Clutter Filter Control
    56,     # Storm Relative Mean
            # Radial Velocity
    78,     # Surface Rainfall Accum.
            # (1 hr)
    79,     # Surface Rainfall Accum.
            # (3 hr)
    80,     # Storm Total Rainfall
            # Accumulation
    94,     # Base Reflectivity Data
            # Array
    99,     # Base Velocity Data
            # Array
    134,    # High Resolution VIL
    135,    # Enhanced Echo Tops
    138,    # Digital Storm Total
            # Precipitation
    159,    # Digital Differential
            # Reflectivity
    161,    # Digital Correlation
            # Coefficient
    163,    # Digital Specific
            # Differential Phase
    165,    # Digital Hydrometeor
            # Classification
    169,    # One Hour Accumulation
    170,    # Digital Accumulation
            # Array
    171,    # Storm Total
            # Accumulation
    172,    # Digital Storm Total
            # Accumulation
    173,    # Digital User-Selectable
            # Accumulation
    174,    # Digital One-Hour
            # Difference Accumulation
    175,    # Digital Storm Total
            # Difference Accumulation
    177,    # Hybrid Hydrometeor
            # Classification
    181,    # Base Reflectivity
    182,    # Base Velocity
    186,    # Base Reflectivity
]

# It should be possible to add support for these NEXRAD Level 3 products
# as they are of Radial Image message format.  No examples of these files
# could be found to test on so support does not exist yet in this module.
#    Code    # Product Name
#    -----   -----------------------
#    16,     # Base Reflectivity
#    17,     # Base Reflectivity
#    18,     # Base Reflectivity
#    21,     # Base Reflectivity
#    22,     # Base Velocity
#    23,     # Base Velocity
#    24,     # Base Velocity
#    26,     # Base Velocity
#    29,     # Base Spectrum Width
#    31,     # User Selectable Storm Total Precipitation Reflectivty
#    33,     # Hybrid Scan Reflectivty
#    55,     # Storm Relative Mean Radial Velocity
#    93,     # ITWS Digital Base Velocity
#    132,    # Clutter Likelihood Reflectivity
#    133,    # Clutter Likelihood Doppler
#    137,    # User Selectable Layer Composite Reflectivity
#    144,    # One-hour Snow Water Equivalent
#    145,    # One-hour Snow Depth
#    146,    # Storm Total Snow Water Equivalent
#    147,    # Storm Total Snow Depth
#    150,    # User Selectable Snow Water Equivalent
#    151,    # User Selectable Snow Depth
#    153,    # Super Resolution Reflectivity Data Array
#    154,    # Super Resolution Velocity Data Array
#    155,    # Super Resolution Spectrum Width Data Array
#    158,    # Differential Reflectivity
#    160,    # Correlation Coefficient
#    162,    # Specific Differential Phase
#    164,    # Hydrometeor Classification
#    194,    # Base Reflectivity Data Array (DoD Version)
#    195,    # Digital Reflectivity, DQA-Edited Data Array
#    199,    # Base Velocity Data Array (DoD Version)
#    183,    # TDWR Base Velocity
#    180,    # TDWR Base Reflectivity
#    185,    # TDWR Base Spectrum Width
#    187,    # TDWR Base Reflectivity

# No support for these NEXRAD Level 3 products is planned as they
# do not have a Radial Image message format.

#   Code    # Product Name             Message format
#   -----   -----------------------    ---------------
#   35,     # Composite Reflectivity  Raster Image/Nongeographic Alpha
#   36,     # Composite Reflectivity  Raster Image/Nongeographic Alpha
#   37,     # Composite Reflectivity  Raster Image/Nongeographic Alpha
#   38,     # Composite Reflectivity  Raster Image/Nongeographic Alpha
#   39,     # Spare
#   40,     # Spare
#   41,     # Echo Tops               Raster Image
#   42,     # Spare
#   43,     # Spare
#   44,     # Spare
#   45,     # Spare
#   46,     # Spare
#   47,     # Spare
#   48,     # VAD Wind Profile        Non-geographic Alphanumeric
#   49,     # Spare
#   50,     # Cross Section           Raster Image (Reflectivity)
#           # (Reflectivity)
#   51,     # Cross Section           Raster Image (Velocity)
#           # (Velocity)
#   52,     # Spare
#   53,     # Spare
#   54,     # Reserved
#   57,     # Vertically Integrated   Raster Image
#           # Liquid
#   58,     # Storm Tracking          Non-geographic Alpha
#           # Information
#   59,     # Hail Index              Non-geographic Alpha
#   60,     # Spare
#   61,     # Tornado Vortex          Geographic and Non-geographic
#           # Signature               Alphanumeric
#   62,     # Storm Structure         Alphanumeric
#   63,     # Layer Composite         Raster Image (Layer 1 Average)
#           # Reflectivity
#   64,     # Layer Composite         Raster Image (Layer 2 Average)
#           # Reflectivity
#   65,     # Layer Composite         Raster Image (Layer 1 Maximum)
#           # Reflectivity
#   66,     # Layer Composite         Raster Image (Layer 2 Maximum)
#           # Reflectivity
#   67,     # Layer Composite         Raster Image
#           # Reflectivity - AP Removed
#   68,     # Spare
#   69,     # Spare
#   70,     # Spare
#   71,     # Spare
#   72,     # Spare
#   73,     # User Alert Message      Alphanumeric
#   74,     # Radar Coded Message     Alphanumeric
#   75,     # Free Text Message       Alphanumeric
#   76,     # Reserved for internal PUP use.
#   81,     # Hourly Digital          Raster Image / Alphanumeric
#           # Precipitation Array
#   82,     # Supplemental            Alphanumeric
#           # Precipitation Data
#   83,     # Spare
#   84,     # Velocity Azimuth        Non-geographic Alphanumeric
#           # Display
#   85,     # Cross Section           Raster Image (Reflectivity)
#           # Reflectivity
#   86,     # Cross Section Velocity  Raster Image (Velocity)
#   87,     # Spare
#   88,     # Spare
#   89,     # Layer Composite         Raster Image - Layer 3 Average
#           # Reflectivity
#   90,     # Layer Composite         Raster Image - Layer 3 Maximum
#           # Reflectivity
#   91-92,  # Reserved for internal PUP and RPG Use
#   95,     # Composite Reflectivity  Raster Image/Nongeographic Alpha
#           # Edited for AP
#   96,     # Composite Reflectivity  Raster Image/Nongeographic Alpha
#           # Edited for AP
#   97,     # Composite Reflectivity  Raster Image/Nongeographic Alpha
#           # Edited for AP
#   98,     # Composite Reflectivity  Raster Image/Nongeographic Alpha
#           # Edited for AP
#   100,    # Site Adaptable parameters for VAD Wind Profile (Product 48)
#   101,    # Storm Track             Alphanumeric Block
#   102,    # Hail Index              Alphanumeric Block
#   103,    # Spare
#   104,    # TVS                     Alphanumeric Block
#   105,    # Site Adaptable Parameters for Combined Shear
#   106,    # Spare
#   107,    # Surface Rainfall (1 hr) Alphanumeric Block
#   108,    # Surface Rainfall (3 hr) Alphanumeric Block
#   109,    # Storm Total Rainfall    Alphanumeric Block
#           # Accumulation
#   110,    # Clutter Likelihood      Alphanumeric Block
#           # Reflectivity
#   111,    # Clutter Likelihood      Alphanumeric Block
#           # Doppler
#   112-131,# Reserved for Future Products
#   136,    # SuperOb Adaptable       Latitude, Longitude
#           #                         (ICD packet code 27)
#   139,    # Spare
#   140,    # Gust Front MIGFA        Generic Data Format
#   141,    # Mesocyclone Detection   Geographic and Non-geographic Alpha
#   143,    # Tornado Vortex          Geographic and Non-geographic Alpha
#           # Signature Rapid Update
#   149,    # Digital Mesocyclone     Generic Data Format
#           # Detection
#   152,    # Archive III Status      Product Generic Data Format
#   156,    # Eddy Dissipation Rate   Digital Radial Data Array
#   157,    # Eddy Dissipation Rate   Digital Radial Data Array
#           # Confidence
#   166,    # Melting Layer           Linked Contour Vectors/
#           #                         Set Color Level
#   176,    # Digital Instantaneous   Generic Radial Product Format
#           # Precipitation Rate
#   178-193,# Reserved for Future Products
#   196-198,# Reserved for Future Products
#   200-210,# Reserved for Future Products
#   211-220,# Reserved for Future Products
#   221-230,# Reserved for Future Products
#   231-240,# Reserved for Future Products
#   241-250,# Reserved for Future Products
#   251-260,# Reserved for Future Products
#   261-270,# Reserved for Future Products
#   271-280,# Reserved for Future Products
#   281-290,# Reserved for Future Products
#   291-296,# Reserved for Internal RPG Use.
#   297-299,# Reserved for Internal RPG Use.
