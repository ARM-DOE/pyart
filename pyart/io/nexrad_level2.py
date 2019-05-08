"""
pyart.io.nexrad_level2
======================

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    NEXRADLevel2File

.. autosummary::
    :toctree: generated/

    _decompress_records
    _get_record_from_buf
    _get_msg31_data_block
    _structure_size
    _unpack_from_buf
    _unpack_structure


"""



# This file is part of the Py-ART, the Python ARM Radar Toolkit
# https://github.com/ARM-DOE/pyart

# Care has been taken to keep this file free from extraneous dependancies
# so that it can be used by other projects with no/minimal modification.

# Please feel free to use this file in other project provided the license
# below is followed. Keeping the above comment lines would also be helpful
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
import warnings

import numpy as np


class NEXRADLevel2File(object):
    """
    Class for accessing data in a NEXRAD (WSR-88D) Level II file.

    NEXRAD Level II files [1]_, also know as NEXRAD Archive Level II or
    WSR-88D Archive level 2, are available from the NOAA National Climate Data
    Center [2]_ as well as on the UCAR THREDDS Data Server [3]_. Files with
    uncompressed messages and compressed messages are supported. This class
    supports reading both "message 31" and "message 1" type files.

    Parameters
    ----------
    filename : str
        Filename of Archive II file to read.

    Attributes
    ----------
    radial_records : list
        Radial (1 or 31) messages in the file.
    nscans : int
        Number of scans in the file.
    scan_msgs : list of arrays
        Each element specifies the indices of the message in the
        radial_records attribute which belong to a given scan.
    volume_header : dict
        Volume header.
    vcp : dict
        VCP information dictionary.
    _records : list
        A list of all records (message) in the file.
    _fh : file-like
        File like object from which data is read.
    _msg_type : '31' or '1':
        Type of radial messages in file.

    References
    ----------
    .. [1] http://www.roc.noaa.gov/WSR88D/Level_II/Level2Info.aspx
    .. [2] http://www.ncdc.noaa.gov/
    .. [3] http://thredds.ucar.edu/thredds/catalog.html

    """

    def __init__(self, filename):
        """ initalize the object. """
        # read in the volume header and compression_record
        if hasattr(filename, 'read'):
            fh = filename
        else:
            fh = open(filename, 'rb')
        size = _structure_size(VOLUME_HEADER)
        self.volume_header = _unpack_structure(fh.read(size), VOLUME_HEADER)
        compression_record = fh.read(COMPRESSION_RECORD_SIZE)

        # read the records in the file, decompressing as needed
        compression_slice = slice(CONTROL_WORD_SIZE, CONTROL_WORD_SIZE + 2)
        compression_or_ctm_info = compression_record[compression_slice]
        if compression_or_ctm_info == b'BZ':
            buf = _decompress_records(fh)
        # The 12-byte compression record previously held the Channel Terminal
        # Manager (CTM) information. Bytes 4 through 6 contain the size of the
        # record (2432) as a big endian unsigned short, which is encoded as
        # b'\t\x80' == struct.pack('>H', 2432).
        # Newer files zero out this section.
        elif compression_or_ctm_info in (b'\x00\x00', b'\t\x80'):
            buf = fh.read()
        else:
            raise IOError('unknown compression record')
        self._fh = fh

        # read the records from the buffer
        self._records = []
        buf_length = len(buf)
        pos = 0
        while pos < buf_length:
            pos, dic = _get_record_from_buf(buf, pos)
            self._records.append(dic)

        # pull out radial records (1 or 31) which contain the moment data.
        self.radial_records = [r for r in self._records
                               if r['header']['type'] == 31]
        self._msg_type = '31'
        if len(self.radial_records) == 0:
            self.radial_records = [r for r in self._records
                                   if r['header']['type'] == 1]
            self._msg_type = '1'
        if len(self.radial_records) == 0:
            raise ValueError('No MSG31 records found, cannot read file')
        elev_nums = np.array([m['msg_header']['elevation_number']
                              for m in self.radial_records])
        self.scan_msgs = [np.where(elev_nums == i + 1)[0]
                          for i in range(elev_nums.max())]
        self.nscans = len(self.scan_msgs)

        # pull out the vcp record
        msg_5 = [r for r in self._records if r['header']['type'] == 5]

        if len(msg_5):
            self.vcp = msg_5[0]
        else:
            # There is no VCP Data.. This is uber dodgy
            warnings.warn("No MSG5 detected. Setting to meaningless data. "
                          "Rethink your life choices and be ready for errors."
                          "Specifically fixed angle data will be missing")

            self.vcp = None
        return

    def close(self):
        """ Close the file. """
        self._fh.close()

    def location(self):
        """
        Find the location of the radar.

        Returns all zeros if location is not available.

        Returns
        -------
        latitude : float
            Latitude of the radar in degrees.
        longitude : float
            Longitude of the radar in degrees.
        height : int
            Height of radar and feedhorn in meters above mean sea level.

        """
        if self._msg_type == '31':
            dic = self.radial_records[0]['VOL']
            height = dic['height'] + dic['feedhorn_height']
            return dic['lat'], dic['lon'], height
        else:
            return 0.0, 0.0, 0.0

    def scan_info(self, scans=None):
        """
        Return a list of dictionaries with scan information.

        Parameters
        ----------
        scans : list ot None
            Scans (0 based) for which ray (radial) azimuth angles will be
            retrieved.  None (the default) will return the angles for all
            scans in the volume.

        Returns
        -------
        scan_info : list, optional
            A list of the scan performed with a dictionary with keys
            'moments', 'ngates', 'nrays', 'first_gate' and 'gate_spacing'
            for each scan.  The 'moments', 'ngates', 'first_gate', and
            'gate_spacing' keys are lists of the NEXRAD moments and gate
            information for that moment collected during the specific scan.
            The 'nrays' key provides the number of radials collected in the
            given scan.

        """
        info = []

        if scans is None:
            scans = range(self.nscans)
        for scan in scans:
            nrays = self.get_nrays(scan)
            if nrays < 2:
                self.nscans -= 1
                continue
            msg31_number = self.scan_msgs[scan][0]
            msg = self.radial_records[msg31_number]
            nexrad_moments = ['REF', 'VEL', 'SW', 'ZDR', 'PHI', 'RHO']
            moments = [f for f in nexrad_moments if f in msg]
            ngates = [msg[f]['ngates'] for f in moments]
            gate_spacing = [msg[f]['gate_spacing'] for f in moments]
            first_gate = [msg[f]['first_gate'] for f in moments]
            info.append({
                'nrays': nrays,
                'ngates': ngates,
                'gate_spacing': gate_spacing,
                'first_gate': first_gate,
                'moments': moments})
        return info

    def get_vcp_pattern(self):
        """
        Return the numerical volume coverage pattern (VCP) or None if unknown.
        """
        if self.vcp is None:
            return None
        else:
            return self.vcp['msg5_header']['pattern_number']

    def get_nrays(self, scan):
        """
        Return the number of rays in a given scan.

        Parameters
        ----------
        scan : int
            Scan of interest (0 based).

        Returns
        -------
        nrays : int
            Number of rays (radials) in the scan.

        """
        return len(self.scan_msgs[scan])

    def get_range(self, scan_num, moment):
        """
        Return an array of gate ranges for a given scan and moment.

        Parameters
        ----------
        scan_num : int
            Scan number (0 based).
        moment : 'REF', 'VEL', 'SW', 'ZDR', 'PHI', or 'RHO'
            Moment of interest.

        Returns
        -------
        range : ndarray
            Range in meters from the antenna to the center of gate (bin).

        """
        dic = self.radial_records[self.scan_msgs[scan_num][0]][moment]
        ngates = dic['ngates']
        first_gate = dic['first_gate']
        gate_spacing = dic['gate_spacing']
        return np.arange(ngates) * gate_spacing + first_gate

    # helper functions for looping over scans
    def _msg_nums(self, scans):
        """ Find the all message number for a list of scans. """
        return np.concatenate([self.scan_msgs[i] for i in scans])

    def _radial_array(self, scans, key):
        """
        Return an array of radial header elements for all rays in scans.
        """
        msg_nums = self._msg_nums(scans)
        temp = [self.radial_records[i]['msg_header'][key] for i in msg_nums]
        return np.array(temp)

    def _radial_sub_array(self, scans, key):
        """
        Return an array of RAD or msg_header elements for all rays in scans.
        """
        msg_nums = self._msg_nums(scans)
        if self._msg_type == '31':
            tmp = [self.radial_records[i]['RAD'][key] for i in msg_nums]
        else:
            tmp = [self.radial_records[i]['msg_header'][key] for i in msg_nums]
        return np.array(tmp)

    def get_times(self, scans=None):
        """
        Retrieve the times at which the rays were collected.

        Parameters
        ----------
        scans : list or None
            Scans (0-based) to retrieve ray (radial) collection times from.
            None (the default) will return the times for all scans in the
            volume.

        Returns
        -------
        time_start : Datetime
            Initial time.
        time : ndarray
            Offset in seconds from the initial time at which the rays
            in the requested scans were collected.

        """
        if scans is None:
            scans = range(self.nscans)
        days = self._radial_array(scans, 'collect_date')
        secs = self._radial_array(scans, 'collect_ms') / 1000.
        offset = timedelta(days=int(days[0]) - 1, seconds=int(secs[0]))
        time_start = datetime(1970, 1, 1) + offset
        time = secs - int(secs[0]) + (days - days[0]) * 86400
        return time_start, time

    def get_azimuth_angles(self, scans=None):
        """
        Retrieve the azimuth angles of all rays in the requested scans.

        Parameters
        ----------
        scans : list ot None
            Scans (0 based) for which ray (radial) azimuth angles will be
            retrieved. None (the default) will return the angles for all
            scans in the volume.

        Returns
        -------
        angles : ndarray
            Azimuth angles in degress for all rays in the requested scans.

        """
        if scans is None:
            scans = range(self.nscans)
        if self._msg_type == '1':
            scale = 180 / (4096 * 8.)
        else:
            scale = 1.
        return self._radial_array(scans, 'azimuth_angle') * scale

    def get_elevation_angles(self, scans=None):
        """
        Retrieve the elevation angles of all rays in the requested scans.

        Parameters
        ----------
        scans : list or None
            Scans (0 based) for which ray (radial) azimuth angles will be
            retrieved. None (the default) will return the angles for
            all scans in the volume.

        Returns
        -------
        angles : ndarray
            Elevation angles in degress for all rays in the requested scans.

        """
        if scans is None:
            scans = range(self.nscans)
        if self._msg_type == '1':
            scale = 180 / (4096 * 8.)
        else:
            scale = 1.
        return self._radial_array(scans, 'elevation_angle') * scale

    def get_target_angles(self, scans=None):
        """
        Retrieve the target elevation angle of the requested scans.

        Parameters
        ----------
        scans : list or None
            Scans (0 based) for which the target elevation angles will be
            retrieved. None (the default) will return the angles for all
            scans in the volume.

        Returns
        -------
        angles : ndarray
            Target elevation angles in degress for the requested scans.

        """
        if scans is None:
            scans = range(self.nscans)
        if self._msg_type == '31':
            if self.vcp is not None:
                cut_parameters = self.vcp['cut_parameters']
            else:
                cut_parameters = [{'elevation_angle': 0.0}] * self.nscans
            scale = 360. / 65536.
            return np.array([cut_parameters[i]['elevation_angle'] * scale
                             for i in scans], dtype='float32')
        else:
            scale = 180 / (4096 * 8.)
            msgs = [self.radial_records[self.scan_msgs[i][0]] for i in scans]
            return np.round(np.array(
                [m['msg_header']['elevation_angle'] * scale for m in msgs],
                dtype='float32'), 1)

    def get_nyquist_vel(self, scans=None):
        """
        Retrieve the Nyquist velocities of the requested scans.

        Parameters
        ----------
        scans : list or None
            Scans (0 based) for which the Nyquist velocities will be
            retrieved. None (the default) will return the velocities for all
            scans in the volume.

        Returns
        -------
        velocities : ndarray
            Nyquist velocities (in m/s) for the requested scans.

        """
        if scans is None:
            scans = range(self.nscans)
        return self._radial_sub_array(scans, 'nyquist_vel') * 0.01

    def get_unambigous_range(self, scans=None):
        """
        Retrieve the unambiguous range of the requested scans.

        Parameters
        ----------
        scans : list or None
            Scans (0 based) for which the unambiguous range will be retrieved.
            None (the default) will return the range for all scans in the
            volume.

        Returns
        -------
        unambiguous_range : ndarray
            Unambiguous range (in meters) for the requested scans.

        """
        if scans is None:
            scans = range(self.nscans)
        # unambiguous range is stored in tenths of km, x100 for meters
        return self._radial_sub_array(scans, 'unambig_range') * 100.

    def get_data(self, moment, max_ngates, scans=None, raw_data=False):
        """
        Retrieve moment data for a given set of scans.

        Masked points indicate that the data was not collected, below
        threshold or is range folded.

        Parameters
        ----------
        moment : 'REF', 'VEL', 'SW', 'ZDR', 'PHI', or 'RHO'
            Moment for which to to retrieve data.
        max_ngates : int
            Maximum number of gates (bins) in any ray.
            requested.
        raw_data : bool
            True to return the raw data, False to perform masking as well as
            applying the appropiate scale and offset to the data.  When
            raw_data is True values of 1 in the data likely indicate that
            the gate was not present in the sweep, in some cases in will
            indicate range folded data.
        scans : list or None.
            Scans to retrieve data from (0 based). None (the default) will
            get the data for all scans in the volume.

        Returns
        -------
        data : ndarray

        """
        if scans is None:
            scans = range(self.nscans)

        # determine the number of rays
        msg_nums = self._msg_nums(scans)
        nrays = len(msg_nums)

        # extract the data
        if moment != 'PHI':
            data = np.ones((nrays, max_ngates), dtype='u1')
        else:
            data = np.ones((nrays, max_ngates), dtype='u2')
        for i, msg_num in enumerate(msg_nums):
            msg = self.radial_records[msg_num]
            if moment not in msg.keys():
                continue
            ngates = min(msg[moment]['ngates'], max_ngates, len(msg[moment]['data']))
            data[i, :ngates] = msg[moment]['data'][:ngates]

        # return raw data if requested
        if raw_data:
            return data

        # mask, scan and offset, assume that the offset and scale
        # are the same in all scans/gates
        for scan in scans:  # find a scan which contains the moment
            msg_num = self.scan_msgs[scan][0]
            msg = self.radial_records[msg_num]
            if moment in msg.keys():
                offset = np.float32(msg[moment]['offset'])
                scale = np.float32(msg[moment]['scale'])
                mask = data <= 1
                scaled_data = (data - offset) / scale
                return np.ma.array(scaled_data, mask=mask)

        # moment is not present in any scan, mask all values
        return np.ma.masked_less_equal(data, 1)


def _decompress_records(file_handler):
    """
    Decompressed the records from an BZ2 compressed Archive 2 file.
    """
    file_handler.seek(0)
    cbuf = file_handler.read()    # read all data from the file
    decompressor = bz2.BZ2Decompressor()
    skip = _structure_size(VOLUME_HEADER) + CONTROL_WORD_SIZE
    buf = decompressor.decompress(cbuf[skip:])
    while len(decompressor.unused_data):
        cbuf = decompressor.unused_data
        decompressor = bz2.BZ2Decompressor()
        buf += decompressor.decompress(cbuf[CONTROL_WORD_SIZE:])

    return buf[COMPRESSION_RECORD_SIZE:]


def _get_record_from_buf(buf, pos):
    """ Retrieve and unpack a NEXRAD record from a buffer. """
    dic = {'header': _unpack_from_buf(buf, pos, MSG_HEADER)}
    msg_type = dic['header']['type']

    if msg_type == 31:
        new_pos = _get_msg31_from_buf(buf, pos, dic)
    elif msg_type == 5:
        # Sometimes we encounter incomplete buffers
        try:
            new_pos = _get_msg5_from_buf(buf, pos, dic)
        except struct.error:
            warnings.warn("Encountered incomplete MSG5. File may be corrupt.",
                          RuntimeWarning)
            new_pos = pos + RECORD_SIZE
    elif msg_type == 29:
        new_pos = _get_msg29_from_buf(buf, pos, dic)
        warnings.warn("Message 29 encountered, not parsing.", RuntimeWarning)
    elif msg_type == 1:
        new_pos = _get_msg1_from_buf(buf, pos, dic)
    else:   # not message 31 or 1, no decoding performed
        new_pos = pos + RECORD_SIZE

    return new_pos, dic

def _get_msg29_from_buf(pos, dic):
    msg_size = dic['header']['size']
    if msg_size == 65535:
        msg_size = dic['header']['segments'] << 16 | dic['header']['seg_num']
    msg_header_size = _structure_size(MSG_HEADER)
    new_pos = pos + msg_header_size + msg_size
    return new_pos

def _get_msg31_from_buf(buf, pos, dic):
    """ Retrieve and unpack a MSG31 record from a buffer. """
    msg_size = dic['header']['size'] * 2 - 4
    msg_header_size = _structure_size(MSG_HEADER)
    new_pos = pos + msg_header_size + msg_size
    mbuf = buf[pos + msg_header_size:new_pos]
    msg_31_header = _unpack_from_buf(mbuf, 0, MSG_31)
    block_pointers = [v for k, v in msg_31_header.items()
                      if k.startswith('block_pointer') and v > 0]
    for block_pointer in block_pointers:
        block_name, block_dic = _get_msg31_data_block(mbuf, block_pointer)
        dic[block_name] = block_dic

    dic['msg_header'] = msg_31_header
    return new_pos


def _get_msg31_data_block(buf, ptr):
    """ Unpack a msg_31 data block into a dictionary. """
    block_name = buf[ptr + 1: ptr + 4].decode('ascii').strip()

    if block_name == 'VOL':
        dic = _unpack_from_buf(buf, ptr, VOLUME_DATA_BLOCK)
    elif block_name == 'ELV':
        dic = _unpack_from_buf(buf, ptr, ELEVATION_DATA_BLOCK)
    elif block_name == 'RAD':
        dic = _unpack_from_buf(buf, ptr, RADIAL_DATA_BLOCK)
    elif block_name in ['REF', 'VEL', 'SW', 'ZDR', 'PHI', 'RHO']:
        dic = _unpack_from_buf(buf, ptr, GENERIC_DATA_BLOCK)
        ngates = dic['ngates']
        ptr2 = ptr + _structure_size(GENERIC_DATA_BLOCK)
        if block_name == 'PHI':
            data = np.frombuffer(buf[ptr2: ptr2 + ngates * 2], '>u2')
        else:
            data = np.frombuffer(buf[ptr2: ptr2 + ngates], '>u1')
        dic['data'] = data
    else:
        dic = {}
    return block_name, dic


def _get_msg1_from_buf(buf, pos, dic):
    """ Retrieve and unpack a MSG1 record from a buffer. """
    msg_header_size = _structure_size(MSG_HEADER)
    msg1_header = _unpack_from_buf(buf, pos + msg_header_size, MSG_1)
    dic['msg_header'] = msg1_header

    sur_nbins = int(msg1_header['sur_nbins'])
    doppler_nbins = int(msg1_header['doppler_nbins'])

    sur_step = int(msg1_header['sur_range_step'])
    doppler_step = int(msg1_header['doppler_range_step'])

    sur_first = int(msg1_header['sur_range_first'])
    doppler_first = int(msg1_header['doppler_range_first'])
    if doppler_first > 2**15:
        doppler_first = doppler_first - 2**16

    if msg1_header['sur_pointer']:
        offset = pos + msg_header_size + msg1_header['sur_pointer']
        data = np.frombuffer(buf[offset:offset+sur_nbins], '>u1')
        dic['REF'] = {
            'ngates': sur_nbins,
            'gate_spacing': sur_step,
            'first_gate': sur_first,
            'data': data,
            'scale': 2.,
            'offset': 66.,
        }
    if msg1_header['vel_pointer']:
        offset = pos + msg_header_size + msg1_header['vel_pointer']
        data = np.frombuffer(buf[offset:offset+doppler_nbins], '>u1')
        dic['VEL'] = {
            'ngates': doppler_nbins,
            'gate_spacing': doppler_step,
            'first_gate': doppler_first,
            'data': data,
            'scale': 2.,
            'offset': 129.0,
        }
        if msg1_header['doppler_resolution'] == 4:
            # 1 m/s resolution velocity, offset remains 129.
            dic['VEL']['scale'] = 1.
    if msg1_header['width_pointer']:
        offset = pos + msg_header_size + msg1_header['width_pointer']
        data = np.frombuffer(buf[offset:offset+doppler_nbins], '>u1')
        dic['SW'] = {
            'ngates': doppler_nbins,
            'gate_spacing': doppler_step,
            'first_gate': doppler_first,
            'data': data,
            'scale': 2.,
            'offset': 129.0,
        }
    return pos + RECORD_SIZE


def _get_msg5_from_buf(buf, pos, dic):
    """ Retrieve and unpack a MSG1 record from a buffer. """
    msg_header_size = _structure_size(MSG_HEADER)
    msg5_header_size = _structure_size(MSG_5)
    msg5_elev_size = _structure_size(MSG_5_ELEV)

    dic['msg5_header'] = _unpack_from_buf(buf, pos + msg_header_size, MSG_5)
    dic['cut_parameters'] = []
    for i in range(dic['msg5_header']['num_cuts']):
        pos2 = pos + msg_header_size + msg5_header_size + msg5_elev_size * i
        dic['cut_parameters'].append(_unpack_from_buf(buf, pos2, MSG_5_ELEV))
    return pos + RECORD_SIZE


def _structure_size(structure):
    """ Find the size of a structure in bytes. """
    return struct.calcsize('>' + ''.join([i[1] for i in structure]))


def _unpack_from_buf(buf, pos, structure):
    """ Unpack a structure from a buffer. """
    size = _structure_size(structure)
    return _unpack_structure(buf[pos:pos + size], structure)


def _unpack_structure(string, structure):
    """ Unpack a structure from a string. """
    fmt = '>' + ''.join([i[1] for i in structure])  # NEXRAD is big-endian
    lst = struct.unpack(fmt, string)
    return dict(zip([i[0] for i in structure], lst))


# NEXRAD Level II file structures and sizes
# The deails on these structures are documented in:
# "Interface Control Document for the Achive II/User" RPG Build 12.0
# Document Number 2620010E
# and
# "Interface Control Document for the RDA/RPG" Open Build 13.0
# Document Number 2620002M
# Tables and page number refer to those in the second document unless
# otherwise noted.
RECORD_SIZE = 2432
COMPRESSION_RECORD_SIZE = 12
CONTROL_WORD_SIZE = 4

# format of structure elements
# section 3.2.1, page 3-2
CODE1 = 'B'
CODE2 = 'H'
INT1 = 'B'
INT2 = 'H'
INT4 = 'I'
REAL4 = 'f'
REAL8 = 'd'
SINT1 = 'b'
SINT2 = 'h'
SINT4 = 'i'

# Figure 1 in Interface Control Document for the Archive II/User
# page 7-2
VOLUME_HEADER = (
    ('tape', '9s'),
    ('extension', '3s'),
    ('date', 'I'),
    ('time', 'I'),
    ('icao', '4s')
)

# Table II Message Header Data
# page 3-7
MSG_HEADER = (
    ('size', INT2),                 # size of data, no including header
    ('channels', INT1),
    ('type', INT1),
    ('seq_id', INT2),
    ('date', INT2),
    ('ms', INT4),
    ('segments', INT2),
    ('seg_num', INT2),
)

# Table XVII Digital Radar Generic Format Blocks (Message Type 31)
# pages 3-87 to 3-89
MSG_31 = (
    ('id', '4s'),                   # 0-3
    ('collect_ms', INT4),           # 4-7
    ('collect_date', INT2),         # 8-9
    ('azimuth_number', INT2),       # 10-11
    ('azimuth_angle', REAL4),       # 12-15
    ('compress_flag', CODE1),       # 16
    ('spare_0', INT1),              # 17
    ('radial_length', INT2),        # 18-19
    ('azimuth_resolution', CODE1),  # 20
    ('radial_spacing', CODE1),      # 21
    ('elevation_number', INT1),     # 22
    ('cut_sector', INT1),           # 23
    ('elevation_angle', REAL4),     # 24-27
    ('radial_blanking', CODE1),     # 28
    ('azimuth_mode', SINT1),        # 29
    ('block_count', INT2),          # 30-31
    ('block_pointer_1', INT4),      # 32-35  Volume Data Constant XVII-E
    ('block_pointer_2', INT4),      # 36-39  Elevation Data Constant XVII-F
    ('block_pointer_3', INT4),      # 40-43  Radial Data Constant XVII-H
    ('block_pointer_4', INT4),      # 44-47  Moment "REF" XVII-{B/I}
    ('block_pointer_5', INT4),      # 48-51  Moment "VEL"
    ('block_pointer_6', INT4),      # 52-55  Moment "SW"
    ('block_pointer_7', INT4),      # 56-59  Moment "ZDR"
    ('block_pointer_8', INT4),      # 60-63  Moment "PHI"
    ('block_pointer_9', INT4),      # 64-67  Moment "RHO"
)


# Table III Digital Radar Data (Message Type 1)
# pages 3-7 to
MSG_1 = (
    ('collect_ms', INT4),           # 0-3
    ('collect_date', INT2),         # 4-5
    ('unambig_range', SINT2),       # 6-7
    ('azimuth_angle', CODE2),       # 8-9
    ('azimuth_number', INT2),       # 10-11
    ('radial_status', CODE2),       # 12-13
    ('elevation_angle', INT2),      # 14-15
    ('elevation_number', INT2),     # 16-17
    ('sur_range_first', CODE2),     # 18-19
    ('doppler_range_first', CODE2),  # 20-21
    ('sur_range_step', CODE2),      # 22-23
    ('doppler_range_step', CODE2),  # 24-25
    ('sur_nbins', INT2),            # 26-27
    ('doppler_nbins', INT2),        # 28-29
    ('cut_sector_num', INT2),       # 30-31
    ('calib_const', REAL4),         # 32-35
    ('sur_pointer', INT2),          # 36-37
    ('vel_pointer', INT2),          # 38-39
    ('width_pointer', INT2),        # 40-41
    ('doppler_resolution', CODE2),  # 42-43
    ('vcp', INT2),                  # 44-45
    ('spare_1', '8s'),              # 46-53
    ('spare_2', '2s'),              # 54-55
    ('spare_3', '2s'),              # 56-57
    ('spare_4', '2s'),              # 58-59
    ('nyquist_vel', SINT2),         # 60-61
    ('atmos_attenuation', SINT2),   # 62-63
    ('threshold', SINT2),           # 64-65
    ('spot_blank_status', INT2),    # 66-67
    ('spare_5', '32s'),             # 68-99
    # 100+  reflectivity, velocity and/or spectral width data, CODE1
)

# Table XI Volume Coverage Pattern Data (Message Type 5 & 7)
# pages 3-51 to 3-54
MSG_5 = (
    ('msg_size', INT2),
    ('pattern_type', CODE2),
    ('pattern_number', INT2),
    ('num_cuts', INT2),
    ('clutter_map_group', INT2),
    ('doppler_vel_res', CODE1),     # 2: 0.5 degrees, 4: 1.0 degrees
    ('pulse_width', CODE1),         # 2: short, 4: long
    ('spare', '10s')                # halfwords 7-11 (10 bytes, 5 halfwords)
)

MSG_5_ELEV = (
    ('elevation_angle', CODE2),  # scaled by 360/65536 for value in degrees.
    ('channel_config', CODE1),
    ('waveform_type', CODE1),
    ('super_resolution', CODE1),
    ('prf_number', INT1),
    ('prf_pulse_count', INT2),
    ('azimuth_rate', CODE2),
    ('ref_thresh', SINT2),
    ('vel_thresh', SINT2),
    ('sw_thresh', SINT2),
    ('zdr_thres', SINT2),
    ('phi_thres', SINT2),
    ('rho_thres', SINT2),
    ('edge_angle_1', CODE2),
    ('dop_prf_num_1', INT2),
    ('dop_prf_pulse_count_1', INT2),
    ('spare_1', '2s'),
    ('edge_angle_2', CODE2),
    ('dop_prf_num_2', INT2),
    ('dop_prf_pulse_count_2', INT2),
    ('spare_2', '2s'),
    ('edge_angle_3', CODE2),
    ('dop_prf_num_3', INT2),
    ('dop_prf_pulse_count_3', INT2),
    ('spare_3', '2s'),
)

# Table XVII-B Data Block (Descriptor of Generic Data Moment Type)
# pages 3-90 and 3-91
GENERIC_DATA_BLOCK = (
    ('block_type', '1s'),
    ('data_name', '3s'),        # VEL, REF, SW, RHO, PHI, ZDR
    ('reserved', INT4),
    ('ngates', INT2),
    ('first_gate', SINT2),
    ('gate_spacing', SINT2),
    ('thresh', SINT2),
    ('snr_thres', SINT2),
    ('flags', CODE1),
    ('word_size', INT1),
    ('scale', REAL4),
    ('offset', REAL4),
    # then data
)

# Table XVII-E Data Block (Volume Data Constant Type)
# page 3-92
VOLUME_DATA_BLOCK = (
    ('block_type', '1s'),
    ('data_name', '3s'),
    ('lrtup', INT2),
    ('version_major', INT1),
    ('version_minor', INT1),
    ('lat', REAL4),
    ('lon', REAL4),
    ('height', SINT2),
    ('feedhorn_height', INT2),
    ('refl_calib', REAL4),
    ('power_h', REAL4),
    ('power_v', REAL4),
    ('diff_refl_calib', REAL4),
    ('init_phase', REAL4),
    ('vcp', INT2),
    ('spare', '2s'),
)

# Table XVII-F Data Block (Elevation Data Constant Type)
# page 3-93
ELEVATION_DATA_BLOCK = (
    ('block_type', '1s'),
    ('data_name', '3s'),
    ('lrtup', INT2),
    ('atmos', SINT2),
    ('refl_calib', REAL4),
)

# Table XVII-H Data Block (Radial Data Constant Type)
# pages 3-93
RADIAL_DATA_BLOCK = (
    ('block_type', '1s'),
    ('data_name', '3s'),
    ('lrtup', INT2),
    ('unambig_range', SINT2),
    ('noise_h', REAL4),
    ('noise_v', REAL4),
    ('nyquist_vel', SINT2),
    ('spare', '2s')
)
