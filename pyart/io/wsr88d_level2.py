"""
"""

import struct
from collections import OrderedDict
import os
import numpy as np
from datetime import datetime, timedelta


class Archive2File():
    """
    Class for reading NEXRAD (WSR-88D) Archive II files.

    Parameters
    ----------
    filename : str
        Filename of Archive II file to read.

    """
    def __init__(self, filename):
        """ initalize the object. """

        # read in the volume header and compression_record
        fh = open(filename, 'rb')
        self.volume_header = _unpack_from_file(fh, VOLUME_HEADER)
        self._compression_record = fh.read(12)
        # TODO check if compressed

        # read records
        self._offsets = []
        self._records = []
        while fh.tell() != os.fstat(fh.fileno()).st_size:
            self._offsets.append(fh.tell())
            self._records.append(_get_record(fh))
        fh.close()

        # pull out records with data
        self.msg31s = [r for r in self._records if r['header']['type'] == 31]
        elev_nums = np.array([m['msg_31_header']['elevation_number']
                             for m in self.msg31s])

        self.scan_msgs = [np.where(elev_nums == i + 1)[0]
                          for i in range(elev_nums.max())]
        self.nscans = len(self.scan_msgs)

    def get_scan_moment(self, scan_num, moment):
        """
        A
        """
        raw_data = self._get_scan_moment_raw(scan_num, moment)

        if moment == 'REF':
            offset = 66.0
            scale = 2.0
        elif moment == 'VEL':
            offset = 129.0
            scale = 1.0         # XXX sometime 2
        elif moment == 'SW':
            offset = 129.0
            scale = 2.0
        elif moment == 'ZDR':
            offset = 128.0
            scale = 16.0
        elif moment == 'PHI':
            offset = 2.0
            scale = 2.8361
        elif moment == 'RHO':
            offset = -60.5
            scale = 300.0
        else:
            offset = 0.0        # these value should never be used
            scale = 1.0
        return (raw_data - offset) / (scale)

    def _get_scan_moment_raw(self, scan_num, moment):
        """
        Retrieve
        """
        msg_for_elev = self.scan_msgs[scan_num]
        nrays = len(msg_for_elev)
        ngates = self.msg31s[msg_for_elev[0]][moment]['ngates']

        data = np.empty((nrays, ngates), dtype='>i1')
        for i, msg_num in enumerate(msg_for_elev):
            data[i] = self.msg31s[msg_num][moment]['data']
        return np.ma.masked_less_equal(data, 2)

    def get_scan_time(self, scan_num):
        """ A """
        days = self._array_from_msg31_headers(scan_num, 'collect_date')
        secs = self._array_from_msg31_headers(scan_num, 'collect_ms') / 1000.

        offset = timedelta(days=int(days[0]) - 1, seconds=int(secs[0]))
        time_start = datetime(1970, 1, 1) + offset

        time = secs - secs[0] + (days - days[0]) * 86400
        return time_start, time

    def get_scan_moment_names(self, scan_num):
        """ A """
        dic = self.msg31s[self.scan_msgs[scan_num][0]]
        moments = ['REF', 'VEL', 'SW', 'ZDR', 'PHI', 'RHO']
        return [k for k in dic.keys() if k in moments]

    def get_scan_range(self, scan_num, moment):
        """ Return an array of gate ranges for a given scan and moment. """
        dic = self.msg31s[self.scan_msgs[scan_num][0]][moment]
        ngates = dic['ngates']
        first_gate = dic['first_gate']
        gate_spacing = dic['gate_spacing']
        return np.arange(ngates) * gate_spacing + first_gate

    def get_location(self):
        """ A """
        dic = self.msg31s[0]['VOL']
        return dic['lat'], dic['lon'], dic['height']

    def get_scan_azimuth_angles(self, scan_num):
        """ A """
        return self._array_from_msg31_headers(scan_num, 'azimuth_angle')

    def get_scan_elevation_angles(self, scan_num):
        """ A """
        return self._array_from_msg31_headers(scan_num, 'elevation_angle')

    def _array_from_msg31_headers(self, scan_num, key):
        """ A """
        return np.array([self.msg31s[i]['msg_31_header'][key]
                         for i in self.scan_msgs[scan_num]])


RECORD_SIZE = 2432
MSG_HEADER_SIZE = 16


def _get_record(fh):
    """ Retrieve and unpack a record from a file. """
    dic = {'header': _unpack_from_file(fh, MSG_HEADER)}
    msg_type = dic['header']['type']
    if msg_type == 31:
        buf = fh.read(dic['header']['size'] * 2 - 4)
        msg_31_header = _unpack_structure(buf[:68], MSG_31)

        block_pointers = [v for k, v in msg_31_header.iteritems()
                          if k.startswith('block_pointer') and v > 0]
        for block_pointer in block_pointers:
            block_name, block_dic = _get_msg31_data_block(buf, block_pointer)
            dic[block_name] = block_dic

        dic['msg_31_header'] = msg_31_header
        #dic['raw_message'] = buf

    else:   # not message 31 or 1
        fh.read(RECORD_SIZE - MSG_HEADER_SIZE)

    return dic


def _get_msg31_data_block(buf, ptr):
    """ unpack a msg_31 data block into a dictionary. """
    block_type = buf[ptr]
    block_name = buf[ptr + 1: ptr + 4].strip()

    if block_name == 'VOL':
        dic = _unpack_structure(buf[ptr:ptr + 44], VOLUME_DATA_BLOCK)
    elif block_name == 'ELV':
        dic = _unpack_structure(buf[ptr:ptr + 12], ELEVATION_DATA_BLOCK)
    elif block_name == 'RAD':
        dic = _unpack_structure(buf[ptr:ptr + 20], RADIAL_DATA_BLOCK)
    elif block_name in ['REF', 'VEL', 'SW', 'ZDR', 'PHI', 'RHO']:
        dic = _unpack_structure(buf[ptr:ptr + 28], GENERIC_DATA_BLOCK)
        ngates = dic['ngates']
        if block_name == 'PHI':
            data = np.fromstring(buf[ptr + 28: ptr + 28 + ngates * 2], '>i2')
        else:
            data = np.fromstring(buf[ptr + 28: ptr + 28 + ngates], '>i1')
        dic['data'] = data
    else:
        dic = {}
    return block_name, dic


####

def _unpack_from_file(f, structure_dict, endiness='>'):
    """ Unpack a structure after reading from a file. """
    fmt = endiness + ''.join(structure_dict.values())
    l = struct.unpack(fmt, f.read(struct.calcsize(fmt)))
    return dict(zip(structure_dict, l))


def _unpack_structure(string, structure_dict):
    """ Unpack a structure """
    fmt = '>' + ''.join(structure_dict.values())
    l = struct.unpack(fmt, string)
    return dict(zip(structure_dict, l))


# Data formats
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

#
#
VOLUME_HEADER = OrderedDict([
    ('tape', '9s'),
    ('extension', '3s'),
    ('date', 'I'),
    ('time', 'I'),
    ('icao', '4s')
])


# rememeber 4-byte control word before this at LDM record start...
# then 8-byes padding???? (12 bytes total)
# Table II Message Header Data
# page 3-8
MSG_HEADER = MESSAGE_HEADER = OrderedDict([
    ('size', INT2),                 # size of data, no including header
    ('channels', INT1),
    ('type', INT1),
    ('seq_id', INT2),
    ('date', INT2),
    ('ms', INT4),
    ('segments', INT2),
    ('seg_num', INT2),
])

# Table XVII Digital Radar Generic Format Blocks (Message Type 31)
# pages 3-93 to 3-95
MSG_31 = MESSAGE_31 = OrderedDict([
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
])

# Table XVII-B Data Block (Descriptor of Generic Data Moment Type)
# pages 3-96 and 3-97
GENERIC_DATA_BLOCK = OrderedDict([
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
])

# Table XVII-E Data Block (Volume Data Constant Type)
# page 3-98
VOLUME_DATA_BLOCK = OrderedDict([
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
])

# Table XVII-F Data Block (Elevation Data Constant Type)
# page 3-99
ELEVATION_DATA_BLOCK = OrderedDict([
    ('block_type', '1s'),
    ('data_name', '3s'),
    ('lrtup', INT2),
    ('atmos', SINT2),
    ('refl_calib', REAL4),
])

# Table XVII-H Data Block (Radial Data Constant Type)
# pages 3-99
RADIAL_DATA_BLOCK = OrderedDict([
    ('block_type', '1s'),
    ('data_name', '3s'),
    ('lrtup', INT2),
    ('unambig_range', SINT2),
    ('noise_h', REAL4),
    ('noise_v', REAL4),
    ('nyquist_vel', SINT2),
    ('spare', '2s')
])
