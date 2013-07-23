"""
pyart.io.nexrad

"""

import bz2
import struct
from datetime import datetime, timedelta

import numpy as np


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
        self.volume_header = _unpack_structure(fh.read(24), VOLUME_HEADER)
        compression_record = fh.read(12)

        # read the records in the file, decompressing as needed
        if compression_record[4:6] == 'BZ':
            buf = _decompress_records(fh)
        elif compression_record[4:6] == '\x00\x00':
            buf = fh.read()
        else:
            raise IOError('unknown compression record')

        fh.close()
        buf_length = len(buf)
        pos = 0
        self._records = []

        # read the records from the buffer
        while pos < buf_length:
            pos, dic = _get_record_from_buf(buf, pos)
            self._records.append(dic)

        # pull out records with data
        self.msg31s = [r for r in self._records if r['header']['type'] == 31]
        elev_nums = np.array([m['msg_31_header']['elevation_number']
                             for m in self.msg31s])

        self.scan_msgs = [np.where(elev_nums == i + 1)[0]
                          for i in range(elev_nums.max())]
        self.nscans = len(self.scan_msgs)

    def get_scan_info(self):
        """ A """
        dic = {'REF': {1: {'ngates': [], 'scans': []},
                       2: {'ngates': [], 'scans': []}},
               'VEL': {1: {'ngates': [], 'scans': []},
                       2: {'ngates': [], 'scans': []}},
               'SW': {1: {'ngates': [], 'scans': []},
                      2: {'ngates': [], 'scans': []}},
               'ZDR': {1: {'ngates': [], 'scans': []},
                       2: {'ngates': [], 'scans': []}},
               'PHI': {1: {'ngates': [], 'scans': []},
                       2: {'ngates': [], 'scans': []}},
               'RHO': {1: {'ngates': [], 'scans': []},
                       2: {'ngates': [], 'scans': []}}}
        for scan in range(self.nscans):
            msg31_number = self.scan_msgs[scan][0]
            msg = self.msg31s[msg31_number]
            res = msg['msg_31_header']['azimuth_resolution']
            for moment in dic.keys():
                if moment in msg.keys():
                    dic[moment][res]['scans'].append(scan)
                    dic[moment][res]['ngates'].append(msg[moment]['ngates'])

        return dic

    def get_data(self, scans, moment, max_gates, raw_data=False):
        """ A """
        msg_nums = np.concatenate([self.scan_msgs[i] for i in scans])
        nrays = len(msg_nums)

        if moment != 'PHI':
            data = np.ones((nrays, max_gates), dtype='i1')
        else:
            data = np.ones((nrays, max_gates), dtype='i2')
        for i, msg_num in enumerate(msg_nums):
            msg = self.msg31s[msg_num]
            data[i, :msg[moment]['ngates']] = msg[moment]['data']

        if raw_data:
            return data
        else:   # mask, scale and offset
            msg = self.msg31s[msg_nums[0]]
            offset = msg[moment]['offset']
            scale = msg[moment]['scale']
            return (np.ma.masked_less_equal(data, 2) - offset) / (scale)

    def get_time_scans(self, scans):
        """ A """
        msg_nums = np.concatenate([self.scan_msgs[i] for i in scans])

        days = np.array([self.msg31s[i]['msg_31_header']['collect_date']
                         for i in msg_nums])
        secs = np.array([self.msg31s[i]['msg_31_header']['collect_ms']
                         for i in msg_nums]) / 1000.

        offset = timedelta(days=int(days[0]) - 1, seconds=int(secs[0]))
        time_start = datetime(1970, 1, 1) + offset

        time = secs - int(secs[0]) + (days - days[0]) * 86400
        return time_start, time

    def get_nrays(self, scan):
        """ A """
        return len(self.scan_msgs[scan])

    def get_azimuth_angles_scans(self, scans):
        """ A """
        msg_nums = np.concatenate([self.scan_msgs[i] for i in scans])
        return np.array([self.msg31s[i]['msg_31_header']['azimuth_angle']
                         for i in msg_nums])

    def get_elevation_angles_scans(self, scans):
        """ A """
        msg_nums = np.concatenate([self.scan_msgs[i] for i in scans])
        return np.array([self.msg31s[i]['msg_31_header']['elevation_angle']
                         for i in msg_nums])

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


RECORD_SIZE = 2432
MSG_HEADER_SIZE = 16


def _decompress_records(fh):
    """
    Decompressed the records from an BZ2 compressed Archive 2 file.
    """
    fh.seek(0)
    cbuf = fh.read()    # read all data from the file
    decompressor = bz2.BZ2Decompressor()
    buf = decompressor.decompress(cbuf[28:])    # skip volume head + size
    while len(decompressor.unused_data):
        cbuf = decompressor.unused_data
        decompressor = bz2.BZ2Decompressor()
        buf += decompressor.decompress(cbuf[4:])  # skip size (4-bytes)

    return buf[12:]


def _get_record_from_buf(buf, pos):
    """ Retrieve and unpack a record from a buffer. """

    packed_header = buf[pos:pos + MSG_HEADER_SIZE]
    dic = {'header': _unpack_structure(packed_header, MSG_HEADER)}
    msg_type = dic['header']['type']

    if msg_type == 31:
        msg_size = dic['header']['size'] * 2 - 4
        new_pos = pos + MSG_HEADER_SIZE + msg_size
        mbuf = buf[pos + MSG_HEADER_SIZE:new_pos]
        msg_31_header = _unpack_structure(mbuf[:68], MSG_31)

        block_pointers = [v for k, v in msg_31_header.iteritems()
                          if k.startswith('block_pointer') and v > 0]
        for block_pointer in block_pointers:
            block_name, block_dic = _get_msg31_data_block(mbuf, block_pointer)
            dic[block_name] = block_dic

        dic['msg_31_header'] = msg_31_header

    else:   # not message 31 or 1, no decoding
        new_pos = pos + RECORD_SIZE

    return new_pos, dic


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

def _unpack_structure(string, structure):
    """ Unpack a structure """
    fmt = '>' + ''.join([i[1] for i in structure])
    l = struct.unpack(fmt, string)
    return dict(zip([i[0] for i in structure], l))

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
VOLUME_HEADER = (
    ('tape', '9s'),
    ('extension', '3s'),
    ('date', 'I'),
    ('time', 'I'),
    ('icao', '4s')
)


# rememeber 4-byte control word before this at LDM record start...
# then 8-byes padding???? (12 bytes total)
# Table II Message Header Data
# page 3-8
MSG_HEADER = MESSAGE_HEADER = (
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
# pages 3-93 to 3-95
MSG_31 = MESSAGE_31 = (
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

# Table XVII-B Data Block (Descriptor of Generic Data Moment Type)
# pages 3-96 and 3-97
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
# page 3-98
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
# page 3-99
ELEVATION_DATA_BLOCK = (
    ('block_type', '1s'),
    ('data_name', '3s'),
    ('lrtup', INT2),
    ('atmos', SINT2),
    ('refl_calib', REAL4),
)

# Table XVII-H Data Block (Radial Data Constant Type)
# pages 3-99
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
