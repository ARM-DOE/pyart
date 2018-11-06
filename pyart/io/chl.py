"""
pyart.io.chl
============

Utilities for reading CSU-CHILL CHL files.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    ChlFile

.. autosummary::
    :toctree: generated/

    read_chl
    _unpack_structure

"""

import struct
from datetime import datetime

import numpy as np

from ..config import FileMetadata, get_fillvalue
from ..core.radar import Radar
from .common import make_time_unit_str, _test_arguments, prepare_for_read


def read_chl(filename, field_names=None, additional_metadata=None,
             file_field_names=None, exclude_fields=None,
             include_fields=None, use_file_field_attributes=True, **kwargs):
    """
    Read a CSU-CHILL CHL file.

    Parameters
    ----------
    filename : str
        Name of CHL file.
    field_names : dict, optional
        Dictionary mapping CHL field names to radar field names. If a
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
        True to use the CHL field names for the field names in the radar
        object. If this case the field_names parameter is ignored.
        The field dictionary will likely only have a 'data' key, unless
        the fields are defined in `additional_metadata`.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters.
    include_fields : list or None, optional
        List of fields to include from the radar object. This is applied
        after the `field_file_names` and `field_names` parameters. Set to
        None to include all fields not in exclude_fields.
    use_file_field_attributes : bool, optional
        True to use information provided by in the file to set the field
        attribute `long_name`, `units`, `valid_max`, and `valid_min`.  False
        will not set these unless they are defined in the configuration file
        or in `additional_metadata`.

    Returns
    -------
    radar : Radar
        Radar object containing data from CHL file.

    """
    # test for non empty kwargs
    _test_arguments(kwargs)

    # create metadata retrival object
    filemetadata = FileMetadata('chl', field_names, additional_metadata,
                                file_field_names, exclude_fields,
                                include_fields)

    # read data
    chl_file = ChlFile(prepare_for_read(filename))

    # time
    time = filemetadata('time')
    tdata = np.array(chl_file.time)
    min_time = np.floor(tdata.min())
    time['data'] = (tdata - min_time).astype('float64')
    time['units'] = make_time_unit_str(datetime.utcfromtimestamp(min_time))

    # range
    _range = filemetadata('range')
    _range['data'] = (np.array(range(chl_file.ngates)) *
                      chl_file.gate_spacing + chl_file.first_gate_offset)
    _range['meters_between_gates'] = np.array(chl_file.gate_spacing)
    _range['meters_to_center_of_first_gate'] = 0.0

    # scan_type
    scan_type = SCAN_MODE_NAMES[chl_file.scan_types[-1]]

    # fields
    fields = {}
    for i, fdata in chl_file.fields.items():

        field_info = chl_file.field_info[i]
        field_name = filemetadata.get_field_name(field_info['name'])
        if field_name is None:
            continue

        field_dic = filemetadata(field_name)
        np.ma.set_fill_value(fdata, get_fillvalue())
        field_dic['data'] = fdata
        field_dic['_FillValue'] = get_fillvalue()

        if use_file_field_attributes:
            field_dic['long_name'] = field_info['descr']
            field_dic['units'] = field_info['units']
            field_dic['valid_max'] = field_info['max_val']
            field_dic['valid_min'] = field_info['min_val']

        fields[field_name] = field_dic

    # metadata
    metadata = filemetadata('metadata')
    metadata['instrument_name'] = chl_file.radar_info['radar_name']
    metadata['original_container'] = 'CHL'

    # longitude, latitude, altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')
    latitude['data'] = np.array([chl_file.radar_info['latitude']], 'f8')
    longitude['data'] = np.array([chl_file.radar_info['longitude']], 'f8')
    altitude['data'] = np.array([chl_file.radar_info['altitude']], 'f8')

    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
    # sweep_end_ray_index
    sweep_number = filemetadata('sweep_number')
    sweep_mode = filemetadata('sweep_mode')
    fixed_angle = filemetadata('fixed_angle')
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')
    sweep_number['data'] = np.arange(chl_file.num_sweeps, dtype='int32')
    sweep_mode['data'] = np.array(
        [SCAN_MODE_NAMES[i] for i in chl_file.scan_types], dtype='S')
    fixed_angle['data'] = np.array(chl_file.fixed_angle, dtype='float32')

    ray_count = chl_file.rays_per_sweep
    ssri = np.cumsum(np.append([0], ray_count[:-1])).astype('int32')
    sweep_start_ray_index['data'] = ssri
    sweep_end_ray_index['data'] = np.cumsum(ray_count).astype('int32') - 1

    # azimuth, elevation
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')
    azimuth['data'] = np.array(chl_file.azimuth, dtype='float32')
    elevation['data'] = np.array(chl_file.elevation, dtype='float32')

    # instrument parameters
    instrument_parameters = None

    chl_file.close()
    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)


class ChlFile(object):
    """
    A file object for CHL data.

    Parameters
    ----------
    filename : str or file-like.
        Name of CHL file to read or a file-like object pointing to the
        beginning of such a file.
    ns_time : bool
        True to determine ray collection times to the nano-second, False
        will only determine times to the second.
    debug : bool
        True to keep packet data in the _packets attribute to aid in
        debugging.

    Attributes
    ----------
    ngates : int
        Number of gates per ray.
    num_sweeps : int
        Number of sweeps in the volume.
    gate_spacing : float
        Spacing in meters between gates.
    first_gate_offset : float
        Distance in meters to the first range gate.
    time : list of ints
        Time in seconds in epoch for each ray in the volume.
    azimuth : list of floats
        Azimuth angle for each ray in the volume in degrees.
    elevation : list of floats
        Elevation angle for each ray in the volume in degrees.
    fixed_angle : list of floats
        Fixed angles for each sweep.
    sweep_number : list of ints
        Sweep numbers reported in file.
    scan_types : list of ints
        Chill defined scan type for each sweep.
    rays_per_sweep : list of ints
        Number of rays in each sweep.
    fields : dict
        Dictionary of field data index by field number.
    radar_info : dict
        Radar information recorded in the file.
    field_info : dict
        Field information (limits, name, etc.) recorded in the file.
    processor_info : dict
        Porcessor information recorded in the file.

    """

    def __init__(self, filename, ns_time=True, debug=False):

        # initalize attributes
        self.ngates = None
        self.num_sweeps = None
        self.gate_spacing = None
        self.time = []
        self.azimuth = []
        self.elevation = []
        self.fixed_angle = []
        self.sweep_number = []
        self.scan_types = []
        self.rays_per_sweep = None
        self.fields = {}
        self.radar_info = None
        self.field_info = {}
        self.processor_info = None
        self.first_gate_offset = None

        # private attributes
        self._dstring = b''     # string containing field data.
        self._bit_mask = None   # bit mask specifying fields present in file.
        self._dtype = None      # NumPy dtype for a single gate (all fields).
        self._ray_bsize = None  # size in bytes of a single ray (all fields).
        self._packets = []      # List of packets, not set if debug is False.
        self._field_nums = None     # Field number available in file.
        self._rays_in_current_sweep = None  # accumulator for counting rays.
        self._fh = None         # file handler
        self._include_ns_time = ns_time

        # read all blocks from the file
        if hasattr(filename, 'read'):
            self._fh = filename
        else:
            self._fh = open(filename, "rb")
        packet = 1
        while packet is not None:
            packet = self._read_block()
            if debug:
                self._packets.append(packet)

        self._extract_fields()
        self.rays_per_sweep.append(self._rays_in_current_sweep)

    def close(self):
        """ Close the file. """
        self._fh.close()

    def _read_block(self):
        """ Read a block from an open CHL file """
        pld = self._fh.read(8)
        if pld == b'':
            return None
        block_id, length = struct.unpack("<2i", pld)
        payload = self._fh.read(length - 8)

        # parse the block
        if block_id == ARCH_ID_FILE_HDR:
            packet = self._parse_file_hdr_block(payload)
        elif block_id == ARCH_ID_FIELD_SCALE:
            packet = self._parse_field_scale_block(payload)
        elif block_id == ARCH_ID_RAY_HDR:
            packet = self._parse_ray_hdr_block(payload)
        elif block_id == HSK_ID_RADAR_INFO:
            packet = self._parse_radar_info_block(payload)
        elif block_id == HSK_ID_PROCESSOR_INFO:
            packet = self._parse_processor_info_block(payload)
        elif block_id == HSK_ID_SCAN_SEG:
            packet = self._parse_scan_seg_block(payload)
        elif block_id == ARCH_ID_SWEEP_BLOCK:
            packet = self._parse_sweep_block(payload)
        else:
            packet = {}
        packet['block_id'] = block_id
        packet['length'] = length
        return packet

    # Block parsers
    def _parse_file_hdr_block(self, payload):
        """ Parse a field_hdr block. """
        return _unpack_structure(payload, ARCH_FILE_HDR_T)

    def _parse_field_scale_block(self, payload):
        """ Parse a field_scale block. Add scale to field_info attr. """
        packet = _unpack_structure(payload, FIELD_SCALE_T)
        packet['name'] = packet['name'].decode('utf-8').rstrip('\x00')
        packet['units'] = packet['units'].decode('utf-8').rstrip('\x00')
        packet['descr'] = packet['descr'].decode('utf-8').rstrip('\x00')
        self.field_info[packet['bit_mask_pos']] = packet
        return packet

    def _parse_radar_info_block(self, payload):
        """ Parse a radar_info block. Update metadata attribute. """
        packet = _unpack_structure(payload, RADAR_INFO_T)
        packet['radar_name'] = (
            packet['radar_name'].decode('utf-8').rstrip('\x00'))
        self.radar_info = packet.copy()
        return packet

    def _parse_processor_info_block(self, payload):
        """ Parse a processor_info block.  Set dr attribute. """
        packet = _unpack_structure(payload, PROCESSOR_INFO)
        self.gate_spacing = packet['gate_spacing']
        self.first_gate_offset = packet['range_offset']
        self.processor_info = packet.copy()
        return packet

    def _parse_scan_seg_block(self, payload):
        """ Parse a scan_seg_block.  Update sweep attributes. """
        packet = _unpack_structure(payload, SCAN_SEG)
        self.sweep_number.append(packet['sweep_num'])
        self.fixed_angle.append(packet['current_fixed_angle'])
        self.scan_types.append(packet['scan_type'])
        if self.rays_per_sweep is None:
            self.rays_per_sweep = []    # first sweep block
            self._rays_in_current_sweep = 0
        else:
            self.rays_per_sweep.append(self._rays_in_current_sweep)
            self._rays_in_current_sweep = 0
        return packet

    def _parse_sweep_block(self, payload):
        """ Parse a sweep block. Set num_sweeps attribute. """
        packet = {}
        packet['num_sweeps'] = struct.unpack('I', payload[0:4])[0]
        packet['swp_offsets'] = struct.unpack(
            str(packet['num_sweeps']) + 'Q', payload[4:])
        self.num_sweeps = packet['num_sweeps']
        return packet

    def _parse_ray_hdr_block(self, payload):
        """ Parse a ray_hdr block. Update associated attributes. """
        packet = _unpack_structure(payload, ARCH_RAY_HEADER)

        if self._bit_mask is None:
            # this is the first ray_hdr block read
            self.ngates = packet['gates']
            self._bit_mask = packet['bit_mask']
            self._field_nums = [b for b in range(38) if self._bit_mask & 2**b]
            self._dtype = ','.join([DATA_FORMAT[self.field_info[i]['format']]
                                    for i in self._field_nums])
            self._ray_bsize = np.dtype(self._dtype).itemsize * packet['gates']
        else:
            # check that the bit_mask and number of gates are constant
            if packet['bit_mask'] != self._bit_mask:
                raise NotImplementedError('bit_mask is not consistent.')
            if packet['gates'] != self.ngates:
                raise NotImplementedError('number of gates vary.')

        # store ray data and pointing data
        self._dstring += self._fh.read(self._ray_bsize)
        if self._include_ns_time:
            self.time.append(packet['time'] + packet['ns_time']/1e9)
        else:
            self.time.append(packet['time'])
        self.azimuth.append(packet['azimuth'])
        self.elevation.append(packet['elevation'])
        self._rays_in_current_sweep += 1

        return packet

    def _extract_fields(self):
        """ Extract field data from _dstring attribute post read. """
        all_data = np.frombuffer(self._dstring, dtype=self._dtype)
        all_data = all_data.reshape(-1, self.ngates)
        for i, field_num in enumerate(self._field_nums):

            fdata = np.ma.masked_values(all_data[all_data.dtype.names[i]], 0)
            # apply scale and offset factors to interger data types
            if issubclass(fdata.dtype.type, np.integer):
                dat_factor = float(self.field_info[field_num]['dat_factor'])
                dat_bias = float(self.field_info[field_num]['dat_bias'])
                fld_factor = float(self.field_info[field_num]['fld_factor'])
                fdata = (fdata * dat_factor + dat_bias) / fld_factor

            self.fields[field_num] = fdata
        return

# CHL packet types
ARCH_FORMAT_VERSION = 0x00010000
ARCH_ID_CONTROL = 0x5aa80001
ARCH_ID_FIELD_SCALE = 0x5aa80002
ARCH_ID_RAY_HDR = 0x5aa80003
ARCH_ID_FILE_HDR = 0x5aa80004
ARCH_ID_SWEEP_BLOCK = 0x5aa80005
HSK_ID_PROCESSOR_INFO = 0x5aa50003
HSK_ID_RADAR_INFO = 0x5aa50001
HSK_ID_SCAN_SEG = 0x5aa50002

# Additional constants
SCAN_MODE_NAMES = ['ppi', 'rhi', 'fixed', 'manual ppi', 'manual rhi', 'idle']
DATA_FORMAT = ['uint8', 'uint64', 'float32', 'uint16']

##############
# Structures #
##############


def _unpack_structure(string, structure):
    """ Unpack a structure """
    fmt = ''.join([i[1] for i in structure])
    tpl = struct.unpack(fmt, string)
    return dict(zip([i[0] for i in structure], tpl))

ARCH_FILE_HDR_T = (
    ('version', 'I'),
    ('creation_version', 'I'),
    ('creator_id', '32s'),
    ('sweep_table_offset', 'Q'),
)

ARCH_RAY_HEADER = (
    ('azimuth', 'f'),
    ('elevation', 'f'),
    ('azimuth_width', 'f'),
    ('elevation_width', 'f'),
    ('gates', 'H'),
    ('beam_index', 'H'),
    ('ns_time', 'I'),
    ('time', 'Q'),
    ('bit_mask', 'Q'),
    ('ray_number', 'I'),
    ('num_pulses', 'I'),
)

FIELD_SCALE_T = (
    ('format', 'i'),
    ('min_val', 'f'),
    ('max_val', 'f'),
    ('bit_mask_pos', 'i'),
    ('type_hint', 'i'),
    ('fld_factor', 'i'),
    ('dat_factor', 'i'),
    ('dat_bias', 'i'),
    ('name', '32s'),
    ('units', '32s'),
    ('descr', '128s')
)

RADAR_INFO_T = (
    ('radar_name', '32s'),
    ('latitude', 'f'),
    ('longitude', 'f'),
    ('altitude', 'f'),
    ('beamwidth', 'f'),
    ('wavelength_cm', 'f'),
    ('unused1', 'f'),
    ('unused2', 'f'),
    ('unused3', 'f'),
    ('unused4', 'f'),
    ('gain_ant_h', 'f'),
    ('gain_ant_v', 'f'),
    ('zdr_cal_base', 'f'),
    ('phidp_rot', 'f'),
    ('base_radar_constant', 'f'),
    ('unused5', 'f'),
    ('power_measurement_loss_h', 'f'),
    ('power_measurement_loss_v', 'f'),
    ('zdr_cal_base_vhs', 'f'),
    ('test_power_h', 'f'),
    ('test_power_v', 'f'),
    ('dc_loss_h', 'f'),
    ('dc_loss_v', 'f'),
)

PROCESSOR_INFO = (
    ('polarization_mode', 'i'),
    ('processing_mode', 'i'),
    ('pulse_type', 'i'),
    ('test_type', 'i'),
    ('integration_cycle_pulses', 'I'),
    ('clutter_filter_number', 'I'),
    ('range_gate_averaging', 'I'),
    ('indexed_beam_width', 'f'),
    ('gate_spacing', 'f'),
    ('prt_usec', 'f'),
    ('range_start', 'f'),
    ('range_stop', 'f'),
    ('max_gate', 'I'),
    ('test_power', 'f'),
    ('unused1', 'f'),
    ('unused2', 'f'),
    ('test_pulse_range', 'f'),
    ('test_pulse_length', 'f'),
    ('prt2', 'f'),
    ('range_offset', 'f'),
)

SCAN_SEG = (
    ('az_manual', 'f'),
    ('el_manual', 'f'),
    ('az_start', 'f'),
    ('el_start', 'f'),
    ('scan_rate', 'f'),
    ('segname', '24s'),
    ('opt', 'i'),
    ('follow_mode', 'i'),
    ('scan_type', 'i'),
    ('scan_flags', 'I'),
    ('volume_num', 'I'),
    ('sweep_num', 'I'),
    ('time_limit', 'I'),
    ('webtilt', 'I'),
    ('left_limit', 'f'),
    ('right_limit', 'f'),
    ('up_limit', 'f'),
    ('down_limit', 'f'),
    ('step', 'f'),
    ('max_sweeps', 'I'),
    ('filter_break_sweep', 'I'),
    ('clutter_filter1', 'I'),
    ('clutter_filter2', 'I'),
    ('project', '16s'),
    ('current_fixed_angle', 'f'),
)
