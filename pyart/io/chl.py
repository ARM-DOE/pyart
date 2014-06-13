"""
pyart.io.chl
============

Utilities for reading CSU-CHILL CHL files.

"""

import struct
import numpy as np

from ..config import FileMetadata, get_fillvalue
from .radar import Radar
from .common import make_time_unit_str
from .common import radar_coords_to_cart

# XXX first gate offset


def read_chl(filename):
    """
    Read a CHL file.

    Parameters
    ----------
    filename : str
        Name of CHL file.

    Returns
    -------
    radar : Radar
        Radar object containing data from CHL file.

    Notes
    -----
    This is still an alpha-level function so use with caution.

    """

    # create metadata retrival object
    filemetadata = FileMetadata('chill')    # XXX additional parameters

    # read data
    chl_file = CHLfile(filename)

    # time
    time = filemetadata('time')
    time['data'] = chl_file.time
    time['units'] = 'seconds since 1970-01-01 00:00 UTC'    # XXX

    # range
    _range = filemetadata('range')
    _range['data'] = np.array(range(chl_file.ngates)) * chl_file.gate_spacing
    _range['meters_between_gates'] = np.array(chl_file.gate_spacing)
    _range['meters_to_center_of_first_gate'] = 0.0      # XXX

    # scan_type
    scan_type = {'data': SCAN_MODE_NAMES[chl_file.scan_types[-1]]}    # XXX

    # fields XXX
    fields = {}
    for i, fdata in chl_file.fields.items():
        field_info = chl_file.field_info[i]
        field_name = CHILL_FIELD_MAPPING[field_info['name']]
        fields[field_name] = {
            'coordinates': 'elevation azimuth range',
            'data': np.ma.masked_array(fdata),
            'long_name': field_info['descr'],
            'standard_name': CHILL_FIELD_STANDARD_NAME[field_info['name']],
            'units': field_info['units'],
            'valid_max': field_info['max_val'],
            'valid_min': field_info['min_val']
        }

    # metadata
    metadata = filemetadata('metadata')
    metadata['instrument_name'] = chl_file.radar_info['radar_name']
    metadata['original_container'] = 'CHL'

    # longitude, latitude, altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')
    latitude['data'] = chl_file.radar_info['latitude']
    longitude['data'] = chl_file.radar_info['longitude']
    altitude['data'] = chl_file.radar_info['altitude']

    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
    # sweep_end_ray_index
    sweep_number = filemetadata('sweep_number')
    sweep_mode = filemetadata('sweep_mode')
    fixed_angle = filemetadata('fixed_angle')
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')
    sweep_number['data'] = {  # XXX
        'data': range(chl_file.num_sweeps), 'long_name': 'Sweep_number',
        'standard_name': 'sweep_number', 'units': 'count'}
    sweep_mode['data'] = SCAN_MODE_NAMES[chl_file.scan_types[-1]]
    fixed_angle['data'] = chl_file.fixed_angle

    ray_count = chl_file.rays_per_sweep
    ssri = np.cumsum(np.append([0], ray_count[:-1])).astype('int32')
    sweep_start_ray_index['data'] = ssri
    # XXX
    sweep_end_ray_index['data'] = np.cumsum(ray_count).astype('int32')  # - 1

    # azimuth, elevation
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')
    azimuth['data'] = chl_file.azimuth
    elevation['data'] = chl_file.elevation

    # instrument parameters
    instrument_parameters = None

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)


class CHLfile(object):
    """
    A file object for CHL data.

    A `CHLFile` object stores metadata and data from a CHL file.  Metadata is
    stored in dictionaries as attributes of the object, field data is
    stored as NumPy ndarrays as attributes with the field name.

    Parameters
    ----------
    filename : str
        Name of CHL file to read.
    """

    def __init__(self, filename, debug=False):

        # public attributes
        self.filename = filename    # filename
        self.ngates = None          # number of gates per ray
        self.num_sweeps = None      # number of sweeps in the volume
        self.gate_spacing = None    # spacing in meters between gates

        self.time = []              # list of time for each ray
        self.azimuth = []           # list of azimuth angle for each ray
        self.elevation = []         # list of elevation angle for each ray

        self.fixed_angle = []       # list of fixed angles for each sweep
        self.sweep_number = []      # sweep numbers reported in file
        self.scan_types = []        # scan type for each sweep
        self.rays_per_sweep = None  # number of rays per sweep

        self.fields = {}            # dictionary of field data

        self.radar_info = None      # dictionary of radar information
        self.field_info = {}        # field information (limits, name, etc.)
        self.processor_info = None

        # private attributes

        self._dstring = ''      # string containing field data
        self._bit_mask = None   # bit_mask
        self._dtype = None      # dtype of a gate
        self._ray_bsize = None  # size in bytes of a single ray (all fields)
        self._packets = []      # list of packets, not set if debug is False.
        self._field_nums = None  # field number available in file
        self._rays_in_current_sweep = None

        # unknown

        # read all blocks from the file
        self._fh = open(self.filename, "rb")
        packet = 1
        while packet is not None:
            packet = self._read_block()
            if debug:
                self._packets.append(packet)
        self._fh.close()

        self._extract_fields()
        self.rays_per_sweep.append(self._rays_in_current_sweep)

    def _read_block(self):
        """ Read a block from an open CHL file """
        pld = self._fh.read(8)
        if pld == '':
            return None
        block_id, length = struct.unpack("<2i", pld)
        payload = self._fh.read(length - 8)

        if block_id is None:  # redundant sanity check
            return None

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
        packet['name'] = packet['name'].rstrip('\x00')
        packet['units'] = packet['units'].rstrip('\x00')
        packet['descr'] = packet['descr'].rstrip('\x00')
        self.field_info[packet['bit_mask_pos']] = packet
        return packet

    def _parse_radar_info_block(self, payload):
        """ Parse a radar_info block. Update metadata attribute. """
        packet = _unpack_structure(payload, RADAR_INFO_T)
        packet['radar_name'] = packet['radar_name'].rstrip('\x00')
        self.radar_info = packet.copy()
        return packet

    def _parse_processor_info_block(self, payload):
        """ Parse a processor_info block.  Set dr attribute. """
        packet = _unpack_structure(payload, PROCESSOR_INFO)
        self.gate_spacing = packet['gate_spacing']
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
        self.time.append(packet['time'])
        self.azimuth.append(packet['azimuth'])
        self.elevation.append(packet['elevation'])
        self._rays_in_current_sweep += 1

        return packet

    def _extract_fields(self):
        """ Extract field data from _dstring attribute post read. """
        all_data = np.fromstring(self._dstring, dtype=self._dtype)
        all_data = all_data.reshape(-1, self.ngates)

        for i, field_num in enumerate(self._field_nums):
            fsl = self.field_info[field_num]
            dat = all_data[all_data.dtype.names[i]].copy()

            # XXX should be masked arrays
            if not issubclass(dat.dtype.type, np.integer):
                dat[dat == 0] = np.nan
            else:
                dat = dat.astype('float')
                dat[dat == 0] = np.nan
                dat = ((dat * fsl['dat_factor'] + fsl['dat_bias']) /
                       float(fsl['fld_factor']))

            self.fields[field_num] = dat
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
    ('gain_ant_h', 'f'),
    ('gain_ant_v', 'f'),
    ('zdr_cal_base', 'f'),
    ('phidp_rot', 'f'),
    ('base_radar_constant', 'f'),
    ('power_measurement_loss_h', 'f'),
    ('power_measurement_loss_v', 'f'),
    ('zdr_cal_base_vhs', 'f'),
    ('test_power_h', 'f'),
    ('test_power_v', 'f'),
    ('dc_loss_h', 'f'),
    ('dc_loss_v', 'f'),
    ('unknown_0', 'f'),
    ('unknown_1', 'f'),
    ('unknown_2', 'f'),
    ('unknown_3', 'f'),
    ('unknown_4', 'f'),
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
    ('test_pulse_range', 'f'),
    ('test_pulse_length', 'f'),
    ('prt2', 'f'),
    ('range_offset', 'f'),
    ('unknown_0', 'f'),
    ('unknown_1', 'f'),
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

# XXX

_FIELD_TABLE = {
    # Chill field name : (Py-ART field name, field long_name attribute)
    'Z': ('DBZ', 'equivalent_reflectivity_factor'),
    'V': ('VEL', 'radial_velocity_of_scatterers_away_from_instrument'),
    'W': ('WIDTH', 'doppler_spectrum_width'),
    'ZDR': ('ZDR', 'log_differential_reflectivity_hv'),
    'LDRH': ('LDRH', 'log_linear_depolarization_ratio_h'),
    'LDRV': ('LDRV', 'log_linear_depolarization_artio_v'),
    '\xce\xa8 DP': ('PHIDP', 'differential_phase_hv'),
    'KDP': ('KDP', 'specific_differential_phase_hv'),
    '\xcf\x81 HV': ('RHOHV', 'cross_correlation_ratio_hv'),
    'NCP': ('NCP', 'normalized_coherent_power'),
    'H Re(lag 1)': ('H Re(lag 1)', 'real_part_of_lag_1_correlation_h'),
    'V Re(lag 2)': ('V Re(lag 2)', 'real_part_of_lag_2_correlation_v'),
    'VAvgQ': ('VAvgQ', 'v_average_quadrature'),
    'V Im(lag 1)': ('V Im(lag 1)', 'imaginary_part_of_v_at_lag_1'),
    'HAvgQ': ('HAvgQ', 'h_average_quadrature'),
    'H Im(lag 2)': ('H Im(lag 2)', 'imaginary_part_lag_2_correlation_h'),
    'V lag 0': ('V lag 0', 'absolute_value_of_lag_0_correlation_v'),
    'H lag 0': ('H lag 0', 'absolute_value_of_lag_0_correlation_h'),
    'H lag 0 cx':
    ('H lag 0 cx', 'absolute_value_of_lag_0_cross_correlation_h'),
    'H Im(lag 1)':
    ('H Im(lag 1)', 'imaginary_part_of_lag_1_correlation_h'),
    'H Re(lag 2)':
    ('H Re(lag 2)', 'real_part_of_lag_2_correlation_h'),
    'V lag 0 cx':
    ('V lag 0 cx', 'absolute_value_of_lag_0_cross_correlation_v'),
    'V Re(lag 1)': ('V Re(lag 1)', 'real_part_of_lag_1_correlation_v'),
    'V Im(lag 2)':
    ('V Im(lag 2)', 'imaginary_part_of_lag_2_correlation_v'),
    'HV lag 0 I':
    ('HV lag 0 I', 'real_part_of_cross_channel_correlation_at_lag_0'),
    'HV lag 0 Q':
    ('HV lag 0 Q',
     'imaginary_part_of_cross_channel_correlation_at_lag_0'),
    'VAvgI': ('VAvgI', 'v_average_inphase'),
    'HAvgI': ('HAvgI', 'h_average_inphase'),
    '\xcf\x81 HCX': ('RHOHCX', 'lag_0_h_co_to_cross_correlation'),
    '\xcf\x81 VCX': ('RHOVCX', 'lag_0_v_co_to_cross_correlation'),
}
CHILL_FIELD_MAPPING = dict((k, v[0]) for k, v in _FIELD_TABLE.items())
CHILL_FIELD_STANDARD_NAME = dict((k, v[1]) for k, v in _FIELD_TABLE.items())
