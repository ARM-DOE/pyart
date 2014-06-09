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
        Radar object containing data from MDV file.

    Notes
    -----
    This is still an alpha-level function so use with caution.

    """

    chl_file = CHLfile(filename)
    return chl_file.return_pyart_radar()

    

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

    def __init__(self, filename):
        self.filename = filename

        self.field_scale_list = {}
        self.num_sweeps=0
        self.radar_info = []
        self.processor_info = []
        self.azimuth = []
        self.elevation = []
        self.fixed_angle = []
        self.sweep_end = []
        self.data_array = []
        self._pyart_fields = {}
        self.time = []
        self.fields = {}
        self._current_ray_num = 0
        self.metadata = {}

        self._chl_arch_open_archive()
        self._chl_close_archive()

        self.sweep_end.append(self._current_ray_num)
        del(self.sweep_end[0])
        self._process_data_blocks()
        self._range = np.array(range(self.num_gates)) * self.dr

        self.sweep_start = [0, ]
        self.sweep_start.extend([idx + 1 for idx in self.sweep_end[0:-1]])

    def return_pyart_radar(self):
        self._range_dict['data'] = self._range
        self._range_dict['meters_between_gates'] = np.array(self.dr),

        self._time_dict['data'] = self.time
        self._sweep_number = {
            'data': range(self.sweep_num), 'long_name': 'Sweep_number',
            'standard_name': 'sweep_number', 'units': 'count'}

        # field processing, we place them into pyart style field dicts
        for i in range(0, len(self.field_scale_list)):
            if(self.field_scale_list[i]['name'] in self.fields.keys()):
                self._pyart_fields[self._variable_name_lookup[self.field_scale_list[i]['name']][0]] = self._field_to_pyart_field(
                    self.fields[self.field_scale_list[i]['name']], self.field_scale_list[i])

        return Radar(
            self._time_dict, self._range_dict, self._pyart_fields, self.metadata, {'data':
                                                                                   self.scan_mode},
            {'data': self._radar_info['latitude']}, {'data': self._radar_info[
                'longitude']}, {'data': self._radar_info['altitude']},
            {'data': self._sweep_number}, {'data': self.scan_mode}, self._pack_fixed_angle(
                self.fixed_angle), {'data': self.sweep_start},
            {'data': self.sweep_end}, self._pack_azimuth(self.azimuth), self._pack_elevation(self.elevation))

    def _chl_arch_open_archive(self):
        self.f = open(self.filename, "rb")
        packet = 1
        while packet is not None:
            packet = self._chl_arch_read_block()

    def _chl_close_archive(self):
        self.f.close()

    # I need to extract this to several shorter funcs
    def _chl_arch_read_block(self):
        pld = self.f.read(8)
        if(pld == ''):
            return None
        id, length = struct.unpack("<2i", pld)
        payload = self.f.read(length - 8)

        if id is None:  # redundant sanity check
            return None

        packet = {}
        if hex(id) == '0x5aa80004':  # arch_file_hdr_t
            packet = dict(
                zip(self.arch_file_hdr_t, struct.unpack(self.arch_file_hdr_fstring, payload)))

        elif hex(id) == '0x5aa80002':  # field_scale_t
            packet = dict(
                zip(self.field_scale_t, struct.unpack(self.field_scale_fstring, payload)))
            self._parse_field_struct_packet(packet)

        elif hex(id) == '0x5aa80003':  # arch_ray_header
            packet = dict(
                zip(self.arch_ray_header, struct.unpack(self.arch_ray_fstring, payload)))
            self._current_ray_num = packet['ray_number']
            fmat_string = self._format_string_from_bitmask(packet['bit_mask'])
            data_packet = self.f.read(
                struct.calcsize(fmat_string) * packet['gates'])
            self.data_array.append(
                struct.unpack(fmat_string * packet['gates'], data_packet))
            # We should only do this once, but this will work for now.
            self.field_count = self._num_fields_from_bitmask(
                packet['bit_mask'])
            self.bit_mask = packet['bit_mask']
            self.time.append(packet['time'])
            self.azimuth.append(packet['azimuth'])
            self.elevation.append(packet['elevation'])
            # This assumes number of gates is constant. Not the best assumption
            # however.
            self.num_gates = packet['gates']

        elif hex(id) == '0x5aa50001':  # radar_info_t
            packet = dict(
                zip(self.radar_info_t, struct.unpack(self.radar_info_fstring, payload)))
            self._radar_info = packet.copy()
            self.metadata['instrument_name'] = packet[
                'radar_name'].rstrip('\x00')
            self.metadata['original_container'] = 'CHL'

        elif hex(id) == '0x5aa50003':  # processor_info
            packet = dict(
                zip(self.processor_info_t, struct.unpack(self.processor_info_fstring, payload)))
            self.dr = packet['gate_spacing']

        elif hex(id) == '0x5aa50002':  # scan_seg
            packet = dict(
                zip(self.scan_seg, struct.unpack(self.scan_seg_fstring, payload)))
            self.sweep_num = packet['sweep_num']
            self.sweep_end.append(self._current_ray_num)
            self.fixed_angle.append(packet['current_fixed_angle'])
            self.scan_mode = self._scan_mode_names[packet['scan_type']]

        elif hex(id) == '0x5aa80005':
            packet['num_sweeps'] = struct.unpack('I', payload[0:4])[0]
            packet['swp_offsets'] = struct.unpack(
                str(packet['num_sweeps']) + 'Q', payload[4:])
            self.num_sweeps = packet['num_sweeps']

        packet['id'] = hex(id)
        packet['length'] = length
        return packet

    def _field_to_pyart_field(self, field, field_scale):
        pyart_field = {
            'coordinates': 'elevation azimuth range',
            'data': np.ma.masked_array(field),
            'long_name': self._variable_name_lookup[field_scale['name']][1],
            # This needs to be fixed eventually
            'standard_name': self._variable_name_lookup[field_scale['name']][0],
            'units': field_scale['units'],
            'valid_max': field_scale['max_val'],
            'valid_min': field_scale['min_val']
        }

        return pyart_field

    def _process_data_blocks(self):
        darray = np.array(self.data_array)
        cvi = 0
        for cv in range(0, len(self.field_scale_list)):
            if(bool(self.bit_mask & 2 ** cv)):
                fsl = self.field_scale_list[cv]

                if(self.field_scale_list[cv]['data_size'] == 'f'):
                    dat = darray[:, cvi::self.field_count].copy()
                    dat[dat == 0] = np.nan
                else:
                    dat = darray[:, cvi::self.field_count].copy()
                    dat[dat == 0] = np.nan
                    # I need to come back and do whole np.ma thing
                    dat = (dat * fsl['dat_factor'] + fsl['dat_bias']) / \
                        float(fsl['fld_factor'])
                self.fields[fsl['name']] = dat
                cvi += 1

    def _format_string_from_bitmask(self, bitmask):
        format_str = ''
        for b in range(0, 38):
            if(bool(bitmask & 2 ** b)):
                format_str += self.field_scale_list[b]['data_size']
        return format_str

    def _num_fields_from_bitmask(self, bitmask):
        count = 0
        for b in range(0, 38):
            if(bool(bitmask & 2 ** b)):
                count += 1
        return count

    def _parse_field_struct_packet(self, packet):
        packet['name'] = packet['name'].rstrip('\x00')
        packet['units'] = packet['units'].rstrip('\x00')
        packet['descr'] = packet['descr'].rstrip('\x00')
        packet['data_size'] = self.format_length_lookup_list[packet['format']]
        self.field_scale_list[packet['bit_mask_pos']] = packet

    format_length_lookup_list = [
        'c',
        'Q',
        'f',
        'H'
    ]

    arch_file_hdr_t = (
        "version",
        "creator_version",
        "creator_id",
        "sweep_table_offset"
    )

    arch_file_hdr_fstring = "2I32sQ"

    arch_housekeeping_t = (
        "id",
        "length",
        "rest_of_packet",
        # Need to finish this
    )

    arch_ray_header = (
        'azimuth',
        'elevation',
        'azimuth_width',
        'elevation_width',
        'gates',
        'beam_index',
        'ns_time',
        'time',
        'bit_mask',
        'ray_number',
        'num_pulses'
    )

    arch_ray_fstring = '4f2HI2Q2I'

    packet_type = {
        '0x00010000': 'ARCH_FORMAT_VERSION',
        '0x5aa80001': 'ARCH_ID_CONTROL',
        '0x5aa80002': 'ARCH_ID_FIELD_SCALE',
        '0x5aa80003': 'ARCH_ID_RAY_HDR',
        '0x5aa80004': 'ARCH_ID_FILE_HDR',
        '0x5aa80005': 'ARCH_ID_SWEEP_BLOCK',
        '0x5aa50003': 'HSK_ID_PROCESSOR_INFO',
        '0x5aa50001': 'HSK_ID_RADAR_INFO',
        '0x5aa50002': 'HSK_ID_SCAN_SEG'
    }

    field_scale_t = (
        'format',
        'min_val',
        'max_val',
        'bit_mask_pos',
        'type_hint',
        'fld_factor',
        'dat_factor',
        'dat_bias',
        'name',
        'units',
        'descr'
    )
    field_scale_fstring = 'i2f5i32s32s128s'

    radar_info_t = (
        'radar_name',
        'latitude',
        'longitude',
        'altitude',
        'beamwidth',
        'wavelength_cm',
        'gain_ant_h',
        'gain_ant_v',
        'zdr_cal_base',
        'phidp_rot',
        'base_radar_constant',
        'power_measurement_loss_h',
        'power_measurement_loss_v',
        'zdr_cal_base_vhs',
        'test_power_h',
        'test_power_v',
        'dc_loss_h',
        'dc_loss_v'
    )
    radar_info_fstring = '32s22f'

    processor_info_t = (
        'polarization_mode',
        'processing_mode',
        'pulse_type',
        'test_type',
        'integration_cycle_pulses',
        'clutter_filter_number',
        'range_gate_averaging',
        'indexed_beam_width',
        'gate_spacing',
        'prt_usec',
        'range_start',
        'range_stop',
        'max_gate',
        'test_power',
        'test_pulse_range',
        'test_pulse_length',
        'prt2',
        'range_offset'
    )
    processor_info_fstring = '4i3I5fI7f'

    scan_seg = (
        'az_manual',
        'el_manual',
        'az_start',
        'el_start',
        'scan_rate',
        'segname',
        'opt',
        'follow_mode',
        'scan_type',
        'scan_flags',
        'volume_num',
        'sweep_num',
        'time_limit',
        'webtilt',
        'left_limit',
        'right_limit',
        'up_limit',
        'down_limit',
        'step',
        'max_sweeps',
        'filter_break_sweep',
        'clutter_filter1',
        'clutter_filter2',
        'project',
        'current_fixed_angle'
    )
    scan_seg_fstring = '5f24s3i5I5f4I16sf'

    _scan_mode_names = ['ppi', 'rhi', 'fixed',
                        'manual ppi', 'manual rhi', 'idle']

    _range_dict = {'axis': 'radial_range_coordinate',
                   'comment': 'Coordinate variable for range. Range to center of each bin.',
                   'long_name': 'range_to_measurement_volume',
                   'meters_to_center_of_first_gate': 0.0,
                   'spacing_is_constant': 'true',
                   'standard_name': 'projection_range_coordinate',
                   'units': 'meters'}

    _time_dict = {'calendar': 'gregorian',
                  'comment': 'Coordinate variable for time. Time at the center of each ray, in fractional seconds since the global variable time_coverage_start',
                  'long_name': 'time_in_seconds_since_volume_start',
                  'standard_name': 'time',
                  'units': 'seconds since 1970-01-01 00:00 UTC'}

    def _pack_azimuth(self, azimuth):
        return {
            'data': self.azimuth,
            'axis': 'radial_azimuth_coordinate',
            'comment': 'Azimuth of antenna relative to true north',
            'long_name': 'azimuth_angle_from_true_north',
            'standard_name': 'beam_azimuth_angle',
            'units': 'degrees'
        }

    def _pack_elevation(self, elevation):
        return {
            'data': self.elevation,
            'axis': 'radial_elevation_coordinate',
            'comment': 'Elevation of antenna relative to the horizontal palne',
            'long_name': 'elevation_angle_from_horizontal_plane',
            'standard_name': 'beam_elevation_angle',
            'units': 'degrees'
        }

    def _pack_fixed_angle(self, fixed_angle):
        return {
            'data': self.fixed_angle,
            'long_name': 'Target angle for sweep',
            'standard_name': 'target_fixed_angle',
            'units': 'degrees'
        }

    _variable_name_lookup = {
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
        'H lag 0 cx': ('H lag 0 cx', 'absolute_value_of_lag_0_cross_correlation_h'),
        'H Im(lag 1)': ('H Im(lag 1)', 'imaginary_part_of_lag_1_correlation_h'),
        'H Re(lag 2)': ('H Re(lag 2)', 'real_part_of_lag_2_correlation_h'),
        'V lag 0 cx': ('V lag 0 cx', 'absolute_value_of_lag_0_cross_correlation_v'),
        'V Re(lag 1)': ('V Re(lag 1)', 'real_part_of_lag_1_correlation_v'),
        'V Im(lag 2)': ('V Im(lag 2)', 'imaginary_part_of_lag_2_correlation_v'),
        'HV lag 0 I': ('HV lag 0 I', 'real_part_of_cross_channel_correlation_at_lag_0'),
        'HV lag 0 Q': ('HV lag 0 Q', 'imaginary_part_of_cross_channel_correlation_at_lag_0'),
        'VAvgI': ('VAvgI', 'v_average_inphase'),
        'HAvgI': ('HAvgI', 'h_average_inphase'),
        '\xcf\x81 HCX': ('RHOHCX', 'lag_0_h_co_to_cross_correlation'),
        '\xcf\x81 VCX': ('RHOVCX', 'lag_0_v_co_to_cross_correlation'),
    }
