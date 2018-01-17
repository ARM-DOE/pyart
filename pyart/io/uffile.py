"""
pyart.io.uffile
===============

Low level class for reading Universal Format (UF) files.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    UFFile
    UFRay

.. autosummary::
    :toctree: generated/

    _structure_size
    _unpack_from_buf
    _unpack_structure

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

import struct
import datetime

import numpy as np


class UFFile(object):
    """
    A class for reading data from Universal Format (UF) files.

    Parameters
    ----------
    filename : str or file-like
        Filename or file-like object containing data in Universal format (UF).

    Attributes
    ----------
    rays : list of UFRay objects
        List of rays within the UF file.
    nrays, nsweeps : int
        Number of rays and sweep in the file.
    ray_sweep_numbers : array
        Sweep number of each ray in the file.
    first_ray_in_sweep, last_ray_in_sweep : array
        Indices of the first and last ray in each sweep.

    """

    def __init__(self, filename):
        """ initialize. """

        # open the file if file object not passed
        if hasattr(filename, 'read'):
            fobj = filename
        else:
            fobj = open(filename, 'rb')
        self._fh = fobj

        # UF files come in three 'flavors' depending upon the size of the
        # padding around each record.  True UF files contain no padding
        # and start with the mandatory header of the first ray.  Other UF
        # files contain a 2 or 4-byte padding immediately before and after
        # each record.  The values in this padding can used to determine the
        # size of each record, but is not used here, rather the size indicated
        # by the 'record_length' structure elements is used.

        # determine padding around records
        buf = fobj.read(8)
        try:
            padding = buf.index(b'UF')
        except ValueError:
            raise IOError('file in not a valid UF file')

        # read in the records, store as a list of rays
        self.rays = []
        while len(buf) == 8:  # read until EOF reached

            # record size stored as a 2-byte int start at byte 2
            record_size = struct.unpack('>h', buf[padding+2:padding+4])[0] * 2

            # read in full record
            bytes_read = len(buf) - padding
            bytes_to_read = record_size - bytes_read
            record = buf[-bytes_read:] + fobj.read(bytes_to_read)

            # convert record into UFRay
            self.rays.append(UFRay(record))

            # read post record padding
            fobj.read(padding)

            # read in the first eight bytes of the next record
            buf = fobj.read(8)

        # determine volume size statistics
        self.nrays = len(self.rays)

        # determine sweep information
        self.ray_sweep_numbers = self._get_ray_sweep_numbers()
        self.nsweeps = len(np.unique(self.ray_sweep_numbers))
        first_ray_in_sweep, last_ray_in_sweep = self._get_sweep_limits()
        self.first_ray_in_sweep = first_ray_in_sweep
        self.last_ray_in_sweep = last_ray_in_sweep

    def close(self):
        """ Close the file. """
        self._fh.close()

    def _get_ray_sweep_numbers(self):
        """ Return an array of the sweep_number stored in each ray. """
        ray_sweep_numbers = np.empty((self.nrays, ), dtype='int32')
        for i, ray in enumerate(self.rays):
            ray_sweep_numbers[i] = ray.mandatory_header['sweep_number']
        return ray_sweep_numbers

    def _get_sweep_limits(self):
        """ Return arrays of indices of first and last ray in each sweep. """
        first_ray_in_sweep = np.empty(self.nsweeps, dtype='int32')
        last_ray_in_sweep = np.empty(self.nsweeps, dtype='int32')
        unique_sweep_numbers = np.unique(self.ray_sweep_numbers)
        for i, sweep_number in enumerate(unique_sweep_numbers):
            matches = np.where(self.ray_sweep_numbers == sweep_number)
            first_ray_in_sweep[i] = matches[0][0]
            last_ray_in_sweep[i] = matches[0][-1]
        return first_ray_in_sweep, last_ray_in_sweep

    def get_field_data(self, field_number):
        """ Return a 2D array of scale/masked field data for the volume. """
        # Assumes that no rays contain more gates than the first ray and
        # that the missing_data_value and scale_factor are identical for all
        # rays.  Additional the order and number of the fields are assumed to
        # be identical between rays.
        first_ray = self.rays[0]
        ngates = len(first_ray.field_raw_data[field_number])
        missing_data_value = first_ray.mandatory_header['missing_data_value']
        scale_factor = first_ray.field_headers[field_number]['scale_factor']

        raw_data = np.empty((self.nrays, ngates), 'int16')
        for i, ray in enumerate(self.rays):
            ray_data = ray.field_raw_data[field_number]
            bins = len(ray_data)
            raw_data[i, :bins] = ray.field_raw_data[field_number]
            raw_data[i, bins:] = missing_data_value

        data = raw_data / float(scale_factor)
        mask = raw_data == missing_data_value
        return np.ma.masked_array(data, mask)

    def get_azimuths(self):
        """ Return an array of azimuth angles for each ray in degrees. """
        azimuth = np.empty((self.nrays, ), dtype='float32')
        for i, ray in enumerate(self.rays):
            azimuth[i] = ray.mandatory_header['azimuth'] / 64.
        return azimuth

    def get_elevations(self):
        """ Return an array of elevation angles for each ray in degrees. """
        elevation = np.empty((self.nrays, ), dtype='float32')
        for i, ray in enumerate(self.rays):
            elevation[i] = ray.mandatory_header['elevation'] / 64.
        return elevation

    def get_sweep_rates(self):
        """ Return an array of sweep rates for each ray in degrees/sec. """
        sweep_rates = np.empty((self.nrays, ), dtype='float32')
        for i, ray in enumerate(self.rays):
            sweep_rates[i] = ray.mandatory_header['sweep_rate'] / 64.
        return sweep_rates

    def get_pulse_widths(self):
        """ Return an array of pulse widths for each ray in meters. """
        pulse_widths = np.empty((self.nrays, ), dtype='float32')
        for i, ray in enumerate(self.rays):
            pulse_widths[i] = ray.field_headers[0]['pulse_width_m']
        return pulse_widths

    def get_prts(self):
        """ Return an array of prts for each ray in microseconds. """
        prts = np.empty((self.nrays, ), dtype='float32')
        for i, ray in enumerate(self.rays):
            prts[i] = ray.field_headers[0]['prt_ms']
        return prts

    def get_nyquists(self):
        """
        Return an array of nyquist velocities for each ray in m/s.

        Returns None if nyquist velocities cannot be determined for all rays.
        """
        field_headers = self.rays[0].field_headers
        try:
            field_idx = ['nyquist' in fh for fh in field_headers].index(True)
        except ValueError:
            return None  # True not in list
        nyquist = np.empty((self.nrays, ), dtype='float32')
        for i, ray in enumerate(self.rays):
            scale = ray.field_headers[field_idx]['scale_factor']
            try:
                nyquist[i] = ray.field_headers[field_idx]['nyquist'] / scale
            except KeyError:
                return None  # nyquist not in field header
        return nyquist

    def get_sweep_fixed_angles(self):
        """ Return an array of fixed angles for each sweep in degrees. """
        fixed = np.empty((self.nsweeps, ), dtype='float32')
        for i, ray_num in enumerate(self.first_ray_in_sweep):
            fixed[i] = self.rays[ray_num].mandatory_header['fixed_angle'] / 64.
        return fixed

    def get_sweep_polarizations(self):
        """ Return an array of polarization modes for each sweep. """
        modes = []
        for ray_num in self.first_ray_in_sweep:
            ray = self.rays[ray_num]
            polarization = ray.field_headers[0]['polarization']
            if polarization > 3:
                polarization = 3
            modes.append(POLARIZATION_STR[polarization])
        return np.array(modes)

    def get_datetimes(self):
        """ Return a list of datetimes for each ray. """
        return [ray.get_datetime() for ray in self.rays]


class UFRay(object):
    """
    A class for reading data from a single ray (record) in a UF file.

    Parameters
    ----------
    record : str
        Byte string containing the binary data for a UF ray.

    Attributes
    ----------
    mandatory_header : dic
        Mandatory header.
    optional_header : dic or None
        Optional header or None if no optional header exists in the record.
    data_header : dic
        Data header.
    field_positions : list
        List of dictionaries containing the data type and data position.
    field_headers : list
        List of field header dictionaries for all fields in the ray.
    field_raw_data : list
        List containing array of raw field data for each field in the ray.
    _buf : str
        Bytes which make up the record.

    """

    def __init__(self, record):
        """ Initalize the object. """

        self._buf = record

        # read in the mandatory header
        self.mandatory_header = _unpack_from_buf(
            self._buf, 0, UF_MANDATORY_HEADER)

        # read in optional header (if present)
        self.optional_header = None
        if self.mandatory_header['offset_optional_header'] != 0:
            offset = (self.mandatory_header['offset_optional_header'] - 1) * 2
            self.optional_header = _unpack_from_buf(
                self._buf, offset, UF_OPTIONAL_HEADER)

        # read in data header
        offset = (self.mandatory_header['offset_data_header'] - 1) * 2
        self.data_header = _unpack_from_buf(self._buf, offset, UF_DATA_HEADER)

        # read in field position information
        self.field_positions = [
            _unpack_from_buf(self._buf, offset + 6 + i*4, UF_FIELD_POSITION)
            for i in range(self.data_header['record_nfields'])]

        # read field headers and data
        self.field_headers = []
        self.field_raw_data = [self.get_field_data(i) for i in
                               range(self.data_header['record_nfields'])]

        return

    def get_field_data(self, field_number):
        """
        Return array of raw data for a particular field in the ray.

        Field header is appended to the list in the field_headers attribute.
        """
        position = self.field_positions[field_number]
        offset = (position['offset_field_header'] - 1) * 2
        field_header = _unpack_from_buf(self._buf, offset, UF_FIELD_HEADER)
        self.field_headers.append(field_header)
        data_offset = (field_header['data_offset'] - 1) * 2

        # read in field specific parameters
        if position['data_type'] in [b'VF', b'VE', b'VR', b'VT', b'VP']:
            if (data_offset - offset) == 42:
                vel_header = _unpack_from_buf(self._buf, offset+38, UF_FSI_VEL)
                field_header.update(vel_header)

        data_str = self._buf[data_offset:data_offset+field_header['nbins']*2]
        raw_data = np.frombuffer(data_str, dtype='>i2')
        return raw_data

    def get_datetime(self):
        """ Return a datetime object for the ray. """
        year = self.mandatory_header['year']
        if year < 1900:
            year += 2000   # years after 2000, 11 -> 2011
        month = self.mandatory_header['month']
        day = self.mandatory_header['day']
        hour = self.mandatory_header['hour']
        minute = self.mandatory_header['minute']
        second = self.mandatory_header['second']
        if hour == 24:
            # Some UF writers incorrectly specify midnight as 24:00:00
            # rather than 00:00:00.  Handle this case explicitly
            assert minute == 0
            assert second == 0
            hour = 23
            minute = 59
            second = 59
            dt = datetime.datetime(year, month, day, hour, minute, second)
            return dt + datetime.timedelta(seconds=1)

        return datetime.datetime(year, month, day, hour, minute, second)

    def get_location(self):
        """ Return the latitude, longitude and height of the ray. """
        lat_deg = self.mandatory_header['latitude_degrees']
        lat_min = self.mandatory_header['latitude_minutes']
        lat_sec = self.mandatory_header['latitude_seconds'] / 64.
        latitude = lat_deg + (lat_min + lat_sec / 60.) / 60.

        lon_deg = self.mandatory_header['longitude_degrees']
        lon_min = self.mandatory_header['longitude_minutes']
        lon_sec = self.mandatory_header['longitude_seconds'] / 64.
        longitude = lon_deg + (lon_min + lon_sec / 60.) / 60.

        height = self.mandatory_header['height_above_sea_level']

        return latitude, longitude, height


def _structure_size(structure):
    """ Find the size of a structure in bytes. """
    return struct.calcsize('>' + ''.join([i[1] for i in structure]))


def _unpack_from_buf(buf, pos, structure):
    """ Unpack a structure from a buffer. """
    size = _structure_size(structure)
    return _unpack_structure(buf[pos:pos + size], structure)


def _unpack_structure(string, structure):
    """ Unpack a structure from a string """
    fmt = '>' + ''.join([i[1] for i in structure])  # UF is big-endian
    lst = struct.unpack(fmt, string)
    return dict(zip([i[0] for i in structure], lst))


# The Universal file format was originally described in the report:
#
# Barnes, Stanley L. Report on a meeting to establish a common Doppler radar
# data exchange format. The Bulletin of the American Meteorological Society,
# Vol 61, No. 11, pp1401-1404. Nov 1980
#
# This report is difficult to obtain as it is not available electronically from
# the American Meteorological Society.
#
# A modern descriptions of the format that can be found at:
# http://www.ral.ucar.edu/projects/titan/docs/radial_formats/UfDoc.txt
#
# Additionally Appendix C of the Vaisala Programmer's Manual contains a
# description of the format.

# UF structures

POLARIZATION_STR = ['horizontal', 'vertical', 'circular', 'elliptical']

INT16 = 'h'

UF_MANDATORY_HEADER = (
    ('uf_string', '2s'),
    ('record_length', INT16),
    ('offset_optional_header', INT16),
    ('offset_local_use_header', INT16),
    ('offset_data_header', INT16),
    ('record_number', INT16),
    ('volume_number', INT16),
    ('ray_number', INT16),
    ('ray_record_number', INT16),
    ('sweep_number', INT16),
    ('radar_name', '8s'),
    ('site_name', '8s'),
    ('latitude_degrees', INT16),
    ('latitude_minutes', INT16),
    ('latitude_seconds', INT16),
    ('longitude_degrees', INT16),
    ('longitude_minutes', INT16),
    ('longitude_seconds', INT16),
    ('height_above_sea_level', INT16),
    ('year', INT16),
    ('month', INT16),
    ('day', INT16),
    ('hour', INT16),
    ('minute', INT16),
    ('second', INT16),
    ('time_zone', '2s'),
    ('azimuth', INT16),
    ('elevation', INT16),
    ('sweep_mode', INT16),
    ('fixed_angle', INT16),
    ('sweep_rate', INT16),
    ('generation_year', INT16),
    ('generation_month', INT16),
    ('generation_day', INT16),
    ('generation_facility_name', '8s'),
    ('missing_data_value', INT16),
)

UF_OPTIONAL_HEADER = (
    ('project_name', '8s'),
    ('baseline_azimuth', INT16),
    ('baseline_elevation', INT16),
    ('volume_hour', INT16),
    ('volume_minute', INT16),
    ('volume_second', INT16),
    ('tape_name', '8s'),
    ('flag', INT16)
)

UF_DATA_HEADER = (
    ('ray_nfields', INT16),
    ('ray_nrecords', INT16),
    ('record_nfields', INT16),
)

UF_FIELD_POSITION = (
    ('data_type', '2s'),
    ('offset_field_header', INT16),
)

UF_FIELD_HEADER = (
    ('data_offset', INT16),
    ('scale_factor', INT16),
    ('range_start_km', INT16),
    ('range_start_m', INT16),
    ('range_spacing_m', INT16),
    ('nbins', INT16),
    ('pulse_width_m', INT16),
    ('beam_width_h', INT16),    # degrees * 64
    ('beam_width_v', INT16),    # degrees * 64
    ('bandwidth', INT16),       # Reciever bandwidth in MHz * 16
    ('polarization', INT16),    # 1: hort, 2: vert 3: circular, 4: ellip
    ('wavelength_cm', INT16),   # cm * 64
    ('sample_size', INT16),
    ('threshold_data', '2s'),
    ('threshold_value', INT16),
    ('scale', INT16),
    ('edit_code', '2s'),
    ('prt_ms', INT16),
    ('bits_per_bin', INT16),    # Must be 16
)

UF_FSI_VEL = (
    ('nyquist', INT16),
    ('spare', INT16),
)

# This structure is defined but not used in Py-ART
# No sample file which contain the structure could be found.
UF_FSI_DM = (
    ('radar_constant', INT16),
    ('noise_power', INT16),
    ('reciever_gain', INT16),
    ('peak_power', INT16),
    ('antenna_gain', INT16),
    ('pulse_duration', INT16),
)
