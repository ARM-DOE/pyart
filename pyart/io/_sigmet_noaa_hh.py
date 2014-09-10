"""
Functions for reading the Sigmet files from NOAA Hurricane hunter
airborne radars
"""

import numpy as np


def decode_noaa_hh_hdr(raw_extended_headers, azimuth, elevation,
                       position_source='irs', heading_source='irs'):
    """
    Extract parameters from NOAA Hurricane hunter radar extended headers.
    """
    xhdr = np.rec.fromstring(raw_extended_headers[..., :68].tostring(),
                             dtype=list(P3_EXTENDED_HEADER))

    # rotation and tilt from azimuth/elevation angles
    rotation_data = 90. - elevation['data'].copy()
    rotation_data[rotation_data < 0] += 360.
    tilt_data = azimuth['data'].copy()
    tilt_data[tilt_data > 180] -= 360.
    rotation = {'data': rotation_data}
    tilt = {'data': tilt_data}

    # airborne parameters
    if heading_source == 'irs':
        heading_data = _bam16_to_angle(xhdr['irs_heading'])
        roll_data = _bam16_to_angle(xhdr['irs_roll'])
        pitch_data = _bam16_to_angle(xhdr['irs_pitch'])
        drift_data = _bam16_to_angle(xhdr['irs_drift'])
    elif heading_source == 'aamps':
        heading_data = _bam16_to_angle(xhdr['aamps_heading'])
        roll_data = _bam16_to_angle(xhdr['aamps_roll'])
        pitch_data = _bam16_to_angle(xhdr['aamps_pitch'])
        drift_data = _bam16_to_angle(xhdr['aamps_drift'])
    else:
        raise ValueError('Unknown heading_source')
    heading = {'data': heading_data}
    roll = {'data': roll_data}
    pitch = {'data': pitch_data}
    drift = {'data': drift_data}

    # positions: latitude, longitude, altitude
    if position_source == 'gps':
        lat_data = _bam32_to_angle(xhdr['gps_lat'])
        lon_data = _bam32_to_angle(xhdr['gps_long'])
        alt_data = xhdr['gps_alt'] / 100.
    elif position_source == 'aamps':
        lat_data = _bam32_to_angle(xhdr['aamps_lat'])
        lon_data = _bam32_to_angle(xhdr['aamps_long'])
        alt_data = xhdr['aamps_alt'] / 100.
    elif position_source == 'irs':
        lat_data = _bam32_to_angle(xhdr['irs_lat'])
        lon_data = _bam32_to_angle(xhdr['irs_long'])
        alt_data = xhdr['gps_alt'] / 100.
    else:
        raise ValueError('Invalid position_source')

    latitude = {'data': lat_data}
    longitude = {'data': lon_data}
    altitude = {'data': alt_data}

    # georeferenced azimuth and elevation
    az, elev = _compute_azimuth_elevation(
        roll_data, pitch_data, heading_data, drift_data, rotation_data,
        tilt_data)
    azimuth['data'] = az
    elevation['data'] = elev

    extended_header_params = {
        'heading': heading,
        'roll': roll,
        'pitch': pitch,
        'drift': drift,
        'rotation': rotation,
        'tilt': tilt,
        'georefs_applied': None}
    return (latitude, longitude, altitude, extended_header_params)


def _compute_azimuth_elevation(roll, pitch, heading, drift,
                               rotation, tilt):
    """ Compute georefered azimuth and elevation angles. """
    # Adapted from SigmetRadxFile::_computeAzEl method of Radx found in
    # SigmetRadxFile.cc
    # Based on transforms in Wen-Chau et al, JTech, 1994, 11, 572-578.
    degree_to_radian = np.pi / 180.

    R = roll * degree_to_radian
    P = pitch * degree_to_radian
    H = heading * degree_to_radian
    D = drift * degree_to_radian
    T = H + D

    sinP = np.sin(P)
    cosP = np.cos(P)
    sinD = np.sin(D)
    cosD = np.cos(D)

    theta_a = rotation * degree_to_radian
    tau_a = tilt * degree_to_radian
    sin_tau_a = np.sin(tau_a)
    cos_tau_a = np.cos(tau_a)
    sin_theta_rc = np.sin(theta_a + R)
    cos_theta_rc = np.cos(theta_a + R)

    xsubt = (cos_theta_rc * sinD * cos_tau_a * sinP +
             cosD * sin_theta_rc * cos_tau_a -
             sinD * cosP * sin_tau_a)

    ysubt = (-cos_theta_rc * cosD * cos_tau_a * sinP +
             sinD * sin_theta_rc * cos_tau_a +
             cosP * cosD * sin_tau_a)

    zsubt = (cosP * cos_tau_a * cos_theta_rc + sinP * sin_tau_a)
    lambda_t = np.arctan2(xsubt, ysubt)
    azimuth = np.fmod(lambda_t + T, 2 * np.pi)
    elevation = np.arcsin(zsubt)

    return azimuth, elevation


def _bam16_to_angle(arr):
    """ Convert a 16-bit Binary angle to a angle """
    return arr * 360. / 65536.0


def _bam32_to_angle(arr):
    return arr * 360. / 4294967296.0

# Structure types
UINT16 = 'H'
INT16 = 'h'
BAM16 = 'h'
INT32 = 'i'
BAM32 = 'i'

P3_EXTENDED_HEADER = (
    ('msecs_since_sweep_start', INT32),
    ('calib_signal_level', INT16),
    ('nbytes_in_header', INT16),
    ('__pad_1', UINT16),
    ('gps_age', UINT16),            # Time in milliseconds since last GPS Input
    ('irs_age', UINT16),            # Time in milliseconds since last IRS Input
    ('aamps_age', UINT16),          # Time in milliseconds since last
                                    # AAMPS Input
    ('gps_lat', BAM32),             # GPS latitude (BAM)
    ('gps_long', BAM32),            # GPS Longitude (BAM)
    ('gps_alt', INT32),             # GPS Altitude (cm)
    ('gps_vel_e', INT32),           # GPS Ground Speed East (cm/second)
    ('gps_vel_n', INT32),           # GPS Ground Speed North (cm/second)
    ('gps_vel_v', INT32),           # GPS Ground Speed Up (cm/second)
    ('irs_lat', BAM32),             # IRS latitude (BAM)
    ('irs_long', BAM32),            # IRS Longitude (BAM)
    ('irs_vel_e', INT32),           # IRS Ground Speed East (cm/second)
    ('irs_vel_n', INT32),           # IRS Ground Speed North (cm/second)
    ('irs_vel_v', INT32),           # IRS Ground Speed Up (cm/second)
    ('irs_pitch', BAM16),           # IRS Pitch (BAM)
    ('irs_roll', BAM16),            # IRS Roll (BAM)
    ('irs_heading', BAM16),         # IRS Heading (BAM)
    ('irs_drift', BAM16),           # IRS Drift (BAM)
    ('irs_tru_track', BAM16),       # IRS True Track (BAM)
    ('irs_pitch_r', BAM16),         # IRS Pitch rate (BAM/sec)
    ('irs_roll_r', BAM16),          # IRS Roll rate (BAM/sec)
    ('irs_yaw_r', BAM16),           # IRS Yaw rate (BAM/sec)
    ('irs_wind_vel', INT32),        # IRS Wind speed (cm/second)
    ('irs_wind_dir', BAM16),        # IRS Wind direction (BAM)
    ('__pad_2', UINT16),
    ('aamps_lat', BAM32),           # AAMPS latitude (BAM)
    ('aamps_long', BAM32),          # AAMPS Longitude (BAM)
    ('aamps_alt', INT32),           # AAMPS Altitude (cm)
    ('aamps_ground_vel', INT32),    # AAMPS Ground Speed East (cm/second)
    ('aamps_time_stamp', INT32),    # AAMPS Timestamp in UTC
                                    # (seconds since the epoch)
    ('aamps_vel_v', INT32),         # AAMPS Vertical Velocity (cm/second)
    ('aamps_pitch', BAM16),         # AAMPS Pitch (BAM)
    ('aamps_roll', BAM16),          # AAMPS Roll (BAM)
    ('aamps_heading', BAM16),       # AAMPS Heading (BAM)
    ('aamps_drift', BAM16),         # AAMPS Drift (BAM)
    ('aamps_track', BAM16),         # AAMPS Track (BAM)
    ('__pad_4', UINT16),
    ('aamps_radalt_val', INT32),    # AAMPS Radar Altitude (cm)
    ('aamps_wind_vel', INT32),      # AAMPS Wind Speed (cm/second)
    ('aamps_wind_dir', BAM16),      # AAMPS Wind direction (BAM)
    ('__pad_5', UINT16),
    ('aamps_wind_vel_v', INT32),    # AAMPS Wind Speed Up (cm/second)
)
