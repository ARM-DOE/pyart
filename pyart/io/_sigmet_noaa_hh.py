"""
Functions needed for reading Sigmet files from the airborne radar located on
NOAA's Hurricane Hunter aircraft.

"""

import numpy as np
from ._sigmetfile import bin4_to_angle, bin2_to_angle


def _decode_noaa_hh_hdr(
        raw_extended_headers, filemetadata, azimuth, elevation,
        position_source='irs', heading_source='irs'):
    """
    Extract data from Sigmet extended headers produced by NOAA
    Hurricane Hunter airborne radars.

    Parameters
    ----------
    raw_extended_headers : ndarray
        Raw Sigmet extended headers.
    filemetadata : FileMetadata
        FileMetadata class from which metadata will be derived.
    azimuth : dict
        Dictionary of azimuth angles recorded in Sigmet file.
    elevation : dict
        Dictionary of elevation angles recorded in Sigmet file.
    position_source: {'irs', 'gps', 'aamps'}, optional
        Instrument from which to derive position parameters.
    heading_source: {'irs', 'aamps'}
        Instrument from which to derive heading parameters.

    Returns
    -------
    latitude : dict
        Dictionary containing latitude data and metadata.
    longitude : dict
        Dictionary containing longitude data and metadata.
    altitude : dict
        Dictionary containing altitude data and metadata.
    heading_params : dict
        Dictionary of dictionary containing aircraft heading data and
        metadata. Contains 'heading', 'roll', pitch', 'drift', 'rotation',
        'tilt' and 'georefs_applied' dictionaries.

    """
    xhdr = np.frombuffer(raw_extended_headers[..., :68].tostring(),
                         dtype=list(NOAA_HH_EXTENDED_HEADER))

    # rotation and tilt from azimuth/elevation angles
    rotation = filemetadata('rotation')
    tilt = filemetadata('tilt')

    rotation_data = 90. - elevation['data'].copy()
    rotation_data[rotation_data < 0] += 360.
    rotation['data'] = rotation_data

    tilt_data = azimuth['data'].copy()
    tilt_data[tilt_data > 180] -= 360.
    tilt['data'] = tilt_data

    # airborne parameters
    heading = filemetadata('heading')
    roll = filemetadata('roll')
    pitch = filemetadata('pitch')
    drift = filemetadata('drift')

    if heading_source == 'irs':
        heading_data = bin2_to_angle(xhdr['irs_heading'])
        roll_data = bin2_to_angle(xhdr['irs_roll'])
        pitch_data = bin2_to_angle(xhdr['irs_pitch'])
        drift_data = bin2_to_angle(xhdr['irs_drift'])
    elif heading_source == 'aamps':
        heading_data = bin2_to_angle(xhdr['aamps_heading'])
        roll_data = bin2_to_angle(xhdr['aamps_roll'])
        pitch_data = bin2_to_angle(xhdr['aamps_pitch'])
        drift_data = bin2_to_angle(xhdr['aamps_drift'])
    else:
        raise ValueError('Unknown heading_source')

    heading['data'] = heading_data
    roll['data'] = roll_data
    pitch['data'] = pitch_data
    drift['data'] = drift_data

    # georeferenced azimuth and elevation
    az, elev = _georeference_yprime(
        roll_data, pitch_data, heading_data, drift_data, rotation_data,
        tilt_data)
    azimuth['data'] = az
    elevation['data'] = elev
    georefs_applied = filemetadata('georefs_applied')
    georefs_applied['data'] = np.ones(az.shape, dtype='int8')

    # positions: latitude, longitude, altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')

    if position_source == 'gps':
        lat_data = bin4_to_angle(xhdr['gps_lat'])
        lon_data = bin4_to_angle(xhdr['gps_long'])
        alt_data = xhdr['gps_alt'] / 100.
    elif position_source == 'aamps':
        lat_data = bin4_to_angle(xhdr['aamps_lat'])
        lon_data = bin4_to_angle(xhdr['aamps_long'])
        alt_data = xhdr['aamps_alt'] / 100.
    elif position_source == 'irs':
        lat_data = bin4_to_angle(xhdr['irs_lat'])
        lon_data = bin4_to_angle(xhdr['irs_long'])
        alt_data = xhdr['gps_alt'] / 100.
    else:
        raise ValueError('Invalid position_source')

    latitude['data'] = lat_data
    longitude['data'] = lon_data
    altitude['data'] = alt_data

    extended_header_params = {
        'heading': heading,
        'roll': roll,
        'pitch': pitch,
        'drift': drift,
        'rotation': rotation,
        'tilt': tilt,
        'georefs_applied': georefs_applied}
    return (latitude, longitude, altitude, extended_header_params)


def _georeference_yprime(roll, pitch, heading, drift, rotation, tilt):
    """
    Compute georeferenced azimuth and elevation angles for a Y-prime radar.

    This is the georeferencing needed for the tail doppler radar on the
    NOAA P3 aircraft.
    """
    # Adapted from Radx's SigmetRadxFile::_computeAzEl method found in
    # SigmetRadxFile.cc
    # Transforms defined in Wen-Chau Lee et al, JTech, 1994, 11, 572-578.

    # Convert to radians and use variable names from Wen-Chau Lee paper
    R = np.radians(roll)        # roll
    P = np.radians(pitch)       # pitch
    H = np.radians(heading)     # heading
    D = np.radians(drift)       # drift
    T = H + D                   # track
    theta_a = np.radians(rotation)
    tau_a = np.radians(tilt)

    # Eq. (9)
    x_t = (np.cos(theta_a + R) * np.sin(D) * np.cos(tau_a) * np.sin(P) +
           np.cos(D) * np.sin(theta_a + R) * np.cos(tau_a) -
           np.sin(D) * np.cos(P) * np.sin(tau_a))

    y_t = (-np.cos(theta_a + R) * np.cos(D) * np.cos(tau_a) * np.sin(P) +
           np.sin(D) * np.sin(theta_a + R) * np.cos(tau_a) +
           np.cos(P) * np.cos(D) * np.sin(tau_a))

    z_t = (np.cos(P) * np.cos(tau_a) * np.cos(theta_a + R) +
           np.sin(P) * np.sin(tau_a))

    # Eq. (12) and discussion after Eq. (17)
    lambda_t = np.arctan2(x_t, y_t)
    azimuth = np.fmod(lambda_t + T, 2 * np.pi)

    # Eq (17)
    elevation = np.arcsin(z_t)

    # convert to degrees and fix range
    azimuth = np.degrees(azimuth)
    azimuth[azimuth < 0] += 360.
    elevation = np.degrees(elevation)
    elevation[elevation > 180] -= 360.
    return azimuth, elevation

# NOAA Hurrican Hunter Sigmet Extended header structure
# scalar definitions
UINT16 = 'H'
INT16 = 'h'
BAM16 = 'h'
INT32 = 'i'
BAM32 = 'i'

NOAA_HH_EXTENDED_HEADER = (
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
