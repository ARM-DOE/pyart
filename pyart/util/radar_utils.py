"""
pyart.util.radar_utils
======================

Functions for working radar instances.

.. autosummary::
    :toctree: generated/

    is_vpt
    to_vpt
    join_radar

"""

from __future__ import print_function

import copy

import numpy as np
from netCDF4 import num2date, date2num

from ..config import get_fillvalue
from . import datetime_utils


def is_vpt(radar, offset=0.5):
    """
    Determine if a Radar appears to be a vertical pointing scan.

    This function only verifies that the object is a vertical pointing scan,
    use the :py:func:`to_vpt` function to convert the radar to a vpt scan
    if this function returns True.

    Parameters
    ----------
    radar : Radar
        Radar object to determine if
    offset : float
        Maximum offset of the elevation from 90 degrees to still consider
        to be vertically pointing.

    Returns
    -------
    flag : bool
        True if the radar appear to be verticle pointing, False if not.

    """
    # check that the elevation is within offset of 90 degrees.
    elev = radar.elevation['data']
    return np.all((elev < 90.0 + offset) & (elev > 90.0 - offset))


def to_vpt(radar, single_scan=True):
    """
    Convert an existing Radar object to represent a vertical pointing scan.

    This function does not verify that the Radar object contains a vertical
    pointing scan.  To perform such a check use :py:func:`is_vpt`.

    Parameters
    ----------
    radar : Radar
        Mislabeled vertical pointing scan Radar object to convert to be
        properly labeled.  This object is converted in place, no copy of
        the existing data is made.
    single_scan : bool, optional
        True to convert the volume to a single scan, any azimuth angle data
        is lost.  False will convert the scan to contain the same number of
        scans as rays, azimuth angles are retained.

    """
    if single_scan:
        nsweeps = 1
        radar.azimuth['data'][:] = 0.0
        seri = np.array([radar.nrays - 1], dtype='int32')
        radar.sweep_end_ray_index['data'] = seri
    else:
        nsweeps = radar.nrays
        # radar.azimuth not adjusted
        radar.sweep_end_ray_index['data'] = np.arange(nsweeps, dtype='int32')

    radar.scan_type = 'vpt'
    radar.nsweeps = nsweeps
    radar.target_scan_rate = None       # no scanning
    radar.elevation['data'][:] = 90.0

    radar.sweep_number['data'] = np.arange(nsweeps, dtype='int32')
    radar.sweep_mode['data'] = np.array(['vertical_pointing'] * nsweeps)
    radar.fixed_angle['data'] = np.ones(nsweeps, dtype='float32') * 90.0
    radar.sweep_start_ray_index['data'] = np.arange(nsweeps, dtype='int32')

    if radar.instrument_parameters is not None:
        for key in ['prt_mode', 'follow_mode', 'polarization_mode']:
            if key in radar.instrument_parameters:
                ip_dic = radar.instrument_parameters[key]
                ip_dic['data'] = np.array([ip_dic['data'][0]] * nsweeps)

    # Attributes that do not need any changes
    # radar.altitude
    # radar.altitude_agl
    # radar.latitude
    # radar.longitude

    # radar.range
    # radar.ngates
    # radar.nrays

    # radar.metadata
    # radar.radar_calibration

    # radar.time
    # radar.fields
    # radar.antenna_transition
    # radar.scan_rate
    return


def join_radar(radar1, radar2):
    """
    Combine two radar instances into one.

    Parameters
    ----------
    radar1 : Radar
        Radar object.
    radar2 : Radar
        Radar object.

    """
    # must have same gate spacing
    new_radar = copy.deepcopy(radar1)
    new_radar.azimuth['data'] = np.append(
        radar1.azimuth['data'], radar2.azimuth['data'])
    new_radar.elevation['data'] = np.append(
        radar1.elevation['data'], radar2.elevation['data'])
    new_radar.fixed_angle['data'] = np.append(
        radar1.fixed_angle['data'], radar2.fixed_angle['data'])
    new_radar.sweep_number['data'] = np.append(
        radar1.sweep_number['data'], radar2.sweep_number['data'])
    new_radar.sweep_start_ray_index['data'] = np.append(
        radar1.sweep_start_ray_index['data'],
        radar2.sweep_start_ray_index['data'] + radar1.nrays)
    new_radar.sweep_end_ray_index['data'] = np.append(
        radar1.sweep_end_ray_index['data'],
        radar2.sweep_end_ray_index['data'] + radar1.nrays)
    new_radar.nsweeps += radar2.nsweeps
    new_radar.sweep_mode['data'] = np.append(
        radar1.sweep_mode['data'], radar2.sweep_mode['data'])

    if len(radar1.range['data']) >= len(radar2.range['data']):
        new_radar.range['data'] = radar1.range['data']
    else:
        new_radar.range['data'] = radar2.range['data']
    new_radar.ngates = len(new_radar.range['data'])

    # to combine times we need to reference them to a standard
    # for this we'll use epoch time
    r1num = datetime_utils.datetimes_from_radar(radar1, epoch=True)
    r2num = datetime_utils.datetimes_from_radar(radar2, epoch=True)
    new_radar.time['data'] = date2num(
        np.append(r1num, r2num), datetime_utils.EPOCH_UNITS)
    new_radar.time['units'] = datetime_utils.EPOCH_UNITS
    new_radar.nrays = len(new_radar.time['data'])

    for var in new_radar.fields.keys():
        sh1 = radar1.fields[var]['data'].shape
        sh2 = radar2.fields[var]['data'].shape
        new_field_shape = (sh1[0] + sh2[0], max(sh1[1], sh2[1]))
        new_field = np.ma.zeros(new_field_shape)
        new_field[:] = np.ma.masked
        new_field.set_fill_value(get_fillvalue())
        new_field[0:sh1[0], 0:sh1[1]] = radar1.fields[var]['data']
        new_field[sh1[0]:, 0:sh2[1]] = radar2.fields[var]['data']
        new_radar.fields[var]['data'] = new_field

    # radar locations
    # TODO moving platforms - any more?
    if (len(radar1.latitude['data']) == 1 &
            len(radar2.latitude['data']) == 1 &
            len(radar1.longitude['data']) == 1 &
            len(radar2.longitude['data']) == 1 &
            len(radar1.altitude['data']) == 1 &
            len(radar2.altitude['data']) == 1):

        lat1 = float(radar1.latitude['data'])
        lon1 = float(radar1.longitude['data'])
        alt1 = float(radar1.altitude['data'])
        lat2 = float(radar2.latitude['data'])
        lon2 = float(radar2.longitude['data'])
        alt2 = float(radar2.altitude['data'])

        if (lat1 != lat2) or (lon1 != lon2) or (alt1 != alt2):
            ones1 = np.ones(len(radar1.time['data']), dtype='float32')
            ones2 = np.ones(len(radar2.time['data']), dtype='float32')
            new_radar.latitude['data'] = np.append(ones1 * lat1, ones2 * lat2)
            new_radar.longitude['data'] = np.append(ones1 * lon1, ones2 * lon2)
            new_radar.latitude['data'] = np.append(ones1 * alt1, ones2 * alt2)
        else:
            new_radar.latitude['data'] = radar1.latitude['data']
            new_radar.longitude['data'] = radar1.longitude['data']
            new_radar.altitude['data'] = radar1.altitude['data']

    else:
        new_radar.latitude['data'] = np.append(radar1.latitude['data'],
                                               radar2.latitude['data'])
        new_radar.longitude['data'] = np.append(radar1.longitude['data'],
                                                radar2.longitude['data'])
        new_radar.altitude['data'] = np.append(radar1.altitude['data'],
                                               radar2.altitude['data'])
    return new_radar
