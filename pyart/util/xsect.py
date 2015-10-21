"""
pyart.util.xsect
================

Function for extracting cross sections from radar volumes.

.. autosummary::
    :toctree: generated/

    cross_section_ppi
    _copy_dic

"""

from copy import copy

import numpy as np

from ..core import Radar


def cross_section_ppi(radar, target_azimuths):
    """
    Extract cross sections from a PPI volume along one or more azimuth angles.

    Parameters
    ----------
    radar : Radar
        Radar volume containing PPI sweeps from which azimuthal
        cross sections will be extracted.
    target_azimuth : list
        Azimuthal angles in degrees where cross sections will be taken.

    Returns
    -------
    radar_rhi : Radar
        Radar volume containing RHI sweeps which contain azimuthal
        cross sections from the original PPI volume.

    """

    # determine which rays from the ppi radar make up the pseudo RHI
    prhi_rays = []
    rhi_nsweeps = len(target_azimuths)
    ppi_nsweeps = radar.nsweeps

    for target_azimuth in target_azimuths:
        for sweep_slice in radar.iter_slice():
            sweep_azimuths = radar.azimuth['data'][sweep_slice]
            ray_number = np.argmin(np.abs(sweep_azimuths - target_azimuth))
            prhi_rays.append(ray_number + sweep_slice.start)

    _range = _copy_dic(radar.range)
    latitude = _copy_dic(radar.latitude)
    longitude = _copy_dic(radar.longitude)
    altitude = _copy_dic(radar.altitude)
    metadata = _copy_dic(radar.metadata)
    scan_type = 'rhi'

    time = _copy_dic(radar.time, excluded_keys=['data'])
    time['data'] = radar.time['data'][prhi_rays].copy()

    azimuth = _copy_dic(radar.azimuth, excluded_keys=['data'])
    azimuth['data'] = radar.azimuth['data'][prhi_rays].copy()

    elevation = _copy_dic(radar.elevation, excluded_keys=['data'])
    elevation['data'] = radar.elevation['data'][prhi_rays].copy()

    fields = {}
    for field_name, orig_field_dic in radar.fields.items():
        field_dic = _copy_dic(orig_field_dic, excluded_keys=['data'])
        field_dic['data'] = orig_field_dic['data'][prhi_rays].copy()
        fields[field_name] = field_dic

    sweep_number = _copy_dic(radar.sweep_number, excluded_keys=['data'])
    sweep_number['data'] = np.arange(rhi_nsweeps, dtype='int32')

    sweep_mode = _copy_dic(radar.sweep_mode, excluded_keys=['data'])
    sweep_mode['data'] = np.array(['rhi']*rhi_nsweeps)

    fixed_angle = _copy_dic(radar.fixed_angle, excluded_keys=['data'])
    fixed_angle['data'] = np.array(target_azimuths, dtype='float32')

    sweep_start_ray_index = _copy_dic(
        radar.sweep_start_ray_index, excluded_keys=['data'])
    ssri = np.arange(rhi_nsweeps, dtype='int32') * ppi_nsweeps
    sweep_start_ray_index['data'] = ssri

    sweep_end_ray_index = _copy_dic(
        radar.sweep_end_ray_index, excluded_keys=['data'])
    seri = np.arange(rhi_nsweeps, dtype='int32')*ppi_nsweeps + ppi_nsweeps-1
    sweep_end_ray_index['data'] = seri

    radar_rhi = Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle,
        sweep_start_ray_index, sweep_end_ray_index,
        azimuth, elevation)

    return radar_rhi


def _copy_dic(orig_dic, excluded_keys=None):
    """ Return a copy of the original dictionary copying each element. """
    if excluded_keys is None:
        excluded_keys = []
    dic = {}
    for k, v in orig_dic.items():
        if k not in excluded_keys:
            dic[k] = copy(v)
    return dic
