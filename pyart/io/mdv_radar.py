"""
pyart.io.mdv_radar
==================

Utilities for reading of MDV radar files.

.. autosummary::
    :toctree: generated/

    read_mdv

"""

import numpy as np
from netCDF4 import date2num

from ..config import FileMetadata, get_fillvalue
from ..core.radar import Radar
from .common import make_time_unit_str, _test_arguments, prepare_for_read
from ..lazydict import LazyLoadDict
from . import mdv_common


def read_mdv(filename, field_names=None, additional_metadata=None,
             file_field_names=False, exclude_fields=None,
             delay_field_loading=False, **kwargs):
    """
    Read a MDV file.

    Parameters
    ----------
    filename : str
        Name of MDV file to read or file-like object pointing to the
        beginning of such a file.
    field_names : dict, optional
        Dictionary mapping MDV data type names to radar field names. If a
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
        True to use the MDV data type names for the field names. If this
        case the field_names parameter is ignored. The field dictionary will
        likely only have a 'data' key, unless the fields are defined in
        `additional_metadata`.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters.
    delay_field_loading : bool
        True to delay loading of field data from the file until the 'data'
        key in a particular field dictionary is accessed.  In this case
        the field attribute of the returned Radar object will contain
        LazyLoadDict objects not dict objects. Not all file types support this
        parameter.

    Returns
    -------
    radar : Radar
        Radar object containing data from MDV file.

    Notes
    -----
    Currently this function can only read polar MDV files with fields
    compressed with gzip or zlib.

    """
    # test for non empty kwargs
    _test_arguments(kwargs)

    # create metadata retrieval object
    filemetadata = FileMetadata('mdv', field_names, additional_metadata,
                                file_field_names, exclude_fields)

    mdvfile = mdv_common.MdvFile(prepare_for_read(filename))

    # value attributes
    az_deg, range_km, el_deg = mdvfile._calc_geometry()
    naz = len(az_deg)
    nele = len(el_deg)
    scan_type = mdvfile.projection

    if scan_type not in ['ppi', 'rhi']:
        raise NotImplementedError('No support for scan_type %s.' % scan_type)

    # time
    time = filemetadata('time')
    units = make_time_unit_str(mdvfile.times['time_begin'])
    time['units'] = units
    time_start = date2num(mdvfile.times['time_begin'], units)
    time_end = date2num(mdvfile.times['time_end'], units)
    time['data'] = np.linspace(time_start, time_end, naz * nele)

    # range
    _range = filemetadata('range')
    _range['data'] = np.array(range_km * 1000.0, dtype='float32')
    _range['meters_to_center_of_first_gate'] = _range['data'][0]
    _range['meters_between_gates'] = (_range['data'][1] - _range['data'][0])

    # fields
    fields = {}
    for mdv_field in set(mdvfile.fields):
        field_name = filemetadata.get_field_name(mdv_field)
        if field_name is None:
            continue

        # create and store the field dictionary
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue()
        dataextractor = mdv_common._MdvVolumeDataExtractor(
            mdvfile, mdvfile.fields.index(mdv_field), get_fillvalue())
        if delay_field_loading:
            field_dic = LazyLoadDict(field_dic)
            field_dic.set_lazy('data', dataextractor)
        else:
            field_dic['data'] = dataextractor()
        fields[field_name] = field_dic

    # metadata
    metadata = filemetadata('metadata')
    for meta_key, mdv_key in mdv_common.MDV_METADATA_MAP.items():
        metadata[meta_key] = mdvfile.master_header[mdv_key]

    # latitude
    latitude = filemetadata('latitude')
    latitude['data'] = np.array([mdvfile.radar_info['latitude_deg']],
                                dtype='float64')
    # longitude
    longitude = filemetadata('longitude')
    longitude['data'] = np.array([mdvfile.radar_info['longitude_deg']],
                                 dtype='float64')
    # altitude
    altitude = filemetadata('altitude')
    altitude['data'] = np.array([mdvfile.radar_info['altitude_km'] * 1000.0],
                                dtype='float64')

    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
    # sweep_end_ray_index
    sweep_number = filemetadata('sweep_number')
    sweep_mode = filemetadata('sweep_mode')
    fixed_angle = filemetadata('fixed_angle')
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')
    len_time = len(time['data'])

    if scan_type == 'ppi':
        nsweeps = nele
        sweep_number['data'] = np.arange(nsweeps, dtype='int32')
        sweep_mode['data'] = np.array(
            nsweeps * ['azimuth_surveillance'], dtype='S')
        fixed_angle['data'] = np.array(el_deg, dtype='float32')
        sweep_start_ray_index['data'] = np.arange(0, len_time, naz,
                                                  dtype='int32')
        sweep_end_ray_index['data'] = np.arange(naz - 1, len_time, naz,
                                                dtype='int32')

    elif scan_type == 'rhi':
        nsweeps = naz
        sweep_number['data'] = np.arange(nsweeps, dtype='int32')
        sweep_mode['data'] = np.array(nsweeps * ['rhi'], dtype='S')
        fixed_angle['data'] = np.array(az_deg, dtype='float32')
        sweep_start_ray_index['data'] = np.arange(0, len_time, nele,
                                                  dtype='int32')
        sweep_end_ray_index['data'] = np.arange(nele - 1, len_time, nele,
                                                dtype='int32')

    # azimuth, elevation
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')

    if scan_type == 'ppi':
        azimuth['data'] = np.tile(az_deg, nele)
        elevation['data'] = np.array(el_deg).repeat(naz)

    elif scan_type == 'rhi':
        azimuth['data'] = np.array(az_deg).repeat(nele)
        elevation['data'] = np.tile(el_deg, naz)

    # instrument parameters
    # we will set 4 parameters in the instrument_parameters dict
    # prt, prt_mode, unambiguous_range, and nyquist_velocity

    # TODO prt mode: Need to fix this.. assumes dual if two prts
    if mdvfile.radar_info['prt2_s'] == 0.0:
        prt_mode_str = 'fixed'
    else:
        prt_mode_str = 'dual'

    prt_mode = filemetadata('prt_mode')
    prt = filemetadata('prt')
    unambiguous_range = filemetadata('unambiguous_range')
    nyquist_velocity = filemetadata('nyquist_velocity')
    beam_width_h = filemetadata('radar_beam_width_h')
    beam_width_v = filemetadata('radar_beam_width_v')

    prt_mode['data'] = np.array([prt_mode_str] * nsweeps, dtype='S')
    prt['data'] = np.array([mdvfile.radar_info['prt_s']] * nele * naz,
                           dtype='float32')

    urange_m = mdvfile.radar_info['unambig_range_km'] * 1000.0
    unambiguous_range['data'] = np.array([urange_m] * naz * nele,
                                         dtype='float32')

    uvel_mps = mdvfile.radar_info['unambig_vel_mps']
    nyquist_velocity['data'] = np.array([uvel_mps] * naz * nele,
                                        dtype='float32')
    beam_width_h['data'] = np.array(
        [mdvfile.radar_info['horiz_beam_width_deg']], dtype='float32')
    beam_width_v['data'] = np.array(
        [mdvfile.radar_info['vert_beam_width_deg']], dtype='float32')

    instrument_parameters = {'prt_mode': prt_mode, 'prt': prt,
                             'unambiguous_range': unambiguous_range,
                             'nyquist_velocity': nyquist_velocity,
                             'radar_beam_width_h': beam_width_h,
                             'radar_beam_width_v': beam_width_v}

    if not delay_field_loading:
        mdvfile.close()
    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)
