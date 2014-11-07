"""
pyart.aux_io.read_gamic
=======================

Utilities for reading gamic hdf5 files.

.. autosummary::
    :toctree: generated/

    read_gamic
    _h5_to_dict
    _h5_moments_to_dict

"""

import datetime

import h5py
import numpy as np

from ..config import FileMetadata
from ..io.common import make_time_unit_str
from ..core.radar import Radar


def read_gamic(filename, field_names=None, additional_metadata=None,
               file_field_names=False, exclude_fields=None):
    """
    Read a gamic hdf5 file.

    Parameters
    ----------
    filename : str
        Name of gamic hdf5 file to read data from.

    Returns
    -------
    radar : Radar
        Radar object.

    Notes
    -----
    First Test.

    """
    # create metadata retrieval object
    filemetadata = FileMetadata('cfradial', field_names, additional_metadata,
                                file_field_names, exclude_fields)

    # open h5 file and get handle
    h5obj = h5py.File(filename, 'r')

    # initialize metadata as dict
    metadata = dict()

    # 4.1 Global attribute -> move to metadata dictionary
    # no global attribs in gamic hdf5

    # 4.2 Dimensions (do nothing) TODO check if n_points present

    # 4.3 Global variable -> move to metadata dictionary
    metadata['volume_number'] = 0

    #map some global vars to possible gamic counterparts
    global_vars = {'platform_type': 'fixed', 'instrument_type': 'radar',
                   'primary_axis': 'axis_z', 'source': 'software',
                   'references': 'template_name',
                   'instrument_name': 'site_name',
                   'institution': 'host_name', 'version': 'sdp_name'}
    # ignore time_* global variables, these are calculated from the time
    # variable when the file is written.
    for var, default_value in global_vars.iteritems():
        try:
            metadata[var] = h5obj['/how'].attrs[default_value]
        except KeyError:
            metadata[var] = default_value

    # 4.4 coordinate variables -> create attribute dictionaries
    time = filemetadata.get_metadata('time')

    # test for possible TimeStamp Issues
    try:
        scan_time = datetime.datetime.strptime(
            h5obj['/scan0/how'].attrs['timestamp'], '%Y-%m-%dT%H:%M:%SZ')
    except ValueError:
        scan_time = datetime.datetime.strptime(
            h5obj['/scan0/how'].attrs['timestamp'], '%Y-%m-%dT%H:%M:%S.%fZ')
    time['units'] = make_time_unit_str(scan_time)

    # get scan0 ray header
    # in volume file there are several scans, scanX
    ray_header = h5obj['/scan0/ray_header']

    # get microseconds since linux epoch
    scan_time_se = (scan_time -
                    datetime.datetime(1970, 1, 1)).total_seconds() * 1e6
    # get timestamp array and subtract epoch, change to seconds again
    time['data'] = (np.array(ray_header['timestamp']) - scan_time_se) / 1e6

    # get range, bin etc
    range_start = h5obj['/scan0/how'].attrs['range_start']
    range_ = h5obj['/scan0/how'].attrs['range']
    bin_count = h5obj['/scan0/how'].attrs['bin_count']
    range_step = h5obj['/scan0/how'].attrs['range_step']
    range_samples = h5obj['/scan0/how'].attrs['range_samples']
    _range = filemetadata.get_metadata('range')
    # create range array
    _range['data'] = np.linspace(
        range_start + (range_step * range_samples / 2.),
        range_ - (range_step * range_samples / 2.), bin_count)

    # 4.5 Ray dimension variables TODO working with this

    # 4.6 Location variables -> create attribute dictionaries
    # gamic location variables are in subgroup /where
    latitude = _h5_to_dict(h5obj, '/where', 'lat', 'latitude', filemetadata)
    longitude = _h5_to_dict(h5obj, '/where', 'lon', 'longitude', filemetadata)
    altitude = _h5_to_dict(h5obj, '/where', 'height', 'altitude',
                           filemetadata)
    altitude_agl = None

    ## 4.7 Sweep variables -> create atrribute dictionaries
    # if only one scan -> one sweep
    # TODO: account for volume scans
    sweep_number = filemetadata.get_metadata('sweep_number')
    sweep_number['data'] = np.array([0])
    sweep_mode = _h5_to_dict(h5obj, '/what', 'object', 'sweep_mode',
                             filemetadata)
    if 'PVOL' in sweep_mode['data']:
        fixed_angle = _h5_to_dict(h5obj, '/scan0/how', 'elevation',
                                  'fixed_angle', filemetadata)
    elif 'RHI' in sweep_mode['data']:
        fixed_angle = _h5_to_dict(h5obj, '/scan0/how', 'azimuth',
                                  'fixed_angle', filemetadata)

    sweep_start_ray_index = filemetadata.get_metadata('sweep_start_ray_index')
    sweep_start_ray_index['data'] = np.array([0])
    sweep_end_ray_index = filemetadata.get_metadata('sweep_end_ray_index')
    sweep_end_ray_index['data'] = [h5obj['/scan0/how'].attrs['ray_count'] - 1]

    # target scan speed is used
    target_scan_rate = _h5_to_dict(h5obj, '/scan0/how', 'scan_speed',
                                   'target_scan_rate', filemetadata)

    # first sweep mode determines scan_type
    mode = sweep_mode['data']
    if "PVOL" in mode:
        scan_type = "ppi"
    elif "sec" in mode:
        scan_type = "sector"
    elif "RHI" in mode:
        scan_type = "rhi"
    else:
        scan_type = "other"

    # 4.8 Sensor pointing variables -> create attribute dictionaries
    azi_start = np.array(ray_header['azimuth_start'])
    azi_stop = np.array(ray_header['azimuth_stop'])
    azimuth = filemetadata.get_metadata('azimuth')
    zero_index = np.where(azi_stop < azi_start)
    azi_stop[zero_index[0]] += 360
    azimuth['data'] = (azi_start + azi_stop) / 2.0

    ele_start = np.array(ray_header['elevation_start'])
    ele_stop = np.array(ray_header['elevation_stop'])
    elevation = filemetadata.get_metadata('elevation')
    elevation['data'] = (ele_start + ele_stop) / 2.0

    # scan speed per ray could be read
    scan_rate = None

    # antenna transistions are not recorder in gamic hdf5
    antenna_transition = None

    # 4.9 Moving platform geo-reference variables
    # TODO moving radar subclass

    # 4.10 Moments field data variables -> field attribute dictionary
    # moments are in subgroup /scanX
    moments = h5obj['/scan0']
    fields = dict(
        [(moments[k].attrs['moment'], _h5_moments_to_dict(v, filemetadata))
         for k, v in moments.iteritems() if 'mom' in k])

    ## 4.5 instrument_parameters sub-convention -> instrument_parameters dict

    ## the meta_group attribute is often set incorrectly so we cannot
    ## use this as a indicator of instrument_parameters
    ## need to discuss gamic hdf5 instrument_parameters
    ##keys = _find_all_meta_group_vars(ncvars, 'instrument_parameters')
    #valid_keys = ['frequency', 'follow_mode', 'pulse_width', 'prt_mode',
                  #'prt', 'prt_ratio', 'polarization_mode', 'nyquist_velocity',
                  #'unambiguous_range', 'n_samples', 'sampling_ration']

    #keys = [k for k in valid_keys if k in ncvars]
    #instrument_parameters = dict((k, _ncvar_to_dict(ncvars[k])) for k in keys)

    ## 4.6 radar_parameters sub-convention -> instrument_parameters dict

    ## the meta_group attribute is often set incorrectly so we cannot
    ## use this as a indicator of instrument_parameters
    ##keys = _find_all_meta_group_vars(ncvars, 'radar_parameters')
    #valid_keys = ['radar_antenna_gain_h', 'radar_antenna_gain_v',
                  #'radar_beam_width_h', 'radar_beam_width_v',
                  #'radar_reciever_bandwidth',
                  #'radar_measured_transmit_power_h',
                  #'radar_measured_transmit_power_v']
    ## these keys are not in CF/Radial 1.2 standard but are common
    #valid_keys += ['radar_rx_bandwidth', 'measured_transmit_power_h',
                   #'measured_transmit_power_v']
    #keys = [k for k in valid_keys if k in ncvars]
    #radar_parameters = dict((k, _ncvar_to_dict(ncvars[k])) for k in keys)
    #instrument_parameters.update(radar_parameters)  # add to instr_params

    #if instrument_parameters == {}:  # if no parameters set to None
    instrument_parameters = None

    # 4.7 lidar_parameters sub-convention -> skip

    ## 4.8 radar_calibration sub-convention -> radar_calibration
    #keys = _find_all_meta_group_vars(ncvars, 'radar_calibration')
    #radar_calibration = dict((k, _ncvar_to_dict(ncvars[k])) for k in keys)
    radar_calibration = None

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters,
        radar_calibration=radar_calibration,
        altitude_agl=altitude_agl,
        scan_rate=scan_rate,
        antenna_transition=antenna_transition)


def _h5_to_dict(h5obj, h5path, h5var, ncstring, fmd):
    """ Convert a HDF Attribute variable to a dictionary. """
    d = fmd.get_metadata(ncstring)
    d['data'] = np.array([h5obj[h5path].attrs[h5var]])
    return d


def _h5_moments_to_dict(h5obj, fmd):
    """
    Convert gamic HDF5 moment Dataset and attached Attributes to a dictionary.
    """
    d = fmd.get_metadata('default')
    d['valid_min'] = h5obj.attrs['dyn_range_min']
    d['valid_max'] = h5obj.attrs['dyn_range_max']
    d['standard_name'] = h5obj.attrs['moment']
    d['long_name'] = h5obj.attrs['moment']
    #d['units'] =

    if h5obj.attrs['format'] == 'UV8':
        div = 256.0
    else:
        div = 65536.0

    d['data'] = (d['valid_min'] + h5obj[...] *
                 (d['valid_max'] - d['valid_min']) / div)
    return d
