"""
pyart.io.read_gamic
=================

Utilities for reading gamic hdf5 files.

.. autosummary::
    :toctree: generated/

    read_gamic
    write_cfradial
    _h5_to_dict
    _h5_moment_to_dict
    _create_ncvar

"""

import getpass
import datetime
import platform

import h5py
import numpy as np
import netCDF4

from ..config import FileMetadata
from .common import stringarray_to_chararray, make_time_unit_str
from .radar import Radar


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
            h5obj['/scan0/how'].attrs['timestamp'], '%Y-%m-%dT%H:%M:%S.000Z')
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


def write_cfradial(filename, radar, format='NETCDF4', time_reference=False):
    """
    Write a Radar object to a CF/Radial compliant netCDF file.

    Parameters
    ----------
    filename : str
        Filename to create.
    radar : Radar
        Radar object.
    format : str, optional
        NetCDF format, one of 'NETCDF4', 'NETCDF4_CLASSIC',
        'NETCDF3_CLASSIC' or 'NETCDF3_64BIT'. See netCDF4 documentation for
        details.
    time_reference : bool
        True to include a time_reference variable, False will not include
        this variable.

    """
    dataset = netCDF4.Dataset(filename, 'w', format=format)

    # determine the maximum string length
    max_str_len = len(radar.sweep_mode['data'][0])
    for k in ['follow_mode', 'prt_mode', 'polarization_mode']:
        if k in radar.instrument_parameters:
            sdim_length = len(radar.instrument_parameters[k]['data'][0])
            max_str_len = max(max_str_len, sdim_length)
    str_len = max(max_str_len, 32)      # minimum string legth of 32

    # create time, range and sweep dimensions
    dataset.createDimension('time', None)
    dataset.createDimension('range', radar.ngates)
    dataset.createDimension('sweep', radar.nsweeps)
    dataset.createDimension('string_length', str_len)

    # global attributes
    # remove global variables from copy of metadata
    metadata_copy = dict(radar.metadata)
    global_variables = ['volume_number', 'platform_type', 'instrument_type',
                        'primary_axis', 'time_coverage_start',
                        'time_coverage_end', 'time_reference']
    for var in global_variables:
        if var in metadata_copy:
            metadata_copy.pop(var)

    # determine the history attribute if it doesn't exist, save for
    # the last attribute.
    if 'history' in metadata_copy:
        history = metadata_copy.pop('history')
    else:
        user = getpass.getuser()
        node = platform.node()
        time_str = datetime.datetime.now().isoformat()
        t = (user, node, time_str)
        history = 'created by %s on %s at %s using Py-ART' % (t)

    dataset.setncatts(metadata_copy)

    if 'Conventions' not in dataset.ncattrs():
        dataset.setncattr('Conventions', "CF/Radial")

    if 'field_names' not in dataset.ncattrs():
        dataset.setncattr('field_names', ', '.join(radar.fields.keys()))

    # history should be the last attribute, ARM standard
    dataset.setncattr('history',  history)

    # standard variables
    _create_ncvar(radar.time, dataset, 'time', ('time', ))
    _create_ncvar(radar.range, dataset, 'range', ('range', ))
    _create_ncvar(radar.azimuth, dataset, 'azimuth', ('time', ))
    _create_ncvar(radar.elevation, dataset, 'elevation', ('time', ))

    # fields
    for field, dic in radar.fields.iteritems():
        _create_ncvar(dic, dataset, field, ('time', 'range'))

    # sweep parameters
    _create_ncvar(radar.sweep_number, dataset, 'sweep_number', ('sweep', ))
    _create_ncvar(radar.fixed_angle, dataset, 'fixed_angle', ('sweep', ))
    _create_ncvar(radar.sweep_start_ray_index, dataset,
                  'sweep_start_ray_index', ('sweep', ))
    _create_ncvar(radar.sweep_end_ray_index, dataset,
                  'sweep_end_ray_index', ('sweep', ))
    _create_ncvar(radar.sweep_mode, dataset, 'sweep_mode',
                  ('sweep', 'string_length'))

    # instrument_parameters
    if 'frequency' in radar.instrument_parameters.keys():
        size = len(radar.instrument_parameters['frequency']['data'])
        dataset.createDimension('frequency', size)

    instrument_dimensions = {
        'frequency': ('frequency'),
        'follow_mode': ('sweep', 'string_length'),
        'pulse_width': ('time', ),
        'prt_mode': ('sweep', 'string_length'),
        'prt': ('time', ),
        'prt_ratio': ('time', ),
        'polarization_mode': ('sweep', 'string_length'),
        'nyquist_velocity': ('time', ),
        'unambiguous_range': ('time', ),
        'n_samples': ('time', ),
        'sampling_ratio': ('time', ),
        'radar_antenna_gain_h': (),
        'radar_antenna_gain_v': (),
        'radar_beam_width_h': (),
        'radar_beam_width_v': (),
        'radar_reciever_bandwidth': (),
        'radar_rx_bandwidth': (),           # non-standard
        'radar_measured_transmit_power_h': ('time', ),
        'radar_measured_transmit_power_v': ('time', ),
        'measured_transmit_power_v': ('time', ),    # non-standard
        'measured_transmit_power_h': ('time', ),    # non-standard
    }
    for k in radar.instrument_parameters.keys():
        if k in instrument_dimensions:
            dim = instrument_dimensions[k]
        else:
            dim = ()
        _create_ncvar(radar.instrument_parameters[k], dataset, k, dim)

    # latitude, longitude, altitude
    # TODO moving platform
    _create_ncvar(radar.latitude, dataset, 'latitude', ())
    _create_ncvar(radar.longitude, dataset, 'longitude', ())
    _create_ncvar(radar.altitude, dataset, 'altitude', ())

    # time_coverage_start and time_coverage_end variables
    time_dim = ('string_length', )
    units = radar.time['units']
    start_dt = netCDF4.num2date(radar.time['data'][0], units)
    end_dt = netCDF4.num2date(radar.time['data'][-1], units)
    start_dic = {'data': np.array(start_dt.isoformat() + 'Z'),
                 'long_name': 'UTC time of first ray in the file',
                 'units': 'unitless'}
    end_dic = {'data': np.array(end_dt.isoformat() + 'Z'),
               'long_name': 'UTC time of last ray in the file',
               'units': 'unitless'}
    _create_ncvar(start_dic, dataset, 'time_coverage_start', time_dim)
    _create_ncvar(end_dic, dataset, 'time_coverage_end', time_dim)
    if time_reference:
        ref_dic = {'data': np.array(radar.time['units'][-20:]),
                   'long_name': 'UTC time reference',
                   'units': 'unitless'}
        _create_ncvar(ref_dic, dataset, 'time_reference', time_dim)

    vol_dic = {'data': np.array([0], dtype='int32'),
               'long_name': 'Volume number',
               'units': 'unitless'}
    _create_ncvar(vol_dic, dataset, 'volume_number', ())
    dataset.close()


def _create_ncvar(dic, dataset, name, dimensions):
    """
    Create and fill a Variable in a netCDF Dataset object.

    Parameters
    ----------
    dic : dict
        Radar dictionary to containing variable data and meta-data
    dataset : Dataset
        NetCDF dataset to create variable in.
    name : str
        Name of variable to create.
    dimension : tuple of str
        Dimension of variable.

    """
    # create array from list, etc.
    data = dic['data']
    if isinstance(data, np.ndarray) is not True:
        print "Warning, converting non-array to array:", name
        data = np.array(data)

    # convert string array to character arrays
    if data.dtype.char is 'S' and data.dtype != 'S1':
        data = stringarray_to_chararray(data)

    # create the dataset variable
    if 'least_significant_digit' in dic:
        lsd = dic['least_significant_digit']
    else:
        lsd = None
    if "_FillValue" in dic:
        fill_value = dic['_FillValue']
    else:
        fill_value = None

    ncvar = dataset.createVariable(name, data.dtype, dimensions,
                                   zlib=True, least_significant_digit=lsd,
                                   fill_value=fill_value)

    # long_name attribute first if present, ARM standard
    if 'long_name' in dic.keys():
        ncvar.setncattr('long_name', dic['long_name'])

    # units attribute second if present, ARM standard
    if 'units' in dic.keys():
        ncvar.setncattr('units', dic['units'])

    # remove _FillValue and replace to make it the third attribute.
    if '_FillValue' in ncvar.ncattrs():
        fv = ncvar._FillValue
        ncvar.delncattr('_FillValue')
        ncvar.setncattr('_FillValue', fv)

    # set all attributes
    for key, value in dic.iteritems():
        if key not in ['data', '_FillValue', 'long_name', 'units']:
            ncvar.setncattr(key, value)

    # set the data
    if data.shape == ():
        data.shape = (1,)
    if data.dtype == 'S1':  # string/char arrays
        ncvar[..., :data.shape[-1]] = data[:]
    else:
        ncvar[:] = data[:]
