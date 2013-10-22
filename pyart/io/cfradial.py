"""
pyart.io.cfradial
=================

Utilities for reading CF/Radial files.

.. autosummary::
    :toctree: generated/

    read_cfradial
    write_cfradial
    _find_all_meta_group_vars
    _ncvar_to_dict
    _stream_ncvar_to_dict
    _stream_to_2d
    _create_ncvar

"""

import getpass
import datetime
import platform

import numpy as np
import netCDF4

from .common import stringarray_to_chararray
from .radar import Radar


def read_cfradial(filename):
    """
    Read a CF/Radial netCDF file.

    Parameters
    ----------
    filename : str
        Name of CF/Radial netCDF file to read data from.

    Returns
    -------
    radar : Radar
        Radar object.

    Notes
    -----
    This function has not been tested on "stream" CF/Radial files.

    """
    ncobj = netCDF4.Dataset(filename)
    ncvars = ncobj.variables

    # 4.1 Global attribute -> move to metadata dictionary
    metadata = dict([(k, getattr(ncobj, k)) for k in ncobj.ncattrs()])

    # 4.2 Dimensions (do nothing) TODO check if n_points present

    # 4.3 Global variable -> move to metadata dictionary
    if 'volume_number' in ncvars:
        metadata['volume_number'] = int(ncvars['volume_number'][:])
    else:
        metadata['volume_number'] = 0

    global_vars = {'platform_type': 'fixed', 'instrument_type': 'radar',
                   'primary_axis': 'axis_z'}
    # ignore time_* global variables, these are calculated from the time
    # variable when the file is written.
    for var, default_value in global_vars.iteritems():
        if k in ncvars:
            metadata[var] = str(netCDF4.chartostring(ncvars[var][:]))
        else:
            metadata[var] = default_value

    # 4.4 coordinate variables -> create attribute dictionaries
    time = _ncvar_to_dict(ncvars['time'])
    _range = _ncvar_to_dict(ncvars['range'])

    # 4.5 Ray dimension variables TODO working with this

    # 4.6 Location variables -> create attribute dictionaries
    latitude = _ncvar_to_dict(ncvars['latitude'])
    longitude = _ncvar_to_dict(ncvars['longitude'])
    altitude = _ncvar_to_dict(ncvars['altitude'])
    if 'altitude_agl' in ncvars:
        altitude_agl = _ncvar_to_dict(ncvars['altitude_agl'])
    else:
        altitude_agl = None

    # 4.7 Sweep variables -> create atrribute dictionaries
    sweep_number = _ncvar_to_dict(ncvars['sweep_number'])
    sweep_mode = _ncvar_to_dict(ncvars['sweep_mode'])
    fixed_angle = _ncvar_to_dict(ncvars['fixed_angle'])
    sweep_start_ray_index = _ncvar_to_dict(ncvars['sweep_start_ray_index'])
    sweep_end_ray_index = _ncvar_to_dict(ncvars['sweep_end_ray_index'])
    if 'target_scan_rate' in ncvars:
        target_scan_rate = _ncvar_to_dict(ncvars['target_scan_rate'])
    else:
        target_scan_rate = None

    # first sweep mode determines scan_type
    mode = str(netCDF4.chartostring(sweep_mode['data'][0]))
    if "sur" in mode:
        scan_type = "ppi"
    elif "sec" in mode:
        scan_type = "sector"
    elif "rhi" in mode:
        scan_type = "rhi"
    else:
        scan_type = "other"

    # 4.8 Sensor pointing variables -> create attribute dictionaries
    azimuth = _ncvar_to_dict(ncvars['azimuth'])
    elevation = _ncvar_to_dict(ncvars['elevation'])
    if 'scan_rate' in ncvars:
        scan_rate = _ncvar_to_dict(ncvars['scan_rate'])
    else:
        scan_rate = None

    if 'antenna_transition' in ncvars:
        antenna_transition = _ncvar_to_dict(ncvars['antenna_transition'])
    else:
        antenna_transition = None

    # 4.9 Moving platform geo-reference variables
    # TODO moving radar subclass

    # 4.10 Moments field data variables -> field attribute dictionary
    if 'ray_start_index' not in ncvars:     # CF/Radial

        # all variables with dimensions of 'time', 'range' are fields
        fields = dict([(k, _ncvar_to_dict(v)) for k, v in ncvars.iteritems()
                       if v.dimensions == ('time', 'range')])

    else:  # stream file
        ngates = ncvars['ray_start_index'][-1] + ncvars['ray_n_gates'][-1]
        sweeps = ncvars['sweep_start_ray_index'][:]
        sweepe = ncvars['sweep_end_ray_index'][:]
        ray_len = ncvars['ray_n_gates'][:]
        maxgates = ncvars['range'].shape[0]
        nrays = ncvars['time'].shape[0]
        ray_start_index = ncvars['ray_start_index'][:]
        keys = [k for k, v in ncvars.iteritems() if v.shape == (ngates,)]

        fields = {}
        for field in keys:
            fields[field] = _stream_ncvar_to_dict(
                ncvars[field], sweeps, sweepe, ray_len, maxgates, nrays,
                ray_start_index)

    # 4.5 instrument_parameters sub-convention -> instrument_parameters dict

    # the meta_group attribute is often set incorrectly so we cannot
    # use this as a indicator of instrument_parameters
    #keys = _find_all_meta_group_vars(ncvars, 'instrument_parameters')
    valid_keys = ['frequency', 'follow_mode', 'pulse_width', 'prt_mode',
                  'prt', 'prt_ratio', 'polarization_mode', 'nyquist_velocity',
                  'unambiguous_range', 'n_samples', 'sampling_ration']
    keys = [k for k in valid_keys if k in ncvars]
    instrument_parameters = dict((k, _ncvar_to_dict(ncvars[k])) for k in keys)

    # 4.6 radar_parameters sub-convention -> instrument_parameters dict

    # the meta_group attribute is often set incorrectly so we cannot
    # use this as a indicator of instrument_parameters
    #keys = _find_all_meta_group_vars(ncvars, 'radar_parameters')
    valid_keys = ['radar_antenna_gain_h', 'radar_antenna_gain_v',
                  'radar_beam_width_h', 'radar_beam_width_v',
                  'radar_reciever_bandwidth',
                  'radar_measured_transmit_power_h',
                  'radar_measured_transmit_power_v']
    # these keys are not in CF/Radial 1.2 standard but are common
    valid_keys += ['radar_rx_bandwidth', 'measured_transmit_power_h',
                   'measured_transmit_power_v']
    keys = [k for k in valid_keys if k in ncvars]
    radar_parameters = dict((k, _ncvar_to_dict(ncvars[k])) for k in keys)
    instrument_parameters.update(radar_parameters)  # add to instr_params

    if instrument_parameters == {}:  # if no parameters set to None
        instrument_parameters = None

    # 4.7 lidar_parameters sub-convention -> skip

    # 4.8 radar_calibration sub-convention -> radar_calibration
    keys = _find_all_meta_group_vars(ncvars, 'radar_calibration')
    radar_calibration = dict((k, _ncvar_to_dict(ncvars[k])) for k in keys)

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


def _find_all_meta_group_vars(ncvars, meta_group_name):
    """
    Return a list of all variables which are in a given meta_group.
    """
    return [k for k, v in ncvars.iteritems() if 'meta_group' in v.ncattrs()
            and v.meta_group == meta_group_name]


def _ncvar_to_dict(ncvar):
    """ Convert a NetCDF Dataset variable to a dictionary. """
    d = dict((k, getattr(ncvar, k)) for k in ncvar.ncattrs())
    d['data'] = ncvar[:]
    return d


def _stream_ncvar_to_dict(ncvar, sweeps, sweepe, ray_len, maxgates, nrays,
                          ray_start_index):
    """ Convert a Stream NetCDF Dataset variable to a dict. """
    d = dict((k, getattr(ncvar, k)) for k in ncvar.ncattrs())
    data = _stream_to_2d(ncvar[:], sweeps, sweepe, ray_len, maxgates, nrays,
                         ray_start_index)
    d['data'] = data
    return d


def _stream_to_2d(data, sweeps, sweepe, ray_len, maxgates, nrays,
                  ray_start_index):
    """ Convert a 1D stream to a 2D array. """
    # XXX clean this up, need to find sample data
    time_range = np.ma.zeros([nrays, maxgates]) - 9999.0
    cp = 0
    for sweep_number in range(len(sweepe)):
        ss = sweeps[sweep_number]
        se = sweepe[sweep_number]
        rle = ray_len[sweeps[sweep_number]]

        if ray_len[ss:se].sum() == rle * (se - ss):
            time_range[ss:se, 0:rle] = (
                data[cp:cp + (se - ss) * rle].reshape(se - ss, rle))
            cp += (se - ss) * rle
        else:
            for rn in range(se - ss):
                time_range[ss + rn, 0:ray_len[ss + rn]] = (
                    data[ray_start_index[ss + rn]:ray_start_index[ss + rn] +
                         ray_len[ss+rn]])
            cp += ray_len[ss:se].sum()
    return time_range


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
