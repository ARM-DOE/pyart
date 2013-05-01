"""
pyart.io.netcdf
===============

Utilities for reading netcdf files.

.. autosummary::
    :toctree: generated/

    read_netcdf
    write_netcdf
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

from radar import Radar


def read_netcdf(filename):
    """
    Read a netCDF file.

    Parameters
    ----------
    filename : str
        Name of netCDF file to read data from.

    Returns
    -------
    radar : Radar
        Radar object.

    Notes
    -----
    This function has not been tested on "stream" netCDF files.

    """
    ncobj = netCDF4.Dataset(filename)
    ncvars = ncobj.variables

    # azimuth, range, elevation, time, nsweeps, ngates
    azimuth = _ncvar_to_dict(ncvars['azimuth'])
    _range = _ncvar_to_dict(ncvars['range'])
    elevation = _ncvar_to_dict(ncvars['elevation'])
    time = _ncvar_to_dict(ncvars['time'])
    nsweeps = len(ncvars['sweep_start_ray_index'])
    ngates = len(ncobj.dimensions['range'])

    # sweep info
    keys = ['sweep_start_ray_index', 'sweep_mode', 'sweep_number',
            'sweep_end_ray_index', 'fixed_angle']
    sweep_info = dict((k, _ncvar_to_dict(ncvars[k])) for k in keys)

    # location
    keys = ['latitude', 'altitude', 'longitude']
    location = dict((k, _ncvar_to_dict(ncvars[k])) for k in keys)

    # inst_params
    keys = ['frequency', 'follow_mode', 'pulse_width', 'prt_mode', 'prt',
            'prt_ratio', 'polarization_mode', 'nyquist_velocity',
            'unambiguous_range', 'n_samples']
    keys = [k for k in keys if k in ncvars.keys()]  # only those present
    inst_params = dict((k, _ncvar_to_dict(ncvars[k])) for k in keys)

    # metadata
    metadata = dict([(k, getattr(ncobj, k)) for k in ncobj.ncattrs()])

    # naz, nele
    ssri = ncvars['sweep_start_ray_index']
    if len(ssri) == 1:
        naz = ssri[0] + 1
    else:
        naz = ssri[1] - ssri[0]
    nele = ssri.shape[0]

    try:
        mode = "".join(ncvars['sweep_mode'][0])
    except TypeError:
        mode = "".join(ncvars['sweep_mode'][0].data)

    if "sur" in mode:
        scan_type = "ppi"
    elif "sec" in mode:
        scan_type = "sec"
    elif "rhi" in mode:
        scan_type = "rhi"
        nele, naz = nele, naz

    sweep_mode = np.array([scan_type] * nsweeps)

    # fields and nrays
    if 'ray_start_index' in ncvars.keys():
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
    else:
        # CF/Radial
        nrays = len(ncobj.dimensions['time'])
        shape = (nrays, ngates)
        keys = [k for k, v in ncvars.iteritems() if v.shape == shape]
        fields = dict((k, _ncvar_to_dict(ncvars[k])) for k in keys)

    # XXX
    tu = 999
    cal = 999
    sweep_number = 999

    return Radar(nsweeps, nrays, ngates, scan_type, naz, nele, _range,
                 azimuth, elevation, tu, cal, time, fields, sweep_info,
                 sweep_mode, sweep_number, location, inst_params, metadata)


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


def write_netcdf(filename, radar, format='NETCDF4'):
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

    """
    dataset = netCDF4.Dataset(filename, 'w', format=format)

    # create time, range and sweep dimensions
    dataset.createDimension('time', None)
    dataset.createDimension('range', radar.ngates)
    dataset.createDimension('sweep', radar.nsweeps)
    dataset.createDimension('string_length_short', 20)

    # global attributes
    dataset.setncatts(radar.metadata)
    if 'history' not in dataset.ncattrs():
        user = getpass.getuser()
        node = platform.node()
        time_str = datetime.datetime.now().isoformat()
        t = (user, node, time_str)
        history_str = 'created by %s on %s at %s using Py-ART' % (t)
        dataset.setncattr('history',  history_str)
    if 'Conventions' not in dataset.ncattrs():
        dataset.setncattr('Conventions', "CF/Radial")

    # standard variables
    _create_ncvar(radar.time, dataset, 'time', ('time', ))
    _create_ncvar(radar.range, dataset, 'range', ('range', ))
    _create_ncvar(radar.azimuth, dataset, 'azimuth', ('time', ))
    _create_ncvar(radar.elevation, dataset, 'elevation', ('time', ))

    # fields
    for field, dic in radar.fields.iteritems():
        _create_ncvar(dic, dataset, field, ('time', 'range'))

    # sweep parameters
    keys = ['fixed_angle', 'sweep_start_ray_index', 'sweep_end_ray_index',
            'sweep_number']
    for key in keys:
        _create_ncvar(radar.sweep_info[key], dataset, key, ('sweep', ))

    sdim_length = len(radar.sweep_info['sweep_mode']['data'][0])
    sdim_string = 'string_length_%d' % (sdim_length)
    dataset.createDimension(sdim_string, sdim_length)
    _create_ncvar(radar.sweep_info['sweep_mode'], dataset, 'sweep_mode',
                  ('sweep', sdim_string))

    # inst_params
    if 'frequency' in radar.inst_params.keys():
        size = len(radar.inst_params['frequency']['data'])
        dataset.createDimension('frequency', size)
    inst_dims = {
        'frequency': ('frequency'),
        'follow_mode': ('sweep', sdim_string),
        'pulse_width': ('time', ),
        'prt_mode': ('sweep', sdim_string),
        'prt': ('time', ),
        'prt_ratio': ('time', ),
        'polarization_mode': ('sweep', sdim_string),
        'nyquist_velocity': ('time', ),
        'unambiguous_range': ('time', ),
        'n_samples': ('time', ),
    }
    for k in radar.inst_params.keys():
        _create_ncvar(radar.inst_params[k], dataset, k, inst_dims[k])

    # location parameters
    # TODO moving platform
    for k in radar.location.keys():
        _create_ncvar(radar.location[k], dataset, k, ())

    # time_coverage_start and time_coverage_end variables
    time_dim = ('string_length_short', )
    units = radar.time['units']
    start_dt = netCDF4.num2date(radar.time['data'][0], units)
    end_dt = netCDF4.num2date(radar.time['data'][-1], units)
    start_dic = {'data': np.array(start_dt.isoformat() + 'Z')}
    end_dic = {'data': np.array(end_dt.isoformat() + 'Z')}
    _create_ncvar(start_dic, dataset, 'time_coverage_start', time_dim)
    _create_ncvar(end_dic, dataset, 'time_coverage_end', time_dim)

    _create_ncvar({'data': np.array([0], dtype='int32')}, dataset,
                  'volume_number', ())

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
        data = netCDF4.stringtochar(data)

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

    # set all attributes
    for key, value in dic.iteritems():
        if key not in ['data', '_FillValue']:
            ncvar.setncattr(key, value)

    # set the data
    if data.shape == ():
        data.shape = (1,)
    ncvar[:] = data[:]
    #if type(data) == np.ma.MaskedArray:
    #    ncvar[:] = data.data
    #else:
    #    ncvar[:] = data
