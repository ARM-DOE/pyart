"""
pyart.io.netcdf
===============

Utilities for reading netcdf files.

.. autosummary::
    :toctree: generated/

    read_netcdf
    _ncvar_to_dict
    _stream_ncvar_to_dict
    _stream_to_2d


"""

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
