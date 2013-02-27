"""
netcdf

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
    
    """
    ncobj = netCDF4.Dataset(filename)
    
    if 'ray_start_index' in ncobj.variables.keys():
        return _read_netcdf_streamcf(ncobj)
    else:
        return _read_netcdf_cf(ncobj)


def _read_netcdf_cf(ncobj):
    """
    Read a 
    
    Parameters
    ----------
    ncobj : Dataset
        NetCDF4 Dataset to read data from

    Returns
    -------
    radar : Radar
        Radar object

    """

    try:
        mode = "".join(ncobj.variables['sweep_mode'][0])
    except TypeError:
        mode = "".join(ncobj.variables['sweep_mode'][0].data)
    print mode, "azimuth_surveillance    "
    
    if "sur" in mode:
        #ppi
        nsweeps = len(ncobj.variables['sweep_start_ray_index'])
        metadata = dict([(key, getattr(ncobj, key)) for key in
                             ncobj.ncattrs()])
        scan_type = "ppi"
        sweep_mode = np.array(['ppi']*nsweeps)
        if len(ncobj.variables['sweep_start_ray_index']) == 1:
            naz = ncobj.variables['sweep_end_ray_index'][0] + 1
        else:
            naz = (ncobj.variables['sweep_start_ray_index'][1] -
                        ncobj.variables['sweep_start_ray_index'][0])
        nele = ncobj.variables['sweep_start_ray_index'].shape[0]
        ngates = len(ncobj.dimensions['range'])
        loc_dict = {}
        for loc_data in ['latitude', 'altitude', 'longitude']:
            loc_dict.update(
                {loc_data: ncvar_to_field(ncobj.variables[loc_data])})
        location = loc_dict
        sweep_dict = {}
        for sweep_data in ['sweep_start_ray_index', 'sweep_mode', 'sweep_number', 'sweep_end_ray_index', 'fixed_angle']:
            sweep_dict.update(
                {sweep_data: ncvar_to_field(ncobj.variables[sweep_data])})
        sweep_info = sweep_dict
        inst_dict = {}
        for inst_data in ['frequency', 'follow_mode', 'pulse_width', 'prt_mode', 'prt', 'prt_ratio', 'polarization_mode', 'nyquist_velocity', 'unambiguous_range', 'n_samples']:
            if inst_data in ncobj.variables.keys():
                inst_dict.update(
                    {inst_data: (
                        ncvar_to_field(ncobj.variables[inst_data]))})
        inst_params = inst_dict
        azimuth = ncvar_to_field(ncobj.variables['azimuth'])
        _range = ncvar_to_field(ncobj.variables['range'])
        elevation = ncvar_to_field(ncobj.variables['elevation'])
        time = ncvar_to_field(ncobj.variables['time'])
        data_fields = create_field_list(
            ncobj.variables, len(ncobj.dimensions['time']), ngates)
        field_dict = {}
        for field in data_fields:
            print field
            my_field = ncvar_to_field(ncobj.variables[field])
            field_dict.update({field: my_field})
        fields = field_dict
    
    if "sec" in mode:
        #sec
        nsweeps = len(ncobj.variables['sweep_start_ray_index'])
        metadata = dict([(key, getattr(ncobj, key)) for key in
                             ncobj.ncattrs()])
        scan_type = "sec"
        sweep_mode = np.array(['sec'] * nsweeps)
        if len(ncobj.variables['sweep_start_ray_index']) == 1:
            naz = ncobj.variables['sweep_end_ray_index'][0] + 1
        else:
            naz = (ncobj.variables['sweep_start_ray_index'][1] -
                        ncobj.variables['sweep_start_ray_index'][0])
        nele = ncobj.variables['sweep_start_ray_index'].shape[0]
        ngates = len(ncobj.dimensions['range'])
        loc_dict = {}
        for loc_data in ['latitude', 'altitude', 'longitude']:
            loc_dict.update(
                {loc_data: ncvar_to_field(ncobj.variables[loc_data])})
        location = loc_dict
        sweep_dict = {}
        for sweep_data in ['sweep_start_ray_index', 'sweep_mode', 'sweep_number', 'sweep_end_ray_index', 'fixed_angle']:
            sweep_dict.update(
                {sweep_data: ncvar_to_field(ncobj.variables[sweep_data])})
        sweep_info = sweep_dict
        inst_dict = {}
        for inst_data in ['frequency', 'follow_mode', 'pulse_width', 'prt_mode', 'prt', 'prt_ratio', 'polarization_mode', 'nyquist_velocity', 'unambiguous_range', 'n_samples']:
            if inst_data in ncobj.variables.keys():
                inst_dict.update(
                    {inst_data: (
                        ncvar_to_field(ncobj.variables[inst_data]))})
        inst_params = inst_dict
        azimuth = ncvar_to_field(ncobj.variables['azimuth'])
        _range = ncvar_to_field(ncobj.variables['range'])
        elevation = ncvar_to_field(ncobj.variables['elevation'])
        time = ncvar_to_field(ncobj.variables['time'])
        data_fields = create_field_list(
            ncobj.variables, len(ncobj.dimensions['time']), ngates)
        field_dict = {}
        for field in data_fields:
            print field
            my_field = ncvar_to_field(ncobj.variables[field])
            field_dict.update({field: my_field})
        fields = field_dict
    
    if "rhi" in mode:
        #rhi
        metadata = dict([(key, getattr(ncobj, key))
                             for key in ncobj.ncattrs()])
        scan_type = "rhi"
        nsweeps = len(ncobj.variables['sweep_start_ray_index'])
        sweep_mode = np.array(['rhi'] * nsweeps)
        if len(ncobj.variables['sweep_start_ray_index']) == 1:
            nele = ncobj.variables['sweep_end_ray_index'][0] + 1
        else:
            nele = (ncobj.variables['sweep_start_ray_index'][1] -
                         ncobj.variables['sweep_start_ray_index'][0])
        naz = ncobj.variables['sweep_start_ray_index'].shape[0]
        ngates = len(ncobj.dimensions['range'])
        loc_dict = {}
        for loc_data in ['latitude', 'altitude', 'longitude']:
            loc_dict.update(
                {loc_data: ncvar_to_field(ncobj.variables[loc_data])})
        location = loc_dict
        sweep_dict = {}
        for sweep_data in ['sweep_start_ray_index', 'sweep_mode', 'sweep_number', 'sweep_end_ray_index', 'fixed_angle']:
            sweep_dict.update(
                {sweep_data: ncvar_to_field(ncobj.variables[sweep_data])})
        sweep_info = sweep_dict
        inst_dict = {}
        for inst_data in ['frequency', 'follow_mode', 'pulse_width', 'prt_mode', 'prt', 'prt_ratio', 'polarization_mode', 'nyquist_velocity', 'unambiguous_range', 'n_samples']:
            if inst_data in ncobj.variables.keys():
                inst_dict.update(
                    {inst_data: (
                        ncvar_to_field(ncobj.variables[inst_data]))})
        inst_params = inst_dict
        azimuth = ncvar_to_field(ncobj.variables['azimuth'])
        _range = ncvar_to_field(ncobj.variables['range'])
        elevation = ncvar_to_field(ncobj.variables['elevation'])
        time = ncvar_to_field(ncobj.variables['time'])
        data_fields = create_field_list(
            ncobj.variables, len(ncobj.dimensions['time']), ngates)
        field_dict = {}
        for field in data_fields:
            print field
            my_field = ncvar_to_field(ncobj.variables[field])
            field_dict.update({field: my_field})
        fields = field_dict
   
    # XXX
    nrays = 99
    tu = 999
    cal = 999
    sweep_number = 999

    return Radar(nsweeps, nrays, ngates, scan_type, naz, nele, _range,
                 azimuth, elevation, tu, cal, time, fields, sweep_info,
                 sweep_mode, sweep_number, location, inst_params, metadata)

def _read_netcdf_streamcf(ncobj):
    try:
        mode = "".join(ncobj.variables['sweep_mode'][0])
    except TypeError:
        mode = "".join(ncobj.variables['sweep_mode'][0].data)
    print mode, "azimuth_surveillance    "
    
    if mode in "azimuth_surveillance    ":
        #ppi
        print "hi"
        metadata = dict(
            [(key, getattr(ncobj, key)) for key in ncobj.ncattrs()])
        scan_type = "ppi"
        naz = (ncobj.variables['sweep_start_ray_index'][1] -
                    ncobj.variables['sweep_start_ray_index'][0])
        nele = ncobj.variables['sweep_start_ray_index'].shape[0]
        ngates = ncobj.variables['range'].shape[0]
        loc_dict = {}
        for loc_data in ['latitude', 'altitude', 'longitude']:
            loc_dict.update(
                {loc_data: ncvar_to_field(ncobj.variables[loc_data])})
        location = loc_dict
        sweep_dict = {}
        for sweep_data in ['sweep_start_ray_index', 'sweep_mode', 'sweep_number', 'sweep_end_ray_index', 'fixed_angle']:
            sweep_dict.update(
                {sweep_data: ncvar_to_field(ncobj.variables[sweep_data])})
        sweep_info = sweep_dict
        azimuth = ncvar_to_field(ncobj.variables['azimuth'])
        _range = ncvar_to_field(ncobj.variables['range'])
        elevation = ncvar_to_field(ncobj.variables['elevation'])
        time = ncvar_to_field(ncobj.variables['time'])
        data_fields = create_field_list_stream(
            ncobj.variables, (ncobj.variables['ray_start_index'][-1] +
                              ncobj.variables['ray_n_gates'][-1]))
        field_dict = {}
        for field in data_fields:
            print field
            my_field = stream_ncvar_to_field(
                ncobj.variables[field],
                ncobj.variables['sweep_start_ray_index'][:],
                ncobj.variables['sweep_end_ray_index'][:],
                ncobj.variables['ray_n_gates'][:],
                ncobj.variables['range'].shape[0],
                ncobj.variables['time'].shape[0],
                ncobj.variables['ray_start_index'][:])
            field_dict.update({field: my_field})
        fields = field_dict

    # XXX
    nrays = 99
    tu = 999
    cal = 999
    sweep_number = 999
    nsweeps = 999
    sweep_mode = 999
    inst_params = {}

    return Radar(nsweeps, nrays, ngates, scan_type, naz, nele, _range,
                 azimuth, elevation, tu, cal, time, fields, sweep_info,
                 sweep_mode, sweep_number, location, inst_params, metadata)


def create_field_list(variables, nrays, ngates):
    print nrays, ngates
    valid_list = []
    for var in variables.keys():
        if variables[var].shape == (nrays, ngates):
            valid_list.append(var)
    return valid_list


def create_field_list_stream(variables, ngates):
    print ngates
    valid_list = []
    for var in variables.keys():
        if variables[var].shape == (ngates,):
            valid_list.append(var)
    return valid_list


def ncvar_to_field(ncvar):
    outdict = {'data': ncvar[:]}
    outdict.update(dict([(key, getattr(ncvar, key)) for key in
                         ncvar.ncattrs()]))
    return outdict


def stream_ncvar_to_field(ncvar, sweeps, sweepe, ray_len, maxgates, nrays,
                          ray_start_index):
    outdict = {'data': stream_to_2d(
        ncvar[:], sweeps, sweepe, ray_len, maxgates, nrays, ray_start_index)}
    outdict.update(dict([(key, getattr(ncvar, key)) for key in
                   ncvar.ncattrs()]))
    return outdict


def stream_to_2d(data, sweeps, sweepe, ray_len, maxgates, nrays,
                 ray_start_index):
    time_range = np.ma.zeros([nrays, maxgates]) - 9999.0
    cp = 0
    for sweep_number in range(len(sweepe)):
        ss = sweeps[sweep_number]
        se = sweepe[sweep_number]
        rls = ray_len[sweeps[sweep_number]]
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
