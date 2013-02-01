"""
Utilities for saving mapped and radar co-ordinate radar data
"""

import sys
import socket
import getpass
import datetime as dt

import netCDF4
import numpy as np


def is_moment(varname, moment_fixes):
    moments = moment_fixes.keys()
    return True in [foo in varname for foo in moments]


def is_radar(varname, radar_list):
    return True in [foo in varname for foo in radar_list]


def fix_variables(cgfile):
    debug = True
    print "go"
    moment_fixes = {
        'DBZ_F': {
            'units': 'dBZ',
            'standard_name': 'equivalent_reflectivity_factor',
            'valid_max': 80.0,
            'valid_min': -45.0},
        'VEL_F': {
            'units': 'm/s',
            'standard_name': (
                'radial_velocity_of_scatterers_away_from_instrument'),
            'valid_max': 45.0,
            'valid_min': -45.0},
        'KDP_F': {
            'units': 'degrees/km',
            'standard_name': 'specific_differential_phase_hv',
            'valid_max': 10.0,
            'valid_min': -10.0},
        'ZDR_F': {
            'units': 'dB',
            'standard_name': 'log_differential_reflectivity_hv',
            'valid_max': 6.0,
            'valid_min': -6.0},
        'RHOHV_F': {
            'units': 'ratio',
            'standard_name': 'cross_correlation_ratio_hv',
            'valid_max': 1.0,
            'valid_min': 0.0},
        'NCP_F': {
            'units': 'ratio',
            'standard_name': 'signal_quality',
            'valid_max': 1.0,
            'valid_min': 0.0},
        'WIDTH_F': {
            'units': 'm/s',
            'standard_name': 'spectrum_width',
            'valid_max': 45.0,
            'valid_min': 0.0},
        'PHIDP_F': {
            'units': 'degrees',
            'standard_name': 'differential_phase_hv',
            'valid_max': 80.0,
            'valid_min': -45.0},
        'VEL_COR': {
            'units': 'm/s',
            'standard_name': (
                'radial_velocity_of_scatterers_away_from_instrument'),
            'valid_max': 45.0,
            'valid_min': -45.0}
    }
    for variable_name in cgfile.variables.keys():
        if debug:
            print "doing ", variable_name
        if is_moment(variable_name, moment_fixes):
            print "doing ", variable_name
            for attr in moment_fixes[variable_name].keys():
                setattr(cgfile.variables[variable_name], attr,
                        moment_fixes[variable_name][attr])


def csapr_standard_names():
    prop_names = {
        'DBZ_F': 'reflectivity_horizontal',
        'VEL_F': 'mean_doppler_velocity',
        'WIDTH_F': 'doppler_spectral_width',
        'ZDR_F': 'diff_reflectivity',
        'RHOHV_F': 'copol_coeff',
        'NCP_F': 'norm_coherent_power',
        'KDP_F': 'diff_phase',
        'PHIDP_F': 'dp_phase_shift',
        'VEL_COR': 'corrected_mean_doppler_velocity',
        'PHIDP_UNF': 'unfolded_dp_phase_shift',
        'KDP_SOB': 'recalculated_diff_phase',
        'DBZ_AC': 'attenuation_corrected_reflectivity_horizontal'}
    return prop_names


def fix_variables_csapr(cgfile):
    debug = True
    print "go"
    csapr_names = csapr_standard_names()
    moment_fixes = {
        'DBZ_F': {
            'units': 'dBZ',
            'standard_name': 'equivalent_reflectivity_factor',
            'long_name': 'equivalent_reflectivity_factor',
            'valid_max': 80.0,
            'valid_min': -45.0},
        'VEL_F': {
            'units': 'm/s',
            'standard_name': (
                'radial_velocity_of_scatterers_away_from_instrument'),
            'long_name': (
                'radial_velocity_of_scatterers_away_from_instrument'),
            'valid_max': 95.0,
            'valid_min': -95.0},
        'KDP_F': {
            'units': 'degrees/km',
            'standard_name': 'specific_differential_phase_hv',
            'long_name': 'specific_differential_phase_hv',
            'valid_max': 20.0,
            'valid_min': -10.0},
        'ZDR_F': {
            'units': 'dB',
            'standard_name': 'log_differential_reflectivity_hv',
            'long_name': 'log_differential_reflectivity_hv',
            'valid_max': 8.0,
            'valid_min': -6.0},
        'RHOHV_F': {
            'units': 'ratio',
            'standard_name': 'cross_correlation_ratio_hv',
            'long_name': 'cross_correlation_ratio_hv',
            'valid_max': 1.0,
            'valid_min': 0.0},
        'NCP_F': {
            'units': 'ratio',
            'standard_name': 'signal_quality',
            'long_name': 'signal_quality',
            'valid_max': 1.0,
            'valid_min': 0.0,
            'comment': 'Also know as Normalized Coherent Power'},
        'WIDTH_F': {
            'units': 'm/s',
            'standard_name': 'spectrum_width',
            'long_name': 'spectrum_width',
            'valid_max': 45.0,
            'valid_min': 0.0},
        'PHIDP_F': {
            'units': 'degrees',
            'standard_name': 'differential_phase_hv',
            'long_name': 'differential_phase_hv',
            'valid_max': 180.0,
            'valid_min': -180.0},
        'VEL_COR': {
            'units': 'm/s',
            'standard_name': (
                'radial_velocity_of_scatterers_away_from_instrument'),
            'long_name': (
                'radial_velocity_of_scatterers_away_from_instrument'),
            'valid_max': 45.0,
            'valid_min': -45.0},
        'PHIDP_UNF': {
            'units': 'degrees',
            'standard_name': 'differential_phase_hv',
            'long_name': 'differential_phase_hv',
            'valid_max': 480.0,
            'valid_min': 0.0},
        'DBZ_AC': {
            'units': 'dBZ',
            'standard_name': 'equivalent_reflectivity_factor',
            'long_name': 'equivalent_reflectivity_factor',
            'valid_max': 80.0,
            'valid_min': -45.0},
        'KDP_SOB': {
            'units': 'degrees/km',
            'standard_name': 'specific_differential_phase_hv',
            'long_name': 'specific_differential_phase_hv',
            'valid_max': 20.0,
            'valid_min': -1.0}
    }
    for long_variable_name in list(set(cgfile.variables.keys()) & set(csapr_names.values())):
        variable_name = csapr_names.keys()[csapr_names.values().index(
            long_variable_name)]
        if debug:
            print "doing ", variable_name
        if is_moment(variable_name, moment_fixes):
            print "doing ", variable_name
            for attr in moment_fixes[variable_name].keys():
                setattr(cgfile.variables[long_variable_name], attr,
                        moment_fixes[variable_name][attr])


def cf_radial_coords(ncobj, radarobj):
    print lats.shape
    latvar[:] = lats
    print lons.shape
    lonvar[:] = lons
    print xar.shape
    xvar[:] = xar
    print yar.shape
    yvar[:] = yar
    print levs.shape
    levvar[:] = levs
    xvar.axis = "X"
    xvar.long_name = "x-coordinate in Cartesian system"
    xvar.standard_name = xvar.long_name
    xvar.units = "m"
    yvar.axis = "Y"
    yvar.long_name = "y-coordinate in Cartesian system"
    yvar.standard_name = yvar.long_name
    yvar.units = "m"
    levvar.long_name = "height"
    levvar.standard_name = levvar.long_name
    levvar.units = "meter"
    lonvar.long_name = "longitude"
    lonvar.standard_name = lonvar.long_name
    lonvar.units = "degrees_east"
    latvar.long_name = "latitude"
    latvar.standard_name = latvar.long_name
    latvar.units = "degrees_north"
    levvar.positive = "up"


def trans_dict_as_ncattr(mydict, ncobj, trans):
    print trans
    for attr in trans:
        ncobj.setncattr(attr, mydict[attr])


def write_radar(ncobj, radarobj, **kwargs):
    """
    Writes a Py-ART antenna coordinate radar object to a CF-Radial
    complaint netcdf file)
    """
    __DoD_version__ = "0.5"

    #Set up some nice formatting tools
    runtime = dict([(key, getattr(dt.datetime.now(), key)) for key in
                   ['year', 'month', 'day', 'hour', 'minute', 'second']])
    runtime.update({'strmon': dt.datetime.now().strftime('%b')})
    runtime.update({'user': getpass.getuser(),
                    'machine': socket.gethostname(), 'exec': sys.argv[0]})

    debug = kwargs.get('debug', False)

    #Set up some nice formatting tools
    runtime = dict([(key, getattr(dt.datetime.now(), key)) for key in
                   ['year', 'month', 'day', 'hour', 'minute', 'second']])
    runtime.update({'strmon': dt.datetime.now().strftime('%b')})
    runtime.update({'user': getpass.getuser(),
                    'machine': socket.gethostname(), 'exec': sys.argv[0]})

    #create dimensions with time first
    ncobj.createDimension('time', len(radarobj.time['data']))
    ncobj.createDimension('range', radarobj.ngates)
    ncobj.createDimension('sweep', radarobj.nsweeps)
    ncobj.createDimension('string_length_24', 24)

    #create axis variables with time first
    time_var = ncobj.createVariable('time', 'double', ('time',))
    time_var[:] = radarobj.time['data']
    trans_dict_as_ncattr(radarobj.time, time_var,
                         list(set(['units', 'comment', 'calendar',
                                   'standard_name', 'long_name']) &
                              set(radarobj.time.keys())))
    range_var = ncobj.createVariable('range', 'float', ('range',))
    range_var[:] = radarobj.range['data']
    trans_dict_as_ncattr(radarobj.range, range_var,
                         list(set(['units', 'comment',  'standard_name',
                                   'long_name',
                                   'meters_to_center_of_first_gate',
                                   'meters_between_gates']) &
                              set(radarobj.range.keys())))
    azimuth_var = ncobj.createVariable('azimuth', 'float', ('time',))
    azimuth_var[:] = radarobj.azimuth['data']
    trans_dict_as_ncattr(radarobj.azimuth, azimuth_var,
                         list(set(['units', 'comment',  'standard_name',
                                   'long_name']) &
                              set(radarobj.azimuth.keys())))
    elevation_var = ncobj.createVariable('elevation', 'float', ('time',))
    elevation_var[:] = radarobj.elevation['data']
    trans_dict_as_ncattr(radarobj.elevation, elevation_var,
                         list(set(['units', 'comment',  'standard_name',
                                   'long_name']) &
                              set(radarobj.elevation.keys())))

    #create moment variables
    moments = radarobj.fields.keys()
    my_nc_vars = {}
    for moment in moments:
        my_nc_vars.update({moment: ncobj.createVariable(
            moment, 'float32', ('time', 'range'), fill_value=-9999)})

    #populate metadata for moment variables
    for moment in moments:
        print moment
        want = list(set(radarobj.fields[moment].keys()) &
                    set(['long_name', 'standard_name', 'units', 'comment',
                         'add_offset', 'scale_factor', 'valid_min',
                         'valid_max', 'axis', 'meta_group', 'coordinates']))
        print want
        trans_dict_as_ncattr(radarobj.fields[moment], my_nc_vars[moment], want)

    #populate sweep parameters
    sweep_params = radarobj.sweep_info.keys()
    sweep_types = {
        'fixed_angle': 'float',
        'sweep_start_ray_index': 'i4',
        'sweep_end_ray_index': 'i4',
        'sweep_mode': 'S1',
        'sweep_number': 'i4'}
    dims = {
        'fixed_angle': 'sweep',
        'sweep_start_ray_index': 'sweep',
        'sweep_end_ray_index': 'sweep',
        'sweep_mode': ('sweep', 'string_length_24'),
        'sweep_number': 'sweep'}
    for item in sweep_params:
        print item, sweep_types[item]
        this_var = ncobj.createVariable(item, sweep_types[item], dims[item])
        print "made"
        trans_dict_as_ncattr(radarobj.sweep_info[item], this_var,
                             list(set(['units', 'comment',  'standard_name',
                                       'long_name']) &
                                  set(radarobj.sweep_info[item].keys())))
        print "Transfered"
        if item != "sweep_mode":
            this_var[:] = radarobj.sweep_info[item]['data']
        else:
            print netCDF4.stringtochar(np.array(
                radarobj.sweep_info[item]['data'])).shape
            this_var[:] = netCDF4.stringtochar(np.array(
                radarobj.sweep_info[item]['data']))
    #populate instr params
    #popuate location params
    if 'sort' in dir(radarobj.location['latitude']['data']):
        # we have an array, aka a moving platform
        for var in radarobj.location.keys():
            platform_is_mobile = 'True'
            this_var = ncobj.createVariable(var, 'double', ('time',))
            trans_dict_as_ncattr(radarobj.location[var], this_var,
                                 list(set(['units', 'comment',
                                           'standard_name', 'long_name']) &
                                      set(radarobj.location[var].keys())))
            this_var[:] = radarobj.location[var]['data']
    else:
        # we have a single float, ie a fixed antenna
        for var in radarobj.location.keys():
            platform_is_mobile = 'False'
            this_var = ncobj.createVariable(var, 'double')
            trans_dict_as_ncattr(radarobj.location[var], this_var,
                                 list(set(['units', 'comment',
                                           'standard_name', 'long_name']) &
                                      set(radarobj.location[var].keys())))
            this_var[:] = radarobj.location[var]['data']

    #append global header data
    trans_dict_as_ncattr(radarobj.metadata, ncobj, radarobj.metadata.keys())
    ncobj.platform_is_mobile = platform_is_mobile
    ncobj.history = "created by user %(user)s on %(machine)s at %(day)d-%(strmon)s-%(year)d,%(hour)d:%(minute)02d:%(second)02d using %(exec)s" % runtime
    ncobj.conventions = "CF/Radial"

    #populate data for moment variables
    for moment in moments:
        my_nc_vars[moment][:] = radarobj.fields[moment]['data']


def write_radar4(ncobj, radarobj, **kwargs):
    """
    Writes a Py-ART antenna coordinate radar object to a CF-Radial complaint
    netcdf file)
    """
    __DoD_version__ = "0.5"

    # Set up some nice formatting tools
    runtime = dict([(key, getattr(dt.datetime.now(), key)) for key in
                   ['year', 'month', 'day', 'hour', 'minute', 'second']])
    runtime.update({'strmon': dt.datetime.now().strftime('%b')})
    runtime.update({'user': getpass.getuser(),
                    'machine': socket.gethostname(), 'exec': sys.argv[0]})
    debug = kwargs.get('debug', False)

    # Set up some nice formatting tools
    runtime = dict([(key, getattr(dt.datetime.now(), key)) for key in
                   ['year', 'month', 'day', 'hour', 'minute', 'second']])
    runtime.update({'strmon': dt.datetime.now().strftime('%b')})
    runtime.update({'user': getpass.getuser(),
                    'machine': socket.gethostname(), 'exec': sys.argv[0]})

    #create dimensions with time first
    ncobj.createDimension('time', len(radarobj.time['data']))
    ncobj.createDimension('range', radarobj.ngates)
    ncobj.createDimension('sweep', radarobj.nsweeps)
    ncobj.createDimension('string_length_24', 24)
    ncobj.createDimension('string_length_32', 32)
    ncobj.createDimension('string_length_12', 12)
    if 'frequency' in radarobj.inst_params.keys():
       ncobj.createDimension('frequency', len(radarobj.inst_params['frequency']['data']))
    
    

    #create axis variables with time first
    time_var = ncobj.createVariable('time', 'double', ('time',))
    time_var[:] = radarobj.time['data']
    trans_dict_as_ncattr(radarobj.time, time_var,
                         list(set(['units', 'comment', 'calendar',
                                   'standard_name', 'long_name']) &
                              set(radarobj.time.keys())))
    range_var = ncobj.createVariable('range', 'float', ('range',))
    range_var[:] = radarobj.range['data']
    trans_dict_as_ncattr(radarobj.range, range_var,
                         list(set(['units', 'comment',  'standard_name',
                                   'long_name',
                                   'meters_to_center_of_first_gate',
                                   'meters_between_gates']) &
                              set(radarobj.range.keys())))
    azimuth_var = ncobj.createVariable('azimuth', 'float', ('time',))
    azimuth_var[:] = radarobj.azimuth['data']
    trans_dict_as_ncattr(radarobj.azimuth, azimuth_var,
                         list(set(['units', 'comment',  'standard_name',
                                   'long_name']) &
                              set(radarobj.azimuth.keys())))
    elevation_var = ncobj.createVariable('elevation', 'float', ('time',))
    elevation_var[:] = radarobj.elevation['data']
    trans_dict_as_ncattr(radarobj.elevation, elevation_var,
                         list(set(['units', 'comment',  'standard_name',
                                   'long_name']) &
                              set(radarobj.elevation.keys())))
    # create moment variables
    moments = radarobj.fields.keys()
    my_nc_vars = {}
    for moment in moments:
        if 'least_significant_digit' in radarobj.fields.keys():
            my_nc_vars.update({moment: ncobj.createVariable(
                moment, 'float32', ('time', 'range'), fill_value=-9999.0,
                zlib=True, least_significant_digit=radarobj.fields[
                    moment]['least_significant_digit'])})
        else:
            my_nc_vars.update({moment: ncobj.createVariable(
                moment, 'float32', ('time', 'range'), fill_value=-9999.0,
                zlib=True)})
	  

    #populate metadata for moment variables
    for moment in moments:
        print moment
        want = list(set(radarobj.fields[moment].keys()) &
                    set(['long_name', 'standard_name', 'units', 'comment',
                         'add_offset', 'scale_factor', 'valid_min',
                         'valid_max', 'axis', 'meta_group', 'coordinates']))
        print want
        trans_dict_as_ncattr(radarobj.fields[moment], my_nc_vars[moment], want)
    #populate data for moment variables
    for moment in moments:
        print my_nc_vars[moment]._FillValue
        print moment
        if '_fill_value' in dir(radarobj.fields[moment]['data']):
            radarobj.fields[moment]['data']._fill_value = None
            # gets around a bug
        my_nc_vars[moment][:] = radarobj.fields[moment]['data']

    #populate sweep parameters
    sweep_params = radarobj.sweep_info.keys()
    sweep_types = {
        'fixed_angle': 'float',
        'sweep_start_ray_index': 'i4',
        'sweep_end_ray_index': 'i4',
        'sweep_mode': 'S1',
        'sweep_number': 'i4'}
    char_shape= netCDF4.stringtochar(np.array(
                radarobj.sweep_info['sweep_mode']['data'])).shape[1]
    string_string='string_length_%(d)d' % {'d':char_shape}
    dims = {'fixed_angle': 'sweep',
            'sweep_start_ray_index': 'sweep',
            'sweep_end_ray_index': 'sweep',
            'sweep_mode': ('sweep',string_string),
            'sweep_number': 'sweep'}
    for item in sweep_params:
        print item, sweep_types[item]
        this_var = ncobj.createVariable(
            item, sweep_types[item], dims[item], zlib=True)
        print "made"
        trans_dict_as_ncattr(radarobj.sweep_info[item], this_var,
                             list(set(['units', 'comment',  'standard_name',
                                       'long_name']) &
                                  set(radarobj.sweep_info[item].keys())))
        print "Transfered"
        if item != "sweep_mode":
            this_var[:] = radarobj.sweep_info[item]['data']
        else:
            print radarobj.sweep_info[item]['data'], 'moom', item
            print netCDF4.stringtochar(np.array(
                radarobj.sweep_info[item]['data'])).shape
            print this_var[:].shape
            this_var[:] = netCDF4.stringtochar(np.array(
                radarobj.sweep_info[item]['data']))

    #populate instr params
    
    #inst_shape= netCDF4.stringtochar(np.array(
    #            radarobj.sweep_info['prt_mode']['data'])).shape[1]
    #inst_str='string_length_%(d)d' % {'d':inst_shape}

    inst_types = {
        'frequency': 'float',
        'follow_mode': 'S1',
        'pulse_width': 'float',
        'prt_mode': 'S1',
        'prt': 'float',
        'prt_ratio': 'float',
        'polarization_mode': 'S1',
        'nyquist_velocity': 'float',
        'unambiguous_range': 'float',
        'n_samples': 'i4'}
    inst_dims = {
        'frequency': 'frequency',
        'follow_mode': ('sweep', string_string),
        'pulse_width': 'time',
        'prt_mode': ('sweep', string_string),
        'prt': 'time',
        'prt_ratio': 'time',
        'polarization_mode': ('sweep', string_string),
        'nyquist_velocity': 'time',
        'unambiguous_range': 'time',
        'n_samples': 'time'}
    for key in radarobj.inst_params.keys():
        this_var = ncobj.createVariable(
            key, inst_types[key], inst_dims[key], zlib=True)
        trans_dict_as_ncattr(radarobj.inst_params[key], this_var,
                             list(set(['units', 'comment',  'standard_name',
                                       'long_name']) &
                                  set(radarobj.inst_params[key].keys())))
        this_var.meta_group = "instrument_parameters"
        if key not in ["follow_mode", "prt_mode", "polarization_mode"]:
	    print key
            this_var[:] = radarobj.inst_params[key]['data']
        else:
            print radarobj.inst_params[key]['data'], 'moom'
            print netCDF4.stringtochar(np.array(
                radarobj.inst_params[key]['data'])).shape
            this_var[:] = netCDF4.stringtochar(np.array(
                radarobj.inst_params[key]['data']))

    #popuate location params
    if 'sort' in dir(radarobj.location['latitude']['data']):
        if len(radarobj.location['latitude']['data']) > 1:
            # we have an array, aka a moving platform
            for var in radarobj.location.keys():
                platform_is_mobile = 'True'
                this_var = ncobj.createVariable(var, 'double', ('time',),
                                                zlib=True)
                trans_dict_as_ncattr(radarobj.location[var], this_var,
                                     list(set(['units', 'comment',
                                              'standard_name', 'long_name']) &
                                          set(radarobj.location[var].keys())))
                this_var[:] = radarobj.location[var]['data']
        else:
            # we have a single float, ie a fixed antenna
            for var in radarobj.location.keys():
                platform_is_mobile = 'False'
                this_var = ncobj.createVariable(var, 'double', zlib=True)
                trans_dict_as_ncattr(radarobj.location[var], this_var,
                                     list(set(['units', 'comment',
                                               'standard_name', 'long_name']) &
                                          set(radarobj.location[var].keys())))
                this_var[:] = radarobj.location[var]['data'][:]
    else:
        # we have a single float, ie a fixed antenna
        for var in radarobj.location.keys():
            platform_is_mobile = 'False'
            this_var = ncobj.createVariable(var, 'double', zlib=True)
            trans_dict_as_ncattr(radarobj.location[var], this_var,
                                 list(set(['units', 'comment',
                                           'standard_name', 'long_name']) &
                                      set(radarobj.location[var].keys())))
            this_var[:] = radarobj.location[var]['data']
    #append global header data
    trans_dict_as_ncattr(radarobj.metadata, ncobj, radarobj.metadata.keys())
    ncobj.platform_is_mobile = platform_is_mobile
    ncobj.history = "created by user %(user)s on %(machine)s at %(day)d-%(strmon)s-%(year)d,%(hour)d:%(minute)02d:%(second)02d using %(exec)s" % runtime
    ncobj.conventions = "CF/Radial"


def is_moment(varname, moment_fixes):
    moments = moment_fixes.keys()
    return True in [foo in varname for foo in moments]


def is_radar(varname, radar_list):
    return True in [foo in varname for foo in radar_list]


def append_global_metatdata_csapr(cgfile):
    cgfile.Conventions = 'CF 1.5'
    cgfile.title = 'Radar moments mapped to a Cartesian grid'
    cgfile.history = 'Dealiasing done with the U Washington 4DD code, attenuation correction using the ZPHI method in high (45+dBz) reflectivities, phidp only elsewhere'
    cgfile.source = 'Mapped moments from the C Band Scanning ARM Radar. '
    cgfile.institution = 'Atmospheric Radiation Measurement Climate Facility, United States Department of Energy'
    cgfile.comment = """Evaluation mapped moments from the C-SAPR radar. Please note that this data is under
    active development and there may be errors or issues in the data. Please report any issues, errors or suggestions
    to the precipitation radar data translator scollis@anl.gov"""
    cgfile.site = 'SGP'


def append_global_metatdata_xsapr(cgfile):
    cgfile.Conventions = 'CF 1.5'
    cgfile.title = 'Radar moments mapped to a Cartesian grid'
    cgfile.history = 'Dealiasing done with the U Washington 4DD code'
    cgfile.source = 'Mapped moments from the X Band Scanning ARM Radar. '
    cgfile.institution = 'Atmospheric Radiation Measurement Climate Facility, United States Department of Energy'
    cgfile.comment = """Evaluation mapped moments from the X-SAPR radar. Please note that this data is under
    active development and there may be errors or issues in the data. Please report any issues, errors or suggestions
    to the precipitation radar data translator scollis@anl.gov"""
    cgfile.site = 'SGP'


def append_coords(cgfile, xar, yar, zar, radar_loc):
    Re = 6371.0 * 1000.0
    rad_at_radar = Re * np.sin(np.pi / 2.0 -
                               np.abs(radar_loc[0] * np.pi / 180.0))
    # ax_radius(float(lat_cpol), units='degrees')
    lons = radar_loc[1] + 360.0 * xar / (rad_at_radar * 2.0 * np.pi)
    lats = radar_loc[0] + 360.0 * yar / (Re * 2.0 * np.pi)
    levs = zar
    levvar = cgfile.createVariable('z_disp', 'float32', ('nz',))
    yvar = cgfile.createVariable('y_disp', 'float32', ('ny',))
    xvar = cgfile.createVariable('x_disp', 'float32', ('nx',))
    latvar = cgfile.createVariable('lat', 'float32', ('nx',))
    lonvar = cgfile.createVariable('lon', 'float32', ('ny',))
    print lats.shape
    latvar[:] = lats
    print lons.shape
    lonvar[:] = lons
    print xar.shape
    xvar[:] = xar
    print yar.shape
    yvar[:] = yar
    print levs.shape
    levvar[:] = levs
    xvar.axis = "X"
    xvar.long_name = "x-coordinate in Cartesian system"
    xvar.standard_name = xvar.long_name
    xvar.units = "m"
    xvar.comment = "X Displacement from the central facility"
    yvar.axis = "Y"
    yvar.comment = "Y Displacement from the central facility"
    yvar.long_name = "y-coordinate in Cartesian system"
    yvar.standard_name = yvar.long_name
    yvar.units = "m"
    levvar.long_name = "height"
    levvar.standard_name = levvar.long_name
    levvar.units = "meter"
    lonvar.long_name = "longitude"
    lonvar.standard_name = lonvar.long_name
    lonvar.units = "degrees_east"
    latvar.long_name = "latitude"
    latvar.standard_name = latvar.long_name
    latvar.units = "degrees_north"
    levvar.positive = "up"


def fix_time(cgfile, time):
    times = cgfile.createVariable('time', 'double', ('time', ))
    times.units = 'seconds since 1970-01-01 00:00:00.0'
    times.calendar = 'gregorian'
    times.standard_name = 'time'
    print "time here is"
    print time
    print netCDF4.date2num(time, units=times.units, calendar=times.calendar)
    times[0] = netCDF4.date2num(time, units=times.units,
                                calendar=times.calendar)
    print times


def dt_to_dict(dt, **kwargs):
    pref = kwargs.get('pref', '')
    return dict([(pref+key, getattr(dt, key)) for key in
                ['year', 'month', 'day', 'hour', 'minute', 'second']])


def dms_to_d(dms):
    return dms[0] + (dms[1] + dms[2] / 60.0) / 60.0


def corner_to_point(corner, point):
    Re = 6371.0 * 1000.0
    Rc = ax_radius(point[0], units='degrees')
    #print Rc/Re
    y = ((point[0] - corner[0]) / 360.0) * np.pi * 2.0 * Re
    x = ((point[1] - corner[1]) / 360.0) * np.pi * 2.0 * Rc
    return x, y


def ax_radius(lat, units='radians'):
    #Determine the radius of a circle of constant longitude at a certain
    #Latitude
    Re = 6371.0 * 1000.0
    if units == 'degrees':
        const = np.pi / 180.0
    else:
        const = 1.0
    R = Re * np.sin(np.pi / 2.0 - np.abs(lat * const))
    return R


def rsl_to_arm_netcdf(rslobj, ofilename, **kwargs):
    """
    Save an enhanced MDV object (containing gridded data) to a
    ARM complaint netcdf file
    DoD Version 0.5
    """
    __DoD_version__ = "0.5"
    debug = kwargs.get('debug', False)

    #Set up some nice formatting tools
    runtime = dict([(key, getattr(dt.datetime.now(), key)) for key in
                    ['year', 'month', 'day', 'hour', 'minute', 'second']])
    runtime.update({'strmon': dt.datetime.now().strftime('%b')})
    runtime.update({'user': getpass.getuser(),
                    'machine': socket.gethostname(), 'exec': sys.argv[0]})
    datastream_format = "%(name)s : %(version)s : %(s_year)04d%(s_month)02d%(s_day)02d.%(s_hour)02d%(s_minute)02d%(s_second)02d-%(e_year)04d%(e_month)02d%(e_day)02d.%(e_hour)02d%(e_minute)02d%(e_second)02d;\n"
    dd2 = dt_to_dict(netCDF4.num2date(
        rslobj.sounding_used.variables['time'][:],
        units=rslobj.sounding_used.variables['time'].units,
        calendar='gregorian')[0], pref='s_')

    print "eat"
    dd2.update(dt_to_dict(netCDF4.num2date(
        rslobj.sounding_used.variables['time'][:],
        units=rslobj.sounding_used.variables['time'].units,
        calendar='gregorian')[-1], pref='e_'))
    dd2.update({'version': rslobj.sounding_used.process_version,
                'name': rslobj.sounding_used.site_id +
                rslobj.sounding_used.dod_version.split('-')[0]})
    ds2 = datastream_format % dd2
    print "cake"
    dd1 = dt_to_dict(rslobj.contents.datetime, pref='s_')
    dd1.update(dt_to_dict(rslobj.contents.datetime, pref='e_'))
    dd1.update({'version': 'raw_stream',
                'name': rslobj.filename.split('/')[-1].lower()[0:3] +
                rslobj.contents.h.radar_name})
    ds1 = datastream_format % dd1
    print ds1
    #open the file, netcdf3 for now
    print "fopo"
    if debug:
        print "opening ", ofilename
    ofile = netCDF4.Dataset(ofilename, 'w', format='NETCDF3_CLASSIC')
    # create the dimensions, time must be first and unlimited,
    # no one letter dimensions
    nz, ny, nx = rslobj.grids[rslobj.grids.keys()[0]].shape
    ofile.createDimension('time', None)
    ofile.createDimension('nz', nz)
    ofile.createDimension('ny', ny)
    ofile.createDimension('nx', nx)
    #create the indexing variables, time must be first,
    fix_time(ofile, rslobj.contents.datetime)
    xar = np.linspace(rslobj.grid_coords['xr'][0],
                      rslobj.grid_coords['xr'][1], nx)
    yar = np.linspace(rslobj.grid_coords['yr'][0],
                      rslobj.grid_coords['yr'][1], ny)
    zar = np.linspace(rslobj.grid_coords['zr'][0],
                      rslobj.grid_coords['zr'][1], nz)
    append_coords(ofile, xar, yar, zar, rslobj.origin)
    #[dms_to_d((rslobj.contents.h.latd, rslobj.contents.h.latm,
    #           rslobj.contents.h.lats)),
    # dms_to_d((rslobj.contents.h.lond, rslobj.contents.h.lonm,
    #           rslobj.contents.h.lons))]
    #write variables
    csapr_names = csapr_standard_names()
    available_pars = list(set(csapr_names.keys()) & set(rslobj.grids.keys()))
    vvars = [ofile.createVariable(csapr_names[par], np.float,
             ('time', 'nz', 'ny', 'nx'), fill_value=-9999) for par in
             available_pars]
    fix_variables_csapr(ofile)

    #write global headers, history must be last
    ofile.process_version = __version__
    ofile.conventions = "CF 1.5"
    ofile.command_line = rslobj.command_line
    ofile.dod_version = __DoD_version__
    sgp_map = {
        'NW': 'I6: Deer Creek, Oklahoma',
        'SW': 'I5: Garber, Oklahoma',
        'SE': 'I4: Billings, Oklahoma'}
    if rslobj.filename.split('/')[-1][1:3] in sgp_map.keys():
        ofile.site_id = sgp_map[rslobj.filename.split('/')[-1][1:3]]
        ofile.facility_id = 'sgp'
    else:
        ofile.site_id = "Alaska, to be fixed"
        ofile.facility_id = 'nsa'
    ofile.pyart_procs = rslobj.pyart_procs
    ofile.input_datastreams_num = "2"
    ofile.input_datastreams = ds1 + ds2
    ofile.history = "created by user %(user)s on %(machine)s at %(day)d-%(strmon)s-%(year)d,%(hour)d:%(minute)02d:%(second)02d using %(exec)s" % runtime
    #write the actual data
    for i in range(len(available_pars)):
        print "writing: ", available_pars[i]
        vvars[i][0, :, :, :] = np.ma.masked_array(
            rslobj.grids[available_pars[i]],
            np.isnan(rslobj.grids[available_pars[i]]))
    #write out and close
    ofile.close()


def mdv_to_arm_netcdf(mdvobj, ofilename, **kwargs):
    """Save an enhanced MDV object (containing gridded data) to a
    ARM complaint netcdf file
    DoD Version 0.5
    """
    __DoD_version__ = "0.5"
    debug = kwargs.get('debug', False)

    #Set up some nice formatting tools
    runtime = dict([(key, getattr(dt.datetime.now(), key)) for key in
                   ['year', 'month', 'day', 'hour', 'minute', 'second']])
    runtime.update({'strmon': dt.datetime.now().strftime('%b')})
    runtime.update({'user': getpass.getuser(),
                    'machine': socket.gethostname(), 'exec': sys.argv[0]})
    datastream_format = "%(name)s : %(version)s : %(s_year)04d%(s_month)02d%(s_day)02d.%(s_hour)02d%(s_minute)02d%(s_second)02d-%(e_year)04d%(e_month)02d%(e_day)02d.%(e_hour)02d%(e_minute)02d%(e_second)02d;\n"
    dd2 = dt_to_dict(netCDF4.num2date(
        mdvobj.sounding_used.variables['time'][:],
        units=mdvobj.sounding_used.variables['time'].units,
        calendar='gregorian')[0], pref='s_')
    dd2.update(dt_to_dict(netCDF4.num2date(
        mdvobj.sounding_used.variables['time'][:],
        units=mdvobj.sounding_used.variables['time'].units,
        calendar='gregorian')[-1], pref='e_'))
    dd2.update({'version': mdvobj.sounding_used.process_version,
                'name': mdvobj.sounding_used.site_id +
                mdvobj.sounding_used.dod_version.split('-')[0]})
    ds2 = datastream_format % dd2
    dd1 = dt_to_dict(mdvobj.times['time_begin'], pref='s_')
    dd1.update(dt_to_dict(mdvobj.times['time_end'], pref='e_'))
    dd1.update({'version': str(mdvobj.master_header['revision_number']),
                'name': mdvobj.master_header['data_set_source'].split(' ')[1] +
                mdvobj.master_header['data_set_name']})
    ds1 = datastream_format % dd1
    print ds1
    #open the file, netcdf3 for now
    if debug:
        print "opening ", ofilename
    ofile = netCDF4.Dataset(ofilename, 'w', format='NETCDF3_CLASSIC')
    # create the dimensions, time must be first and unlimited,
    # no one letter dimensions
    nz, ny, nx = mdvobj.grids[mdvobj.grids.keys()[0]].shape
    ofile.createDimension('time', None)
    ofile.createDimension('nz', nz)
    ofile.createDimension('ny', ny)
    ofile.createDimension('nx', nx)

    #create the indexing variables, time must be first,
    print "Fixing time"
    print mdvobj.times['time_begin']
    fix_time(ofile, mdvobj.times['time_begin'])
    print "time fixed"
    xar = np.linspace(mdvobj.grid_coords['xr'][0],
                      mdvobj.grid_coords['xr'][1], nx)
    yar = np.linspace(mdvobj.grid_coords['yr'][0],
                      mdvobj.grid_coords['yr'][1], ny)
    zar = np.linspace(mdvobj.grid_coords['zr'][0],
                      mdvobj.grid_coords['zr'][1], nz)
    append_coords(ofile, xar, yar, zar, mdvobj.grid_origin)

    #write variables
    csapr_names = csapr_standard_names()
    available_pars = list(set(csapr_names.keys()) & set(mdvobj.grids.keys()))
    vvars = [ofile.createVariable(csapr_names[par], np.float,
             ('time', 'nz', 'ny', 'nx'), fill_value=-9999) for par in
             available_pars]
    fix_variables_csapr(ofile)

    #write global headers, history must be last
    ofile.process_version = __version__
    ofile.conventions = "CF 1.5"
    ofile.command_line = mdvobj.command_line
    ofile.dod_version = __DoD_version__
    if mdvobj.master_header['data_set_source'].split(' ')[1].lower() == 'sgp':
        ofile.site_id = "I7: Nardin, Oklahoma"
    else:
        ofile.site_id = "I1: Lombrum, Manus Island"
    ofile.facility_id = mdvobj.master_header[
        'data_set_source'].split(' ')[1].lower()
    ofile.pyart_procs = mdvobj.pyart_procs
    ofile.input_datastreams_num = "2"
    ofile.input_datastreams = ds1 + ds2
    ofile.history = "created by user %(user)s on %(machine)s at %(day)d-%(strmon)s-%(year)d,%(hour)d:%(minute)02d:%(second)02d using %(exec)s" % runtime

    #write the actual data
    for i in range(len(available_pars)):
        vvars[i][0, :, :, :] = np.ma.masked_array(
            mdvobj.grids[available_pars[i]],
            np.isnan(mdvobj.grids[available_pars[i]]))
    #write out and close
    ofile.close()


def save_netcdf_cube(grids, output_file, myfile, xr, yr, zr, origin, **kwargs):
    meta = kwargs.get('meta', 'csapr_ac')
    print "opening ", output_file
    ofile = netCDF4.Dataset(output_file, 'w', format='NETCDF3_CLASSIC')
    nz, ny, nx = grids[grids.keys()[0]].shape
    ofile.createDimension('x', nx)
    ofile.createDimension('y', ny)
    ofile.createDimension('z', nz)
    ofile.createDimension('time', 1)
    csapr_names = csapr_standard_names()
    available_pars = list(set(csapr_names.keys()) & set(grids.keys()))
    print available_pars
    vvars = [ofile.createVariable(csapr_names[par], np.float,
             ('time', 'z', 'y', 'x'), fill_value=-999) for par in
             available_pars]
    for i in range(len(available_pars)):
        vvars[i][0, :, :, :] = np.ma.masked_array(
            grids[available_pars[i]], np.isnan(grids[available_pars[i]]))
    fix_variables_csapr(ofile)
    if meta == 'csapr_ac':
        append_global_metatdata_csapr(ofile)
    elif meta == 'xsapr':
        append_global_metatdata_xsapr(ofile)
    xar = np.linspace(xr[0], xr[1], nx)
    yar = np.linspace(yr[0], yr[1], ny)
    zar = np.linspace(zr[0], zr[1], nz)
    if "radar_info" in dir(myfile):  # this is a mdv object
        radar_loc = [myfile.radar_info['latitude_deg'],
                     myfile.radar_info['longitude_deg']]
        fix_time(ofile, myfile.times['time_begin'])
    else:
        # this is an RSL object
        radar_loc = [dms_to_d((myfile.contents.h.latd,
                     myfile.contents.h.latm, myfile.contents.h.lats)),
                     dms_to_d((myfile.contents.h.lond, myfile.contents.h.lonm,
                               myfile.contents.h.lons))]
        fix_time(ofile, myfile.contents.datetime)
    append_coords(ofile, xar, yar, zar, origin)
    ofile.close()


def noncf_append_global_metatdata(cgfile):
    cgfile.Conventions = 'None, this is not a cfradial file'
    cgfile.title = 'Radar moments converted to a NetCDF cube (nsweeps, nrays, ngates)'
    cgfile.history = 'Saved by save_mdv_ncf'
    cgfile.source = 'Moments from the C Band Scanning ARM Radar. '
    cgfile.institution = 'Atmospheric Radiation Measurement Climate Facility, United States Department of Energy'
    cgfile.comment = """Moments from the C-SAPR radar. Please note that this data is under
    active development and there may be errors or issues in the data. Please report any issues, errors or suggestions
    to the precipitation radar data translator scollis@anl.gov"""
    cgfile.site = 'SGP'


def copy_meta(myfile, ofile):
    for key in myfile.radar_info.keys():
        setattr(ofile, key, myfile.radar_info[key])


def copy_axis(myfile, cgfile):
    xvar = cgfile.createVariable('x', 'float32',
                                 ('time', 'nsweeps', 'nrays', 'ngates'))
    yvar = cgfile.createVariable('y', 'float32',
                                 ('time', 'nsweeps', 'nrays', 'ngates'))
    zvar = cgfile.createVariable('z', 'float32',
                                 ('time', 'nsweeps', 'nrays', 'ngates'))
    xvar[0, :, :, :] = myfile.carts['x']
    yvar[0, :, :, :] = myfile.carts['y']
    zvar[0, :, :, :] = myfile.carts['z']
    azvar = cgfile.createVariable('azimtuth', 'float32', ('nrays',))
    rnvar = cgfile.createVariable('range', 'float32', ('ngates',))
    elvar = cgfile.createVariable('elevation', 'float32', ('nsweeps',))
    azvar[:] = myfile.az_deg
    rnvar[:] = myfile.range_km * 1000.0
    elvar[:] = myfile.el_deg
    xvar.axis = "X"
    xvar.long_name = "x-coordinate in Cartesian system"
    xvar.units = "m"
    yvar.axis = "Y"
    yvar.long_name = "y-coordinate in Cartesian system"
    yvar.units = "m"
    zvar.long_name = "height"
    zvar.units = "meter"
    zvar.positive = "up"
    azvar.long_name = "Azimuth of antenna"
    azvar.units = "degrees"
    elvar.long_name = "Elevation of antenna"
    elvar.units = "degrees"
    rnvar.long_name = "Range from antenna"
    rnvar.units = "meter"


def save_pyGrid(ncfobj, pygrid):
    """Saves a pyGrid object to a CF and ARM standard netcdf file

    usage: save_pyGrid(netCDF4_object, pyGrid_object)

    """

    #first create the time dimension

    ncfobj.createDimension('time', None)

    #Axes dimensions
    #grab the dimensions from the first moment field

    nz, ny, nx = pygrid.fields[pygrid.fields.keys()[0]]['data'].shape
    ncfobj.createDimension('nz', nz)
    ncfobj.createDimension('ny', ny)
    ncfobj.createDimension('nx', nx)

    #axes
    #Generate the variables

    vvars = [ncfobj.createVariable(
        key, np.float, ('time', 'nz', 'ny', 'nx'),
        fill_value=pygrid.fields[key]['_FillValue']) for key in
        pygrid.fields.keys()]

    # loop and populate attributes

    for field in pygrid.fields.keys():
        for meta in pygrid.fields[field].keys():
            if meta != 'data' and meta != '_FillValue':
                setattr(ncfobj.variables[field], meta,
                        pygrid.fields[field][meta])
    akeys = pygrid.axes.keys()
    akeys.sort()  # makes sure time comes first

    #For ARM compliance we want alt, lat and lon to be at the end

    for mkey in ['lat', 'lon', 'alt']:
        try:
            akeys.remove(mkey)
            akeys.append(mkey)
        except ValueError:
            print(mkey, " not existing")

    dims_lookup = {'x_disp': 'nx', 'y_disp': 'ny', 'z_disp': 'nz',
                   'time_end': 'time', 'time_start': 'time',
                   'lat': 'time', 'lon': 'time', 'alt': 'time'}
    avars = [ncfobj.createVariable(key, np.float, (dims_lookup[key], ))
             for key in akeys]

    # loop and populate attributes

    for axis in akeys:
        metakeys = pygrid.axes[axis].keys()

        #again, reorder to meet ARM standards..

        for mkey in ['units', 'long_name']:
            try:
                metakeys.remove(mkey)
                metakeys.insert(0, mkey)
            except ValueError:
                print(mkey, " not existing")
        for meta in metakeys:
            if meta != 'data':
                setattr(ncfobj.variables[axis],
                        meta, pygrid.axes[axis][meta])

    # global metadata

    if 'process_version' in pygrid.metadata.keys():
        ncfobj.process_version = pygrid.metadata['process_version']
    for meta in pygrid.metadata.keys():
        if meta != 'history' or meta != 'process_version':
            setattr(ncfobj, meta, pygrid.metadata[meta])
    ncfobj.history = pygrid.metadata['history']

    #now populate data.. we leave this until last to speed up..

    for i in range(len(akeys)):
        print(akeys[i], pygrid.axes[akeys[i]]['data'].shape, avars[i].shape)
        avars[i][:] = pygrid.axes[akeys[i]]['data']
    for i in range(len(pygrid.fields.keys())):
        vvars[i][0, :, :, :] = pygrid.fields[
            pygrid.fields.keys()[i]]['data'][:, :, :]


def save_mdv_ncf(myfile, output_file, parms, **kwargs):
    debug = kwargs.get('debug', False)
    if debug:
        print "opening ", output_file
    ofile = netCDF4.Dataset(output_file, 'w', format='NETCDF3_CLASSIC')
    if debug:
        print 'Appending global metadata'
    noncf_append_global_metatdata(ofile)
    if debug:
        print 'Appending radar metadata'
    copy_meta(myfile, ofile)
    nsweeps, nrays, ngates = (myfile.read_a_field(0)).shape
    if debug:
        print "Creating dimensions"
    ofile.createDimension('nsweeps', nsweeps)
    ofile.createDimension('nrays', nrays)
    ofile.createDimension('ngates', ngates)
    ofile.createDimension('time', 1)
    if debug:
        print "Creating axes"
    copy_axis(myfile, ofile)
    if debug:
        print "Creating variables"
    csapr_names = csapr_standard_names()
    available_pars = list(set(csapr_names.keys()) & set(parms))
    vvars = [ofile.createVariable(csapr_names[par], np.float,
             ('time', 'nsweeps', 'nrays', 'ngates'), fill_value=-999) for
             par in available_pars]
    if debug:
        print 'fixing variables'
    fix_variables_csapr(ofile)
    if debug:
        print "Filling in data"
    for i in range(len(parms)):
        if debug:
            print "Doing ", parms[i]
        mdv_data = myfile.read_a_field(myfile.fields.index(parms[i]))
        vvars[i][0, :, :, :] = np.ma.masked_array(
            mdv_data, np.isnan(mdv_data))
    if debug:
        print "fixing time"
    fix_time(ofile, myfile.times['time_begin'])
    ofile.close()
