"""
A general central radial scanning (or dwelling) instrument class

USE
---
xsapr = py4dd.RSL_anyformat_to_radar(filename)
myradar=radar.Radar(xsapr)

"""

import sys
import os
from datetime import datetime
import copy

# TODO change this to import numpy
from numpy import tile, array, isnan, where, ma, linspace, arange, zeros, \
    float32, abs, empty, append, max
from netCDF4 import date2num



def dms_to_d(dms):
    """ Degrees, minutes, seconds to degrees """
    return dms[0] + (dms[1] + dms[2] / 60.0) / 60.0


def csapr_standard_names():
    prop_names = {'DBZ_F': 'reflectivity_horizontal',
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
                  'DBZ_AC': 'attenuation_corrected_reflectivity_horizontal', }
    return prop_names


def dt_to_dict(dt, **kwargs):
    pref = kwargs.get('pref', '')
    return dict([(pref+key, getattr(dt, key)) for key in
                ['year', 'month', 'day', 'hour', 'minute', 'second']])


def defaut_mdv_metadata_map():
    """
    produce the default mappings from mdv name space to cf-radial name space
    """
    mdm = {'instrument_name': 'data_set_source', 'source': 'data_set_info'}
    return mdm


def create_cube_array(volume):
    ppi = zeros([volume.h.nsweeps, volume.sweeps[0].h.nrays,
                volume.sweeps[0].rays[0].h.nbins],
                dtype=float32) + 1.31072000e+05
    for levnum in range(volume.h.nsweeps):
        for raynum in range(volume.sweeps[0].h.nrays):
            data = volume.sweeps[levnum].rays[raynum].data
            ppi[levnum, raynum, 0:len(data)] = data
    return ppi


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


def stream_to_2d(data, sweeps, sweepe, ray_len, maxgates, nrays,
                 ray_start_index):
    time_range = ma.zeros([nrays, maxgates]) - 9999.0
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


def stream_ncvar_to_field(ncvar, sweeps, sweepe, ray_len, maxgates, nrays,
                          ray_start_index):
    outdict = {'data': stream_to_2d(
        ncvar[:], sweeps, sweepe, ray_len, maxgates, nrays, ray_start_index)}
    outdict.update(dict([(key, getattr(ncvar, key)) for key in
                   ncvar.ncattrs()]))
    return outdict


class Radar:
    """
    A class for storing antenna coordinate radar data which will interact
    nicely with CF-Radial files and other pyart code
    """
    def __init__(self, nsweeps, nrays, ngates, scan_type, naz, nele, _range,
                 azimuth, elevation, tu, cal, time, fields, sweep_info,
                 sweep_mode, sweep_number, location, inst_params, metadata):
        
        self.azimuth = azimuth
        self.elevation = elevation
        self.fields = fields
        self.location = location
        self.metadata  = metadata
        self.naz = naz
        self.nele = nele
        self.ngates = ngates
        self.nsweeps = nsweeps
        self.range = _range
        self.scan_type = scan_type
        self.sweep_info = sweep_info
        self.sweep_mode = sweep_mode
        self.sweep_number = sweep_number
        self.time = time
        self.inst_params = inst_params
        
        self.cal = cal
        self.nrays = nrays
        self.tu = tu
    
    """
    def __init__(self, radarobj, **kwargs):
    
        #check format of file object
        if 'vlevel_headers' in dir(radarobj):  # this is a mdv file
            self.mdv2rad(radarobj, **kwargs)
        elif 'contents' in dir(radarobj):  # this is a ctype, rsl object?
            if 'h' in dir(radarobj.contents):  # yes a rsl object
                self.rsl2rad(radarobj, **kwargs)
        elif 'variables' in dir(radarobj):  # netCDF file
            if 'ray_start_index' in radarobj.variables.keys():
                self.streamcf2rad(radarobj, **kwargs)
            else:
                self.cf2rad(radarobj, **kwargs)
    """
    def cf2rad(self, ncobj):
        try:
            mode = "".join(ncobj.variables['sweep_mode'][0])
        except TypeError:
            mode = "".join(ncobj.variables['sweep_mode'][0].data)
        print mode, "azimuth_surveillance    "
        if "sur" in mode:
            #ppi
            self.nsweeps = len(ncobj.variables['sweep_start_ray_index'])
            self.metadata = dict([(key, getattr(ncobj, key)) for key in
                                 ncobj.ncattrs()])
            self.scan_type = "ppi"
            self.sweep_mode = array(['ppi']*self.nsweeps)
            if len(ncobj.variables['sweep_start_ray_index']) == 1:
                self.naz = ncobj.variables['sweep_end_ray_index'][0] + 1
            else:
                self.naz = (ncobj.variables['sweep_start_ray_index'][1] -
                            ncobj.variables['sweep_start_ray_index'][0])
            self.nele = ncobj.variables['sweep_start_ray_index'].shape[0]
            self.ngates = len(ncobj.dimensions['range'])
            loc_dict = {}
            for loc_data in ['latitude', 'altitude', 'longitude']:
                loc_dict.update(
                    {loc_data: ncvar_to_field(ncobj.variables[loc_data])})
            self.location = loc_dict
            sweep_dict = {}
            for sweep_data in ['sweep_start_ray_index', 'sweep_mode', 'sweep_number', 'sweep_end_ray_index', 'fixed_angle']:
                sweep_dict.update(
                    {sweep_data: ncvar_to_field(ncobj.variables[sweep_data])})
            self.sweep_info = sweep_dict
            inst_dict = {}
            for inst_data in ['frequency', 'follow_mode', 'pulse_width', 'prt_mode', 'prt', 'prt_ratio', 'polarization_mode', 'nyquist_velocity', 'unambiguous_range', 'n_samples']:
                if inst_data in ncobj.variables.keys():
                    inst_dict.update(
                        {inst_data: (
                            ncvar_to_field(ncobj.variables[inst_data]))})
            self.inst_params = inst_dict
            self.azimuth = ncvar_to_field(ncobj.variables['azimuth'])
            self.range = ncvar_to_field(ncobj.variables['range'])
            self.elevation = ncvar_to_field(ncobj.variables['elevation'])
            self.time = ncvar_to_field(ncobj.variables['time'])
            data_fields = create_field_list(
                ncobj.variables, len(ncobj.dimensions['time']), self.ngates)
            field_dict = {}
            for field in data_fields:
                print field
                my_field = ncvar_to_field(ncobj.variables[field])
                field_dict.update({field: my_field})
            self.fields = field_dict
        if "sec" in mode:
            #sec
            self.nsweeps = len(ncobj.variables['sweep_start_ray_index'])
            self.metadata = dict([(key, getattr(ncobj, key)) for key in
                                 ncobj.ncattrs()])
            self.scan_type = "sec"
            self.sweep_mode = array(['sec']*self.nsweeps)
            if len(ncobj.variables['sweep_start_ray_index']) == 1:
                self.naz = ncobj.variables['sweep_end_ray_index'][0] + 1
            else:
                self.naz = (ncobj.variables['sweep_start_ray_index'][1] -
                            ncobj.variables['sweep_start_ray_index'][0])
            self.nele = ncobj.variables['sweep_start_ray_index'].shape[0]
            self.ngates = len(ncobj.dimensions['range'])
            loc_dict = {}
            for loc_data in ['latitude', 'altitude', 'longitude']:
                loc_dict.update(
                    {loc_data: ncvar_to_field(ncobj.variables[loc_data])})
            self.location = loc_dict
            sweep_dict = {}
            for sweep_data in ['sweep_start_ray_index', 'sweep_mode', 'sweep_number', 'sweep_end_ray_index', 'fixed_angle']:
                sweep_dict.update(
                    {sweep_data: ncvar_to_field(ncobj.variables[sweep_data])})
            self.sweep_info = sweep_dict
            inst_dict = {}
            for inst_data in ['frequency', 'follow_mode', 'pulse_width', 'prt_mode', 'prt', 'prt_ratio', 'polarization_mode', 'nyquist_velocity', 'unambiguous_range', 'n_samples']:
                if inst_data in ncobj.variables.keys():
                    inst_dict.update(
                        {inst_data: (
                            ncvar_to_field(ncobj.variables[inst_data]))})
            self.inst_params = inst_dict
            self.azimuth = ncvar_to_field(ncobj.variables['azimuth'])
            self.range = ncvar_to_field(ncobj.variables['range'])
            self.elevation = ncvar_to_field(ncobj.variables['elevation'])
            self.time = ncvar_to_field(ncobj.variables['time'])
            data_fields = create_field_list(
                ncobj.variables, len(ncobj.dimensions['time']), self.ngates)
            field_dict = {}
            for field in data_fields:
                print field
                my_field = ncvar_to_field(ncobj.variables[field])
                field_dict.update({field: my_field})
            self.fields = field_dict
        if "rhi" in mode:
            #rhi
            self.metadata = dict([(key, getattr(ncobj, key))
                                 for key in ncobj.ncattrs()])
            self.scan_type = "rhi"
            self.nsweeps = len(ncobj.variables['sweep_start_ray_index'])
            self.sweep_mode = array(['rhi']*self.nsweeps)
            if len(ncobj.variables['sweep_start_ray_index']) == 1:
                self.nele = ncobj.variables['sweep_end_ray_index'][0] + 1
            else:
                self.nele = (ncobj.variables['sweep_start_ray_index'][1] -
                             ncobj.variables['sweep_start_ray_index'][0])
            self.naz = ncobj.variables['sweep_start_ray_index'].shape[0]
            self.ngates = len(ncobj.dimensions['range'])
            loc_dict = {}
            for loc_data in ['latitude', 'altitude', 'longitude']:
                loc_dict.update(
                    {loc_data: ncvar_to_field(ncobj.variables[loc_data])})
            self.location = loc_dict
            sweep_dict = {}
            for sweep_data in ['sweep_start_ray_index', 'sweep_mode', 'sweep_number', 'sweep_end_ray_index', 'fixed_angle']:
                sweep_dict.update(
                    {sweep_data: ncvar_to_field(ncobj.variables[sweep_data])})
            self.sweep_info = sweep_dict
            inst_dict = {}
            for inst_data in ['frequency', 'follow_mode', 'pulse_width', 'prt_mode', 'prt', 'prt_ratio', 'polarization_mode', 'nyquist_velocity', 'unambiguous_range', 'n_samples']:
                if inst_data in ncobj.variables.keys():
                    inst_dict.update(
                        {inst_data: (
                            ncvar_to_field(ncobj.variables[inst_data]))})
            self.inst_params = inst_dict
            self.azimuth = ncvar_to_field(ncobj.variables['azimuth'])
            self.range = ncvar_to_field(ncobj.variables['range'])
            self.elevation = ncvar_to_field(ncobj.variables['elevation'])
            self.time = ncvar_to_field(ncobj.variables['time'])
            data_fields = create_field_list(
                ncobj.variables, len(ncobj.dimensions['time']), self.ngates)
            field_dict = {}
            for field in data_fields:
                print field
                my_field = ncvar_to_field(ncobj.variables[field])
                field_dict.update({field: my_field})
            self.fields = field_dict

    def streamcf2rad(self, ncobj):
        try:
            mode = "".join(ncobj.variables['sweep_mode'][0])
        except TypeError:
            mode = "".join(ncobj.variables['sweep_mode'][0].data)
        print mode, "azimuth_surveillance    "
        if mode in "azimuth_surveillance    ":
            #ppi
            print "hi"
            self.metadata = dict(
                [(key, getattr(ncobj, key)) for key in ncobj.ncattrs()])
            self.scan_type = "ppi"
            self.naz = (ncobj.variables['sweep_start_ray_index'][1] -
                        ncobj.variables['sweep_start_ray_index'][0])
            self.nele = ncobj.variables['sweep_start_ray_index'].shape[0]
            self.ngates = ncobj.variables['range'].shape[0]
            loc_dict = {}
            for loc_data in ['latitude', 'altitude', 'longitude']:
                loc_dict.update(
                    {loc_data: ncvar_to_field(ncobj.variables[loc_data])})
            self.location = loc_dict
            sweep_dict = {}
            for sweep_data in ['sweep_start_ray_index', 'sweep_mode', 'sweep_number', 'sweep_end_ray_index', 'fixed_angle']:
                sweep_dict.update(
                    {sweep_data: ncvar_to_field(ncobj.variables[sweep_data])})
            self.sweep_info = sweep_dict
            self.azimuth = ncvar_to_field(ncobj.variables['azimuth'])
            self.range = ncvar_to_field(ncobj.variables['range'])
            self.elevation = ncvar_to_field(ncobj.variables['elevation'])
            self.time = ncvar_to_field(ncobj.variables['time'])
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
            self.fields = field_dict

    def get_mdv_meta(self, radarobj, field):
        debug = True
        print "go"
        csapr_names = csapr_standard_names()
        moment_fixes = {
            'DBZ_F': {
                'units': 'dBZ',
                'standard_name': 'equivalent_reflectivity_factor',
                'long_name': 'equivalent_reflectivity_factor',
                'valid_max': 80.0,
                'valid_min': -45.0,
                'least_significant_digit': 2},

            'VEL_F': {
                'units': 'm/s',
                'standard_name': (
                    'radial_velocity_of_scatterers_away_from_instrument'),
                'long_name': (
                    'radial_velocity_of_scatterers_away_from_instrument'),
                'valid_max': 95.0,
                'valid_min': -95.0,
                'least_significant_digit': 2},

            'KDP_F': {
                'units': 'degrees/km',
                'standard_name': 'specific_differential_phase_hv',
                'long_name': 'specific_differential_phase_hv',
                'valid_max': 20.0,
                'valid_min': -10.0,
                'least_significant_digit': 2},

            'ZDR_F': {
                'units': 'dB',
                'standard_name': 'log_differential_reflectivity_hv',
                'long_name': 'log_differential_reflectivity_hv',
                'valid_max': 8.0,
                'valid_min': -6.0,
                'least_significant_digit': 3},

            'RHOHV_F': {
                'units': 'ratio',
                'standard_name': 'cross_correlation_ratio_hv',
                'long_name': 'cross_correlation_ratio_hv',
                'valid_max': 1.0,
                'valid_min': 0.0,
                'least_significant_digit': 5},

            'NCP_F': {
                'units': 'ratio',
                'standard_name': 'signal_quality',
                'long_name': 'signal_quality',
                'valid_max': 1.0,
                'valid_min': 0.0,
                'comment': 'Also know as Normalized Coherent Power',
                'least_significant_digit': 5},

            'WIDTH_F': {
                'units': 'm/s',
                'standard_name': 'spectrum_width',
                'long_name': 'spectrum_width',
                'valid_max': 45.0,
                'valid_min': 0.0,
                'least_significant_digit': 2},

            'PHIDP_F': {
                'units': 'degrees',
                'standard_name': 'differential_phase_hv',
                'long_name': 'differential_phase_hv',
                'valid_max': 180.0,
                'valid_min': -180.0,
                'least_significant_digit': 2},

            'VEL_COR': {
                'units': 'm/s',
                'standard_name': (
                    'radial_velocity_of_scatterers_away_from_instrument'),
                'long_name': (
                    'radial_velocity_of_scatterers_away_from_instrument'),
                'valid_max': 45.0,
                'valid_min': -45.0,
                'least_significant_digit': 2},

            'PHIDP_UNF': {
                'units': 'degrees',
                'standard_name': 'differential_phase_hv',
                'long_name': 'differential_phase_hv',
                'valid_max': 480.0,
                'valid_min': 0.0,
                'least_significant_digit': 2},

            'DBZ_AC': {
                'units': 'dBZ',
                'standard_name': 'equivalent_reflectivity_factor',
                'long_name': 'equivalent_reflectivity_factor',
                'valid_max': 80.0,
                'valid_min': -45.0,
                'least_significant_digit': 2},

            'KDP_SOB': {
                'units': 'degrees/km',
                'standard_name': 'specific_differential_phase_hv',
                'long_name': 'specific_differential_phase_hv',
                'valid_max': 20.0,
                'valid_min': -1.0,
                'least_significant_digit': 3}}
        return moment_fixes[field]


def join_radar(radar1, radar2):

    #must have same gate spacing
    new_radar = copy.deepcopy(radar1)
    new_radar.azimuth['data'] = append(radar1.azimuth['data'],
                                       radar2.azimuth['data'])
    new_radar.elevation['data'] = append(radar1.elevation['data'],
                                         radar2.elevation['data'])

    if len(radar1.range['data']) >= len(radar2.range['data']):
        new_radar.range['data'] = radar1.range['data']
    else:
        new_radar.range['data'] = radar3.range['data']
    new_radar.time['data'] = append(radar1.time['data'], radar2.time['data'])

    for var in new_radar.fields.keys():
        sh1 = radar1.fields[var]['data'].shape
        sh2 = radar2.fields[var]['data'].shape
        print sh1, sh2
        new_field = ma.zeros([sh1[0] + sh2[0], max([sh1[1], sh2[1]])]) - 9999.0
        new_field[0:sh1[0], 0:sh1[1]] = radar1.fields[var]['data']
        new_field[sh1[0]:, 0:sh2[1]] = radar2.fields[var]['data']
        new_radar.fields[var]['data'] = new_field

    # This will not work for two already moving platforms..
    # need to enhance later
    if radar1.location['latitude']['data'] != radar2.location['latitude']['data'] or radar1.location['longitude']['data'] != radar2.location['longitude']['data'] or radar1.location['altitude']['data'] != radar2.location['altitude']['data']:
        for key in radar1.location.keys():
            new_radar.location[key]['data'] = append(
                zeros(len(radar1.time['data'])) +
                radar1.location[key]['data'],
                zeros(len(radar2.time['data'])) +
                radar2.location[key]['data'])
    return new_radar
