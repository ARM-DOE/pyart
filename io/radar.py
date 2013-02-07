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

import py4dd


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

###############################################
# RSL functions                               #
# should be moved to rsl.py or a similar file #
###############################################


def get_avail_moments(volumes):
    """ Return a list of fields present in a RSL Volumes object """
    av = []
    for i in range(len(volumes)):
        if volumes[i] is not None:
            av.append(py4dd.fieldTypes().list[i])
    return av


# move to rsl.py
def create_cube_array_lim(volume, nsweeps, nrays):
    """ Extract a field from an RSL Volume.

    Parameters
    ----------
    volume : RSL_Volume
        RSL Volume from which to extract the field.
    nsweeps : int
        Number of valid (non-null) sweeps in the volume.
    nrays : int
        Number of valid (non-null rays in each Sweep in the volume.

    Returns
    -------
    data : array, (nsweep, nrays, *), dtype=float32
        Three dimensional array holding the extracted field.

    """
    # XXX move to rsl.py
    ppi = zeros([nsweeps, nrays, volume.sweeps[0].rays[0].h.nbins],
                dtype=float32) + 1.31072000e+05
    for levnum in range(nsweeps):
        rays = volume.sweeps[levnum].rays
        len_d = len(rays[0].data)
        for raynum in range(nrays):
            ppi[levnum, raynum, 0:len_d] = rays[raynum].data
    return ppi


def rsl_extract_inst_param(sample_volume, nsweeps, nrays, valid_nyq_vel):
    """ Extract instrument parameters from RSL Volume """
    # XXX move to rsl.py
    total_rays = nsweeps * nrays
    pm_data = empty(nsweeps, dtype='|S24')
    nv_data = empty(total_rays, dtype='float64')
    pr_data = empty(total_rays, dtype='float64')
    ur_data = empty(total_rays, dtype='float64')

    for i, sweep in enumerate(sample_volume.sweeps):
        for j, ray in enumerate(sweep.rays):
            if j == 0:
                pm_data[i] = prtmode(ray.h)
            idx = j * nsweeps + i

            if valid_nyq_vel:
                nv_data[idx] = ray.nyq_vel
            else:
                nv_data[idx] = ray.wavelength * ray.prf / 4.0
            pr_data[idx] = 1. / ray.prf
            ur_data[idx] = ray.unam_rng * 1000.0

    return {
        'prt_mode': {
            'data': pm_data,
            'comments': 'Pulsing mode Options are: "fixed", "staggered", "dual" Assumed "fixed" if missing.'},

        'nyquist_velocity': {
            'data': nv_data,
            'units': 'm/s',
            'comments': "unamb velocity"},

        'prt': {
            'data': pr_data,
            'units': 'seconds',
            'comments': "Pulse repetition time. For staggered prt, also see prt_ratio."},

        'unambiguous_range': {
            'data': ur_data,
            'units': 'meters',
            'comment': 'Unambiguous range'}}


def prtmode(h):
    # XXX move to rsl.py
    # TODO prt mode: Need to fix this.. assumes dual if two prts
    if h.prf2 != h.prf:
        mode = 'dual                    '
    else:
        mode = 'fixed                   '
    return mode


# XXX move to rsl.py
def extract_rsl_pointing(volume, nsweeps, nrays):
    """ Extract the azimuth and elevation parameters from a RSL Volume.

    Parameters
    ----------
    volume : RSL_Volume
        RSL Volume from which to extract the azimuth and elevation from.
    nsweeps : int
        Number of valid (non-null) sweeps in the volume
    nrays : int
        Number of valid (non-null rays in each Sweep in the volume.

    Returns
    -------
    azimuth : array, (nsweeps, nrays),  dtype=float
        Array containing azimuth values in degrees.
    elevation : array, (nsweeps, nrays), dtype=float
        Array containing elevation values in degrees.
    """
    azimuth = zeros([nsweeps, nrays], dtype=float)
    elevation = zeros([nsweeps, nrays], dtype=float)
    for i, sweep in enumerate(volume.sweeps):
        for j, ray in enumerate(sweep.rays):
            azimuth[i, j] = ray.h.azimuth
            elevation[i, j] = ray.h.elev
    return azimuth, elevation


def rsl_header_to_dict(header):
    all_keys = dir(header)
    my_dict = {}
    for key in all_keys:
        if key[0] != '_':
            my_dict.update({key: getattr(header, key)})
    return my_dict


def ray_header_time_to_dict(h):
    return {'year': h.year, 'month': h.month, 'day': h.day,
            'hour': h.hour, 'minute': h.minute, 'second': h.sec}


def ray_header_time_to_datetime(h):
    return datetime(h.year, h.month, h.day, h.hour, h.minute, int(h.sec))


### END


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

    def rsl2rad(self, radarobj, **kwargs):
        """ Load data and parameters from RSL Radar object"""

        # TODO
        # An issue that needs to be resolved is that this code likes all
        # sweeps to have the same number of rays.. so for now we take
        # min(nrays) across sweeps and drop rays out side of this...
        # this is an "easy" issue to resolve caused by the fact I have been
        # treating things as cubes and then flattening them
        # what needs to be done is to make the field['data'] be masked arrays
        # and mask out location where the ray is Null

        # additional metadata which will overwrite data from the radar header.
        # this helps when you know there are issues with the meta
        add_meta = kwargs.get('add_meta', {})

        # determine which fields should be transfered to the radar object
        # only rtransfer the fields which we have valid names.
        available_fields = get_avail_moments(radarobj.contents.volumes)
        name_transfer = {'ZT': 'DBZ',
                         'VR': 'VEL_F',
                         'DR': 'ZDR',
                         'KD': 'KDP_F',
                         'SQ': 'NCP_F',
                         'PH': 'PHIDP_F',
                         'VE': 'VEL_COR',
                         'RH': 'RHOHV_F',
                         'DZ': 'DBZ_F',
                         'CZ': 'DBZ',
                         'SW': 'WIDTH',
                         'ZD': 'ZDR_F',
                         'TI': 'VEL'}
        name_transfer_back = dict((v, k) for (k, v) in
                                  name_transfer.iteritems())
        # the next line requires python 2.7+ dict comprehensions
        #name_transfer_back = {v: k for v, k in name_transfer.iteritems()}
        fields = [name_transfer[key] for key in available_fields]
        todo_fields = set(fields) & set(csapr_standard_names().keys())

        # extract a sample volume, sweep and ray
        sample_volume = radarobj.contents.volumes[
            py4dd.fieldTypes().list.index(available_fields[0])]
        sample_sweep = sample_volume.sweeps[0]
        sample_ray = sample_sweep.rays[0]

        # determine the shape parameters of the fields
        self.nsweeps = sample_volume.h.nsweeps
        rays = array([sample_volume.sweeps[i].h.nrays for i in
                      range(self.nsweeps)])
        self.nrays = rays.min()  # see TODO above
        self.ngates = sample_ray.h.nbins

        # set scan_type, naz, and nele
        if sample_sweep.h.azimuth == -999.0:
            self.scan_type = 'ppi'
            self.naz = self.nrays
            self.nele = self.nsweeps
        else:
            self.scan_type = 'rhi'
            self.naz = self.nsweeps
            self.nele = self.nrays

        # extract the elevation and azimuth attributes
        azimuth, elevation = extract_rsl_pointing(sample_volume, self.nsweeps,
                                                  self.nrays)

        # the range array which describes the range of all beams
        self.range = {
            'data': sample_ray.dists,
            'units': 'meters',
            'standard_name': 'projection_range_coordinate',
            'long_name': 'range_to_measurement_volume',
            'comment': (
                'Coordinate variable for range. Range to center of each bin.'),
            'spacing_is_constant': 'true',
            'meters_to_center_of_first_gate': sample_ray.h.range_bin1,
            'meters_between_gates': sample_ray.h.gate_size
        }

        # the flat azimuth array which describes the azimuth of each beam
        self.azimuth = {
            'data': azimuth.flatten(),
            'units': 'degrees',
            'comment': 'Azimuth of antenna relative to true north',
            'long_name': 'azimuth_angle_from_true_north',
            'standard_name': 'beam_azimuth_angle', }

        self.elevation = {
            'data': elevation.flatten(),
            'units': 'degrees',
            'standard_name': 'beam_elevation_angle',
            'comment': 'Elevation of antenna relative to the horizontal plane',
            'long_name': 'elevation_angle_from_horizontal_plane', }

        # set the time attributes
        last_ray = sample_volume.sweeps[-1].rays[-1]
        t_start = ray_header_time_to_datetime(sample_ray.h)
        t_end = ray_header_time_to_datetime(last_ray.h)
        t_span = (t_end - t_start).seconds
        self.tu = "seconds since " + t_start.strftime("%Y-%m-%d %H:%M:%S.0")
        self.cal = "gregorian"
        self.time = {
            'data': linspace(0, t_span, self.nrays * self.nsweeps),
            'units': self.tu,
            'calendar': self.cal,
            'comment': 'Coordinate variable for time. Time at the center of each ray, in fractional seconds since the global variable time_coverage_start',
            'standard_name': 'time',
            'long_name': 'time in seconds since volume start'}

        # extract the fields
        fields_dict = {}
        for field in todo_fields:
            #create a dictionary tree for all data fields
            print "Doing ", field
            rsl_field = [key for key, value in name_transfer.iteritems() if
                         value == field][0]
            print "Corresponds to ",  rsl_field

            volume = radarobj.contents.volumes[
                py4dd.fieldTypes().list.index(rsl_field)]
            data = create_cube_array_lim(volume, self.nsweeps, self.nrays)
            data[where(isnan(data))] = -9999.0
            data[where(data == 131072)] = -9999.0
            meta = self.get_mdv_meta(radarobj, field)  # fetch metadata
            fielddict = {'data': ma.masked_equal(data, -9999.0).reshape(
                data.shape[0] * data.shape[1], data.shape[2]),  '_FillValue': -9999.0}
            fielddict.update(meta)
            fields_dict.update({csapr_standard_names()[field]: fielddict})
        self.fields = fields_dict

        # set the sweep parameters
        if self.scan_type == 'ppi':
            self.nsweeps = self.nele
            sweep_number = {'data': range(self.nsweeps), 'units': 'count',
                            'long_name': 'sweep_number'}
            sweep_mode = {
                'data': self.nsweeps * ['azimuth_surveillance    '],
                'long_name': 'sweep_mode',
                'units': 'uniteless',
                'comment': 'Options are:"sector","coplane",rhi","vertical_pointing","idle","azimuth_surveillance","elevation_surveillance","sunscan","pointing","manual_ppi","manual_rhi"'}
            fixed_angle = {
                'data': array([sample_volume.sweeps[i].h.elev
                               for i in range(self.nsweeps)]),
                'long_name': 'target_angle_for_sweep',
                'units': 'degrees',
                'standard_name': 'target_fixed_angle'}
            sweep_start_ray_index = {
                'data': arange(0, len(self.time['data']), self.naz),
                'long_name': 'index of first ray in sweep, 0-based',
                'units': 'count'}
            sweep_end_ray_index = {
                'data': arange(self.naz - 1, len(self.time['data']), self.naz),
                'long_name': 'index of last ray in sweep, 0-based',
                'units': 'count'}
        elif self.scan_type == 'rhi':
            print "WASSL!"
            self.nsweeps = self.naz
            sweep_number = {'data': range(self.nsweeps), 'units': 'count',
                            'long_name': 'sweep_number'}
            sweep_mode = {
                'data': self.nsweeps * ['rhi                     '],
                'long_name': 'sweep_mode',
                'units': 'uniteless',
                'comment': 'Options are:"sector","coplane",rhi","vertical_pointing","idle","azimuth_surveillance","elevation_surveillance","sunscan","pointing","manual_ppi","manual_rhi"'}
            fixed_angle = {
                'data': array([sample_volume.sweeps[i].h.azimuth
                               for i in range(self.nsweeps)]),
                'long_name': 'target_angle_for_sweep',
                'units': 'degrees',
                'standard_name': 'target_fixed_angle'}
            sweep_start_ray_index = {
                'data': arange(0, len(self.time['data']), self.nele),
                'long_name': 'index of first ray in sweep, 0-based',
                'units': 'count'}
            sweep_end_ray_index = {
                'data': arange(self.nele - 1, len(self.time['data']),
                               self.nele),
                'long_name': 'index of last ray in sweep, 0-based',
                'units': 'count'}

        self.sweep_info = {
            'sweep_number': sweep_number,
            'sweep_mode': sweep_mode,
            'fixed_angle': fixed_angle,
            'sweep_start_ray_index': sweep_start_ray_index,
            'sweep_end_ray_index': sweep_end_ray_index}
        self.sweep_mode = array([self.scan_type] * self.nsweeps)
        self.sweep_number = linspace(0, self.nsweeps-1, self.nsweeps)
        metadata = {'original_container': 'rsl'}
        need_from_rsl_header = {
            'name': 'instrument_name', 'project': 'project', 'state': 'state',
            'country': 'country'}
        # rsl_name: radar meta name
        rsl_dict = rsl_header_to_dict(radarobj.contents.h)
        for rsl_key in need_from_rsl_header.keys():
            metadata.update({need_from_rsl_header[rsl_key]: rsl_dict[rsl_key]})
        metadata.update(add_meta)
        self.metadata = metadata
        #now for location variables
        radar_loc = [dms_to_d((
            radarobj.contents.h.latd, radarobj.contents.h.latm,
            radarobj.contents.h.lats)),
            dms_to_d((radarobj.contents.h.lond, radarobj.contents.h.lonm,
                     radarobj.contents.h.lons))]
        lat = {'data': radar_loc[0], 'standard_name': 'Latitude',
               'units': 'degrees_north'}
        lon = {'data': radar_loc[1], 'standard_name': 'Longitude',
               'units': 'degrees_east'}
        elv = {'data': radarobj.contents.h.height,
               'standard_name': 'Altitude', 'units': 'meters'}
        self.location = {'latitude': lat, 'longitude': lon, 'altitude': elv}

        # set instrument parameters attribute
        valid_nyq_vel = (abs(sample_ray.nyq_vel) > 0.1)
        self.inst_params = rsl_extract_inst_param(
            sample_volume, self.nsweeps, self.nrays, valid_nyq_vel)

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

    def mdv2rad(self, radarobj):

        # We only want to transfer fields that we have valid names for...
        valid_fields = csapr_standard_names()
        todo_fields = set(radarobj.fields) & set(csapr_standard_names().keys())
        # the intersection of the set of available and valid names
        flat_dict = {}  # a holder to be updated
        self.naz = len(radarobj.az_deg)
        self.nele = len(radarobj.el_deg)
        self.ngates = len(radarobj.range_km)

        #The flat azimuth array which describes the azimuth of each beam
        self.azimuth = {
            'data': tile(radarobj.az_deg, self.nele),
            'units': 'degrees',
            'comment': 'Azimuth of antenna relative to true north',
            'long_name': 'azimuth_angle_from_true_north',
            'standard_name': 'beam_azimuth_angle'}

        #the range array which describes the range of all beams (note in this
        self.range = {
            'data': array(radarobj.range_km*1000.0),
            'units': 'meters',
            'standard_name': 'projection_range_coordinate',
            'long_name': 'range_to_measurement_volume',
            'comment': (
                'Coordinate variable for range. Range to center of each bin.'),
            'spacing_is_constant': 'true',
            'meters_to_center_of_first_gate': (radarobj.range_km[0]) * 1000.0,
            'meters_between_gates': (radarobj.range_km[1] -
                                     radarobj.range_km[0])*1000.0}

        self.elevation = {
            'data': array(radarobj.el_deg).repeat(self.naz),
            'units': 'degrees',
            'standard_name': 'beam_elevation_angle',
            'comment': 'Elevation of antenna relative to the horizontal plane',
            'long_name': 'elevation_angle_from_horizontal_plane'}
        self.scan_type = radarobj.scan_type

        #append time
        tu = "seconds since %(year)d-%(month)02d-%(day)02d %(hour)02d:%(minute)02d:%(second)02d.0" % dt_to_dict(radarobj.times['time_begin'])
        cal = "gregorian"
        time_array = linspace(date2num(radarobj.times['time_begin'], tu, cal),
                              date2num(radarobj.times['time_end'], tu, cal),
                              self.naz * self.nele)
        self.time = {
            'data': time_array, 'units': tu, 'calendar': cal,
            'comment': 'Coordinate variable for time. Time at the center of each ray, in fractional seconds since the global variable time_coverage_start',
            'standard_name': 'time',
            'long_name': 'time in seconds since volume start'}
        for field in todo_fields:
            # create a dictionary tree for all data fields
            print "Doing ", field
            # grab data from MDV object
            data = radarobj.read_a_field(radarobj.fields.index(field))
            data[where(isnan(data))] = -9999.0
            data[where(data == 131072)] = -9999.0
            meta = self.get_mdv_meta(radarobj, field)  # fetch metadata
            fielddict = {
                'data': ma.masked_equal(data, -9999.0).reshape(
                    data.shape[0] * data.shape[1], data.shape[2]),
                '_FillValue': -9999.0}
            fielddict.update(meta)
            flat_dict.update({csapr_standard_names()[field]: fielddict})
        self.fields = flat_dict
        #sweep based stuff now:
        if radarobj.scan_type == 'ppi':
            self.nsweeps = self.nele
            sweep_number = {'data': range(self.nsweeps), 'units': 'count',
                            'long_name': 'sweep_number'}
            sweep_mode = {
                'data': self.nsweeps*['azimuth_surveillance    '],
                'long_name': 'sweep_mode',
                'units': 'uniteless',
                'comment': 'Options are:"sector","coplane",rhi","vertical_pointing","idle","azimuth_surveillance","elevation_surveillance","sunscan","pointing","manual_ppi","manual_rhi"'}
            fixed_angle = {
                'data': array(radarobj.el_deg),
                'long_name': 'target_angle_for_sweep',
                'units': 'degrees',
                'standard_name': 'target_fixed_angle'}
            sweep_start_ray_index = {
                'data': arange(0, len(self.time['data']), self.naz),
                'long_name': 'index of first ray in sweep, 0-based',
                'units': 'count'}
            sweep_end_ray_index = {
                'data': arange(self.naz-1, len(self.time['data']), self.naz),
                'long_name': 'index of last ray in sweep, 0-based',
                'units': 'count'}
        elif radarobj.scan_type == 'rhi':
            self.nsweeps = self.naz
            sweep_number = {
                'data': range(self.nsweeps),
                'units': 'count',
                'long_name': 'sweep_number'}
            sweep_mode = {
                'data': self.nsweeps * ['rhi                     '],
                'long_name': 'sweep_mode',
                'units': 'uniteless',
                'comment': 'Options are:"sector","coplane",rhi","vertical_pointing","idle","azimuth_surveillance","elevation_surveillance","sunscan","pointing","manual_ppi","manual_rhi"'}
            fixed_angle = {
                'data': array(radarobj.az_deg),
                'long_name': 'target_angle_for_sweep',
                'units': 'degrees',
                'standard_name': 'target_fixed_angle'}
            sweep_start_ray_index = {
                'data': arange(0, len(self.time['data']), self.nele),
                'long_name': 'index of first ray in sweep, 0-based',
                'units': 'count'}
            sweep_end_ray_index = {
                'data': arange(self.nele - 1, len(self.time['data']),
                               self.nele),
                'long_name': 'index of last ray in sweep, 0-based',
                'units': 'count'}

        elif radarobj.scan_type == 'vpr':
            self.nsweeps = 1
        self.sweep_info = {
            'sweep_number': sweep_number,
            'sweep_mode': sweep_mode,
            'fixed_angle': fixed_angle,
            'sweep_start_ray_index': sweep_start_ray_index,
            'sweep_end_ray_index': sweep_end_ray_index}
        self.sweep_mode = array([radarobj.scan_type] * self.nsweeps)
        self.sweep_number = linspace(0, self.nsweeps-1, self.nsweeps)
        mapme = defaut_mdv_metadata_map()
        metadata = {}
        masterdata = dict([(key, radarobj.master_header[mapme[key]])
                           for key in mapme.keys()])
        metadata.update(masterdata)
        self.metadata = metadata
        #now for location variables
        lat = {
            'data': radarobj.radar_info['latitude_deg'],
            'standard_name': 'Latitude',
            'units': 'degrees_north'}
        lon = {
            'data': radarobj.radar_info['longitude_deg'],
            'standard_name': 'Longitude',
            'units': 'degrees_east'}
        elv = {
            'data': radarobj.radar_info['altitude_km']*1000.0,
            'standard_name': 'Altitude',
            'units': 'meters'}
        self.location = {'latitude': lat, 'longitude': lon, 'altitude': elv}
        # now for instrument parameters..
        # sorry but I am just going to brute force this!
        # prt mode: Need to fix this.. assumes dual if two prts
        if radarobj.radar_info['prt2_s'] == 0.0:
            prt_mode = 'fixed                   '
        else:
            prt_mode = 'dual                    '
        inst_params = {
            'prt_mode': {
                'data': array([prt_mode]*self.nele),
                'comments': 'Pulsing mode Options are: "fixed", "staggered", "dual" Assumed "fixed" if missing.'},

            'prt': {
                'data': array([radarobj.radar_info['prt_s']] *
                              self.nele * self.naz),
                'units': 'seconds',
                'comments': "Pulse repetition time.For staggered prt, also see prt_ratio."},

            'unambiguous_range': {
                'data': array([radarobj.radar_info['unambig_range_km'] *
                              1000.0] * self.naz * self.nele),
                'units': 'meters',
                'comment': 'Unambiguous range'},

            'nyquist_velocity': {
                'data': array([radarobj.radar_info['unambig_vel_mps']] *
                              self.naz*self.nele),
                'units': 'm/s',
                'comments': "unamb velocity"}
        }
        self.inst_params = inst_params

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
