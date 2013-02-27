""" Python wrapper around the RSL library """

from datetime import datetime

import numpy as np

import _rsl
from radar import Radar
from common import dms_to_d, csapr_standard_names, get_mdv_meta


def read_rsl(filename, add_meta={}):
    """
    Read a file supported by RSL

    Parameters
    ----------
    filename : str
        Name of file whose format is supported by RSL
    add_meta : dict
        Dictionary containing additional metadata to add to the created
        Radar object.  This will overwrite metadata extracted from the file.

    Returns
    -------
    radar : Radar
        Radar object.

    """
    radarobj = _rsl.RSL_anyformat_to_radar(filename)

    # TODO
    # An issue that needs to be resolved is that this code likes all
    # sweeps to have the same number of rays.. so for now we take
    # min(nrays) across sweeps and drop rays out side of this...
    # this is an "easy" issue to resolve caused by the fact I have been
    # treating things as cubes and then flattening them
    # what needs to be done is to make the field['data'] be masked arrays
    # and mask out location where the ray is Null

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
        _rsl.fieldTypes().list.index(available_fields[0])]
    sample_sweep = sample_volume.sweeps[0]
    sample_ray = sample_sweep.rays[0]

    # determine the shape parameters of the fields
    nsweeps = sample_volume.h.nsweeps
    rays = np.array([sample_volume.sweeps[i].h.nrays for i in range(nsweeps)])
    nrays = rays.min()  # see TODO above
    ngates = sample_ray.h.nbins

    # set scan_type, naz, and nele
    if sample_sweep.h.azimuth == -999.0:
        scan_type = 'ppi'
        naz = nrays
        nele = nsweeps
    else:
        scan_type = 'rhi'
        naz = nsweeps
        nele = nrays

    # extract the elevation and azimuth attributes
    azimuth, elevation = extract_rsl_pointing(sample_volume, nsweeps, nrays)

    # the range array which describes the range of all beams
    _range = {
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
    azimuth = {
        'data': azimuth.flatten(),
        'units': 'degrees',
        'comment': 'Azimuth of antenna relative to true north',
        'long_name': 'azimuth_angle_from_true_north',
        'standard_name': 'beam_azimuth_angle', }

    elevation = {
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
    tu = "seconds since " + t_start.strftime("%Y-%m-%dT%H:%M:%S.0Z")
    cal = "gregorian"
    time = {
        'data': np.linspace(0, t_span, nrays * nsweeps),
        'units': tu,
        'calendar': cal,
        'comment': ('Coordinate variable for time. '
                    'Time at the center of each ray, in fractional seconds '
                    'since the global variable time_coverage_start'),
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
            _rsl.fieldTypes().list.index(rsl_field)]
        data = create_cube_array_lim(volume, nsweeps, nrays)
        data[np.where(np.isnan(data))] = -9999.0
        data[np.where(data == 131072)] = -9999.0
        meta = get_mdv_meta(radarobj, field)  # fetch metadata
        fielddict = {'data': np.ma.masked_equal(data, -9999.0).reshape(
            data.shape[0] * data.shape[1], data.shape[2]),
            '_FillValue': -9999.0}
        fielddict.update(meta)
        fields_dict.update({csapr_standard_names()[field]: fielddict})
    fields = fields_dict

    # set the sweep parameters
    if scan_type == 'ppi':
        nsweeps = nele
        sweep_number = {'data': range(nsweeps), 'units': 'count',
                        'long_name': 'sweep_number'}
        sweep_mode = {
            'data': nsweeps * ['azimuth_surveillance    '],
            'long_name': 'sweep_mode',
            'units': 'uniteless',
            'comment': ('Options are: "sector", "coplane", "rhi", '
                        '"vertical_pointing", "idle", '
                        '"azimuth_surveillance", "elevation_surveillance", '
                        '"sunscan", "pointing", "manual_ppi", "manual_rhi"')}
        fixed_angle = {
            'data': np.array([sample_volume.sweeps[i].h.elev
                              for i in range(nsweeps)]),
            'long_name': 'target_angle_for_sweep',
            'units': 'degrees',
            'standard_name': 'target_fixed_angle'}
        sweep_start_ray_index = {
            'data': np.arange(0, len(time['data']), naz),
            'long_name': 'index of first ray in sweep, 0-based',
            'units': 'count'}
        sweep_end_ray_index = {
            'data': np.arange(naz - 1, len(time['data']), naz),
            'long_name': 'index of last ray in sweep, 0-based',
            'units': 'count'}
    elif scan_type == 'rhi':
        print "WASSL!"
        nsweeps = naz
        sweep_number = {'data': range(nsweeps), 'units': 'count',
                        'long_name': 'sweep_number'}
        sweep_mode = {
            'data': nsweeps * ['rhi                     '],
            'long_name': 'sweep_mode',
            'units': 'uniteless',
            'comment': ('Options are: "sector", "coplane", "rhi", '
                        '"vertical_pointing", "idle", '
                        '"azimuth_surveillance", "elevation_surveillance", '
                        '"sunscan", "pointing", "manual_ppi", "manual_rhi"')}
        fixed_angle = {
            'data': np.array([sample_volume.sweeps[i].h.azimuth
                              for i in range(nsweeps)]),
            'long_name': 'target_angle_for_sweep',
            'units': 'degrees',
            'standard_name': 'target_fixed_angle'}
        sweep_start_ray_index = {
            'data': np.arange(0, len(time['data']), nele),
            'long_name': 'index of first ray in sweep, 0-based',
            'units': 'count'}
        sweep_end_ray_index = {
            'data': np.arange(nele - 1, len(time['data']), nele),
            'long_name': 'index of last ray in sweep, 0-based',
            'units': 'count'}

    sweep_info = {
        'sweep_number': sweep_number,
        'sweep_mode': sweep_mode,
        'fixed_angle': fixed_angle,
        'sweep_start_ray_index': sweep_start_ray_index,
        'sweep_end_ray_index': sweep_end_ray_index}
    sweep_mode = np.array([scan_type] * nsweeps)
    sweep_number = np.linspace(0, nsweeps-1, nsweeps)
    metadata = {'original_container': 'rsl'}
    need_from_rsl_header = {
        'name': 'instrument_name', 'project': 'project', 'state': 'state',
        'country': 'country'}
    # rsl_name: radar meta name
    rsl_dict = rsl_header_to_dict(radarobj.contents.h)
    for rsl_key in need_from_rsl_header.keys():
        metadata.update({need_from_rsl_header[rsl_key]: rsl_dict[rsl_key]})
    metadata.update(add_meta)
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
    location = {'latitude': lat, 'longitude': lon, 'altitude': elv}

    # set instrument parameters attribute
    valid_nyq_vel = (abs(sample_ray.nyq_vel) > 0.1)
    inst_params = rsl_extract_inst_param(
        sample_volume, nsweeps, nrays, valid_nyq_vel)

    return Radar(nsweeps, nrays, ngates, scan_type, naz, nele, _range,
                 azimuth, elevation, tu, cal, time, fields, sweep_info,
                 sweep_mode, sweep_number, location, inst_params, metadata)

###############################################
# RSL functions                               #
###############################################


def get_avail_moments(volumes):
    """ Return a list of fields present in a RSL Volumes object """
    av = []
    for i in range(len(volumes)):
        if volumes[i] is not None:
            av.append(_rsl.fieldTypes().list[i])
    return av


def create_cube_array(volume):
    ppi = np.zeros([volume.h.nsweeps, volume.sweeps[0].h.nrays,
                    volume.sweeps[0].rays[0].h.nbins],
                   dtype=np.float32) + 1.31072000e+05
    for levnum in range(volume.h.nsweeps):
        for raynum in range(volume.sweeps[0].h.nrays):
            data = volume.sweeps[levnum].rays[raynum].data
            ppi[levnum, raynum, 0:len(data)] = data
    return ppi


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
    ppi = np.zeros([nsweeps, nrays, volume.sweeps[0].rays[0].h.nbins],
                   dtype='float32') + 1.31072000e+05
    for levnum in range(nsweeps):
        rays = volume.sweeps[levnum].rays
        len_d = len(rays[0].data)
        for raynum in range(nrays):
            ppi[levnum, raynum, 0:len(rays[raynum].data)] = rays[raynum].data
    return ppi


def rsl_extract_inst_param(sample_volume, nsweeps, nrays, valid_nyq_vel):
    """ Extract instrument parameters from RSL Volume """
    total_rays = nsweeps * nrays
    pm_data = np.empty(nsweeps, dtype='|S24')
    nv_data = np.empty(total_rays, dtype='float64')
    pr_data = np.empty(total_rays, dtype='float64')
    ur_data = np.empty(total_rays, dtype='float64')

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
            'comments': ('Pulsing mode Options are: "fixed", "staggered", '
                         '"dual". Assumed "fixed" if missing.')},

        'nyquist_velocity': {
            'data': nv_data,
            'units': 'm/s',
            'comments': "unamb velocity"},

        'prt': {
            'data': pr_data,
            'units': 'seconds',
            'comments': ("Pulse repetition time. For staggered prt, "
                         "also see prt_ratio.")},

        'unambiguous_range': {
            'data': ur_data,
            'units': 'meters',
            'comment': 'Unambiguous range'}}


def prtmode(h):
    # TODO prt mode: Need to fix this.. assumes dual if two prts
    if h.prf2 != h.prf:
        mode = 'dual                    '
    else:
        mode = 'fixed                   '
    return mode


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
    azimuth = np.zeros([nsweeps, nrays], dtype=float)
    elevation = np.zeros([nsweeps, nrays], dtype=float)
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
