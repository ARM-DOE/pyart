"""
pyart.io.rsl
============

Python wrapper around the RSL library.

.. autosummary::
    :toctree: generated/

    read_rsl
    get_avail_moments
    create_cube_array_lim
    rsl_extract_inst_param
    prtmode
    extract_rsl_pointing
    rsl_header_to_dict
    ray_header_time_to_datetime

"""

from datetime import datetime

import numpy as np

import _rsl
from radar import Radar
from common import dms_to_d, COMMON2STANDARD, get_metadata, make_tu_str


def read_rsl(filename, add_meta=None):
    """
    Read a file supported by RSL

    Parameters
    ----------
    filename : str or RSL_radar
        Name of file whose format is supported by RSL, or a RSL_Radar object.
    add_meta : dict or None
        Dictionary containing additional metadata to add to the created
        Radar object.  This will overwrite metadata extracted from the file.
        None is add no additional metadata.

    Returns
    -------
    radar : Radar
        Radar object.

    """
    if type(filename) != str:
        radarobj = filename
    else:
        radarobj = _rsl.RSL_anyformat_to_radar(filename)

    # TODO
    # An issue that needs to be resolved is that this code likes all
    # sweeps to have the same number of rays.. so for now we take
    # min(nrays) across sweeps and drop rays out side of this...
    # this is an "easy" issue to resolve caused by the fact I have been
    # treating things as cubes and then flattening them
    # what needs to be done is to make the field['data'] be masked arrays
    # and mask out location where the ray is Null

    # extract a sample volume, sweep and ray
    available_fields = get_avail_moments(radarobj.contents.volumes)
    first_field_idx = _rsl.fieldTypes().list.index(available_fields[0])
    sample_volume = radarobj.contents.volumes[first_field_idx]
    sample_sweep = sample_volume.sweeps[0]
    sample_ray = sample_sweep.rays[0]

    # determine the shape parameters of the fields
    nsweeps = sample_volume.h.nsweeps
    rays = np.array([sample_volume.sweeps[i].h.nrays for i in range(nsweeps)])
    nrays = rays.min()  # see TODO above
    ngates = sample_ray.h.nbins

    # scan_type, naz, and nele
    if sample_sweep.h.azimuth == -999.0:
        scan_type = 'ppi'
        naz = nrays
        nele = nsweeps
    else:
        scan_type = 'rhi'
        naz = nsweeps
        nele = nrays

    # range
    _range = get_metadata('range')
    _range['data'] = sample_ray.dists
    _range['spacing_is_constant'] = 'true'
    _range['meters_to_center_of_first_gate'] = sample_ray.h.range_bin1
    _range['meters_between_gates'] = sample_ray.h.gate_size

    # azimuth and elevation
    _azimuth, _elevation = extract_rsl_pointing(sample_volume, nsweeps, nrays)
    azimuth = get_metadata('azimuth')
    azimuth['data'] = _azimuth.flatten()
    elevation = get_metadata('elevation')
    elevation['data'] = _elevation.flatten()

    # time
    last_ray = sample_volume.sweeps[-1].rays[-1]
    t_start = ray_header_time_to_datetime(sample_ray.h)
    t_end = ray_header_time_to_datetime(last_ray.h)
    t_span = (t_end - t_start).seconds
    tu = make_tu_str(t_start)
    cal = "gregorian"
    time = get_metadata('time')
    time['data'] = np.linspace(0, t_span, nrays * nsweeps)
    time['units'] = tu
    time['calendar'] = cal

    # fields
    fields = {}
    rsl2common = {'ZT': 'DBZ',
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
    common2rsl = dict((v, k) for k, v in rsl2common.iteritems())

    available_fields_common = [rsl2common[key] for key in available_fields]
    # transfer only those which are available and have a standard name
    for field in set(available_fields_common) & set(COMMON2STANDARD.keys()):

        # extract the field, mask and reshape
        rsl_field = common2rsl[field]
        volume = radarobj.contents.volumes[
            _rsl.fieldTypes().list.index(rsl_field)]
        data = create_cube_array_lim(volume, nsweeps, nrays)
        data[np.where(np.isnan(data))] = -9999.0
        data[np.where(data == 131072)] = -9999.0
        data = np.ma.masked_equal(data, -9999.0)
        data.shape = (data.shape[0] * data.shape[1], data.shape[2])

        # create the field dictionary
        fielddict = get_metadata(field)
        fielddict['data'] = data
        fielddict['_FillValue'] = -9999.0
        fields[COMMON2STANDARD[field]] = fielddict

    # sweep parameters
    sweep_number = get_metadata('sweep_number')
    sweep_mode = get_metadata('sweep_mode')
    fixed_angle = get_metadata('fixed_angle')
    sweep_start_ray_index = get_metadata('sweep_start_ray_index')
    sweep_end_ray_index = get_metadata('sweep_end_ray_index')
    len_time = len(time['data'])

    if scan_type == 'ppi':
        nsweeps = nele
        sweep_number['data'] = range(nsweeps)
        sweep_mode['data'] = nsweeps * ['azimuth_surveillance    ']
        fixed_angle['data'] = np.array([sample_volume.sweeps[i].h.elev
                                        for i in range(nsweeps)])
        sweep_start_ray_index['data'] = np.arange(0, len_time, naz)
        sweep_end_ray_index['data'] = np.arange(naz - 1, len_time, naz)

    elif scan_type == 'rhi':
        nsweeps = naz
        sweep_number['data'] = range(nsweeps)
        sweep_mode['data'] = nsweeps * ['rhi                     ']
        fixed_angle['data'] = np.array([sample_volume.sweeps[i].h.azimuth
                                        for i in range(nsweeps)])
        sweep_start_ray_index['data'] = np.arange(0, len_time, nele)
        sweep_end_ray_index['data'] = np.arange(nele - 1, len_time, nele)

    sweep_info = {
        'sweep_number': sweep_number,
        'sweep_mode': sweep_mode,
        'fixed_angle': fixed_angle,
        'sweep_start_ray_index': sweep_start_ray_index,
        'sweep_end_ray_index': sweep_end_ray_index}
    sweep_mode = np.array([scan_type] * nsweeps)
    sweep_number = np.linspace(0, nsweeps - 1, nsweeps)

    # metadata
    metadata = {'original_container': 'rsl'}
    rsl_dict = rsl_header_to_dict(radarobj.contents.h)
    need_from_rsl_header = {
        'name': 'instrument_name', 'project': 'project', 'state': 'state',
        'country': 'country'}  # rsl_name : radar_metadata_name
    for rsl_key, metadata_key in need_from_rsl_header.iteritems():
        metadata[metadata_key] = rsl_dict[rsl_key]
    if add_meta is not None:
        metadata.update(add_meta)

    # location
    lat = get_metadata('latitude')
    lon = get_metadata('longitude')
    elv = get_metadata('altitude')
    latd = radarobj.contents.h.latd
    latm = radarobj.contents.h.latm
    lats = radarobj.contents.h.lats
    lat['data'] = dms_to_d((latd, latm, lats))

    lond = radarobj.contents.h.lond
    lonm = radarobj.contents.h.lonm
    lons = radarobj.contents.h.lons
    lon['data'] = dms_to_d((lond, lonm, lons))

    elv['data'] = radarobj.contents.h.height
    location = {'latitude': lat, 'longitude': lon, 'altitude': elv}

    # set instrument parameters attribute
    valid_nyq_vel = abs(sample_ray.nyq_vel) > 0.1
    inst_params = rsl_extract_inst_param(
        sample_volume, nsweeps, nrays, valid_nyq_vel)

    return Radar(nsweeps, nrays, ngates, scan_type, naz, nele, _range,
                 azimuth, elevation, tu, cal, time, fields, sweep_info,
                 sweep_mode, sweep_number, location, inst_params, metadata)


#####################
# private functions #
#####################

def get_avail_moments(volumes):
    """ Return a list of fields present in a RSL Volumes object """
    av = []
    for i in xrange(len(volumes)):
        if volumes[i] is not None:
            av.append(_rsl.fieldTypes().list[i])
    return av


def create_cube_array_lim(volume, nsweeps, nrays):
    """
    Extract a field from an RSL Volume.

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
    data : array, (nsweep, nrays, nbins), dtype=float32
        Three dimensional array holding the extracted field.

    """
    nbins = volume.sweeps[0].rays[0].h.nbins
    ppi = np.zeros([nsweeps, nrays, nbins], dtype='float32') + 1.31072000e+05
    for levnum in xrange(nsweeps):
        rays = volume.sweeps[levnum].rays
        for raynum in range(nrays):
            ppi[levnum, raynum, 0:len(rays[raynum].data)] = rays[raynum].data
    return ppi


def rsl_extract_inst_param(sample_volume, nsweeps, nrays, valid_nyq_vel):
    """ Return an instrument parameter dictionary from an RSL Volume. """
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
    inst_param = {}
    inst_param['prt_mode'] = get_metadata('prt_mode')
    inst_param['nyquist_velocity'] = get_metadata('nyquist_velocity')
    inst_param['prt'] = get_metadata('prt')
    inst_param['unambiguous_range'] = get_metadata('unambiguous_range')
    inst_param['prt_mode']['data'] = pm_data
    inst_param['nyquist_velocity']['data'] = nv_data
    inst_param['prt']['data'] = pr_data
    inst_param['unambiguous_range']['data'] = ur_data
    return inst_param


def prtmode(h):
    """ Return the prt mode from a RSL ray header. """
    # TODO prt mode: Need to fix this.. assumes dual if two prts
    if h.prf2 != h.prf:
        mode = 'dual                    '
    else:
        mode = 'fixed                   '
    return mode


def extract_rsl_pointing(volume, nsweeps, nrays):
    """
    Extract the azimuth and elevation parameters from a RSL Volume.

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
    azimuth = np.zeros([nsweeps, nrays], dtype=np.float)
    elevation = np.zeros([nsweeps, nrays], dtype=np.float)
    for i, sweep in enumerate(volume.sweeps):
        for j, ray in enumerate(sweep.rays):
            azimuth[i, j] = ray.h.azimuth
            elevation[i, j] = ray.h.elev
    return azimuth, elevation


def rsl_header_to_dict(header):
    """ Return a dictionary containing the attributes in a RSL header. """
    all_keys = dir(header)
    d = {}
    for key in all_keys:
        if key[0] != '_':
            d[key] = getattr(header, key)
    return d


def ray_header_time_to_datetime(h):
    """ Return a datetime object from a RSL ray header. """
    return datetime(h.year, h.month, h.day, h.hour, h.minute, int(h.sec))
