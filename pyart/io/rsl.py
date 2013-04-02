"""
pyart.io.rsl
============

Python wrapper around the RSL library.

.. autosummary::
    :toctree: generated/

    read_rsl
    create_cube_array_lim
    ray_header_time_to_datetime


"""

from datetime import datetime

import numpy as np

import _rsl_interface
from pyart.io.radar import Radar
from pyart.io.common import dms_to_d, get_metadata, make_tu_str


def read_rsl(filename, add_meta=None):
    """
    Read a file supported by RSL

    Parameters
    ----------
    filename : str or RSL_radar
        Name of file whose format is supported by RSL.
    add_meta : dict or None
        Dictionary containing additional metadata to add to the created
        Radar object.  This will overwrite metadata extracted from the file.
        None is add no additional metadata.

    Returns
    -------
    radar : Radar
        Radar object.

    """
    rslfile = _rsl_interface.RslFile(filename)
    available_vols = rslfile.available_moments()
    first_volume = rslfile.get_volume(available_vols[0])
    first_sweep = first_volume.get_sweep(0)
    first_ray = first_sweep.get_ray(0)

    # TODO
    # An issue that needs to be resolved is that this code likes all
    # sweeps to have the same number of rays.. so for now we take
    # min(nrays) across sweeps and drop rays out side of this...
    # this is an "easy" issue to resolve caused by the fact I have been
    # treating things as cubes and then flattening them
    # what needs to be done is to make the field['data'] be masked arrays
    # and mask out location where the ray is Null

    # determine the shape parameters of the fields
    nsweeps = first_volume.nsweeps
    nrays = min(first_volume.get_nray_list())
    ngates = first_ray.nbins

    # scan_type, naz, and nele
    if first_sweep.azimuth == -999.0:
        scan_type = 'ppi'
        naz = nrays
        nele = nsweeps
    else:
        scan_type = 'rhi'
        naz = nsweeps
        nele = nrays

    # range
    _range = get_metadata('range')
    range_bin1 = first_ray.range_bin1
    gate_size = first_ray.gate_size
    _range['data'] = range_bin1 + gate_size * np.arange(ngates)
    _range['spacing_is_constant'] = 'true'
    _range['meters_to_center_of_first_gate'] = range_bin1
    _range['meters_between_gates'] = gate_size

    # azimuth and elevation
    _azimuth, _elevation = first_volume.get_azimuth_and_elev_array()
    azimuth = get_metadata('azimuth')
    azimuth['data'] = _azimuth.flatten()
    elevation = get_metadata('elevation')
    elevation['data'] = _elevation.flatten()

    # time
    last_sweep = first_volume.get_sweep(nsweeps - 1)
    last_ray = last_sweep.get_ray(nrays - 1)
    t_start = first_ray.get_datetime()
    t_end = last_ray.get_datetime()
    t_span = (t_end - t_start).seconds
    tu = make_tu_str(t_start)
    cal = "gregorian"
    time = get_metadata('time')
    time['data'] = np.linspace(0, t_span, nrays * nsweeps)
    time['units'] = tu
    time['calendar'] = cal

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
        fixed_angle['data'] = first_volume.get_sweep_elevs()
        sweep_start_ray_index['data'] = np.arange(0, len_time, naz)
        sweep_end_ray_index['data'] = np.arange(naz - 1, len_time, naz)

    elif scan_type == 'rhi':
        nsweeps = naz
        sweep_number['data'] = range(nsweeps)
        sweep_mode['data'] = nsweeps * ['rhi                     ']
        fixed_angle['data'] = first_volume.get_sweep_azimuths()
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
    rsl_dict = rslfile.get_radar_header()
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
    latd = rsl_dict['latd']
    latm = rsl_dict['latm']
    lats = rsl_dict['lats']
    lat['data'] = dms_to_d((latd, latm, lats))

    lond = rsl_dict['lond']
    lonm = rsl_dict['lonm']
    lons = rsl_dict['lons']
    lon['data'] = dms_to_d((lond, lonm, lons))

    elv['data'] = rsl_dict['height']
    location = {'latitude': lat, 'longitude': lon, 'altitude': elv}

    # set instrument parameters attribute
    inst_params = {}
    inst_params['prt_mode'] = get_metadata('prt_mode')
    inst_params['nyquist_velocity'] = get_metadata('nyquist_velocity')
    inst_params['prt'] = get_metadata('prt')
    inst_params['unambiguous_range'] = get_metadata('unambiguous_range')

    pm_data, nv_data, pr_data, ur_data = first_volume.get_instr_params()
    inst_params['prt_mode']['data'] = pm_data
    inst_params['nyquist_velocity']['data'] = nv_data.flatten()
    inst_params['prt']['data'] = pr_data.flatten()
    inst_params['unambiguous_range']['data'] = ur_data.flatten()

    # fields
    # transfer only those which are available and have a standard name
    fields = {}
    available_vols = rslfile.available_moments()
    good_vols = VOLUMENUM2STANDARDNAME.keys()
    volumes_to_extract = [i for i in available_vols if i in good_vols]

    for volume_num in volumes_to_extract:

        # extract the field, mask and reshape
        data = rslfile.get_volume_array(volume_num)
        data[np.where(np.isnan(data))] = -9999.0
        data[np.where(data == 131072)] = -9999.0
        data = np.ma.masked_equal(data, -9999.0)
        data.shape = (data.shape[0] * data.shape[1], data.shape[2])

        # create the field dictionary
        standard_field_name = VOLUMENUM2STANDARDNAME[volume_num]
        fielddict = get_metadata(standard_field_name)
        fielddict['data'] = data
        fielddict['_FillValue'] = -9999.0
        fields[standard_field_name] = fielddict

    return Radar(nsweeps, nrays, ngates, scan_type, naz, nele, _range,
                 azimuth, elevation, tu, cal, time, fields, sweep_info,
                 sweep_mode, sweep_number, location, inst_params, metadata)


#####################
# private functions #
#####################

#   id  abbreviation    common_name     standard_name
#   0   DZ              DBZ_F           reflectivity_horizontal
#   1   VR              VEL_F           mean_doppler_velocity
#   2   SW              WIDTH
#   3   CZ              DBZ
#   4   ZT              DBZ
#   5   DR              ZDR
#   6   LR
#   7   ZD              ZDR_F           diff_reflectivity
#   8   DM
#   9   RH              RHOHV_F         copol_coeff
#   10  PH              PHIDP_F         dp_phase_shift
#   11  XZ
#   12  CD
#   13  MZ
#   14  MD
#   15  ZE
#   16  VE              VEL_COR         corrected_mean_doppler_velocity
#   17  KD              KDP_F           diff_phase
#   18  TI              VEL
#   19  DX
#   20  CH
#   21  AH
#   22  CV
#   23  AV
#   24  SQ              NCP_F           norm_coherent_power

VOLUMENUM2STANDARDNAME = {
    0: 'reflectivity_horizontal',
    1: 'mean_doppler_velocity',
    7: 'diff_reflectivity',
    9: 'copol_coeff',
    10: 'dp_phase_shift',
    16: 'corrected_mean_doppler_velocity',
    17: 'diff_phase',
    24: 'norm_coherent_power'}


# XXX these should be removed when the ctypes RSL wrapper is dropped


def ray_header_time_to_datetime(h):
    """ Return a datetime object from a RSL ray header. """
    return datetime(h.year, h.month, h.day, h.hour, h.minute, int(h.sec))


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
