"""
pyart.io.rsl
============

Python wrapper around the RSL library.

.. autosummary::
    :toctree: generated/

    read_rsl

"""

# Nothing from this module is imported into pyart.io if RSL is not installed.
import numpy as np

from . import _rsl_interface
from .radar import Radar
from .common import dms_to_d, get_metadata, make_time_unit_str


def read_rsl(filename, radar_format=None, callid=None, add_meta=None):
    """
    Read a file supported by RSL

    Parameters
    ----------
    filename : str or RSL_radar
        Name of file whose format is supported by RSL.
    radar_format : str or None
        Format of the radar file.  Must be 'wsr88d' or None.
    callid : str or None
        Four letter NEXRAD radar Call ID, only used when radar_format is
        'wsr88d'.
    add_meta : dict or None
        Dictionary containing additional metadata to add to the created
        Radar object.  This will overwrite metadata extracted from the file.
        None is add no additional metadata.

    Returns
    -------
    radar : Radar
        Radar object.

    """
    fillvalue = -9999.0
    rslfile = _rsl_interface.RslFile(filename, radar_format, callid)
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

    # time
    time = get_metadata('time')

    t_start = first_ray.get_datetime()

    last_sweep = first_volume.get_sweep(nsweeps - 1)
    last_ray = last_sweep.get_ray(nrays - 1)
    t_end = last_ray.get_datetime()

    t_span = (t_end - t_start).seconds
    time['data'] = np.linspace(0, t_span, nrays * nsweeps)
    time['units'] = make_time_unit_str(t_start)

    # range
    _range = get_metadata('range')
    gate0 = first_ray.range_bin1
    gate_size = first_ray.gate_size
    _range['data'] = gate0 + gate_size * np.arange(ngates, dtype='float32')
    _range['meters_to_center_of_first_gate'] = _range['data'][0]
    _range['meters_between_gates'] = np.array(gate_size, dtype='float32')

    # fields
    # transfer only those which are available and have a standard name
    fields = {}
    available_vols = rslfile.available_moments()
    good_vols = VOLUMENUM2STANDARDNAME.keys()
    volumes_to_extract = [i for i in available_vols if i in good_vols]

    for volume_num in volumes_to_extract:

        # extract the field, mask and reshape
        data = rslfile.get_volume_array(volume_num)
        data[np.where(np.isnan(data))] = fillvalue
        data[np.where(data == 131072)] = fillvalue
        data = np.ma.masked_equal(data, fillvalue)
        data.shape = (data.shape[0] * data.shape[1], data.shape[2])

        # create the field dictionary
        standard_field_name = VOLUMENUM2STANDARDNAME[volume_num]
        fielddict = get_metadata(standard_field_name)
        fielddict['data'] = data
        fielddict['_FillValue'] = fillvalue
        fields[standard_field_name] = fielddict

    # metadata
    metadata = {'original_container': 'rsl'}
    rsl_dict = rslfile.get_radar_header()
    need_from_rsl_header = {
        'name': 'instrument_name', 'project': 'project', 'state': 'state',
        'country': 'country'}  # rsl_name : radar_metadata_name
    for rsl_key, metadata_key in need_from_rsl_header.iteritems():
        metadata[metadata_key] = rsl_dict[rsl_key]

    # additional required CF/Radial metadata set to blank strings
    metadata['title'] = ''
    metadata['institution'] = ''
    metadata['references'] = ''
    metadata['source'] = ''
    metadata['comment'] = ''

    if add_meta is not None:
        metadata.update(add_meta)

    # latitude
    latitude = get_metadata('latitude')
    lat = dms_to_d((rsl_dict['latd'], rsl_dict['latm'], rsl_dict['lats']))
    latitude['data'] = np.array([lat], dtype='float64')

    # longitude
    longitude = get_metadata('longitude')
    lon = dms_to_d((rsl_dict['lond'], rsl_dict['lonm'], rsl_dict['lons']))
    longitude['data'] = np.array([lon], dtype='float64')

    # altitude
    altitude = get_metadata('altitude')
    altitude['data'] = np.array([rsl_dict['height']], dtype='float64')

    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
    # sweep_end_ray_index
    sweep_number = get_metadata('sweep_number')
    sweep_mode = get_metadata('sweep_mode')
    fixed_angle = get_metadata('fixed_angle')
    sweep_start_ray_index = get_metadata('sweep_start_ray_index')
    sweep_end_ray_index = get_metadata('sweep_end_ray_index')
    len_time = len(time['data'])

    if scan_type == 'ppi':
        nsweeps = nele
        sweep_number['data'] = np.arange(nsweeps, dtype='int32')
        sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])
        fixed_angle['data'] = first_volume.get_sweep_elevs()
        sweep_start_ray_index['data'] = np.arange(0, len_time, naz,
                                                  dtype='int32')
        sweep_end_ray_index['data'] = np.arange(naz - 1, len_time, naz,
                                                dtype='int32')

    elif scan_type == 'rhi':
        nsweeps = naz
        sweep_number['data'] = np.arange(nsweeps, dtype='int32')
        sweep_mode['data'] = np.array(nsweeps * ['rhi'])
        fixed_angle['data'] = first_volume.get_sweep_azimuths()
        sweep_start_ray_index['data'] = np.arange(0, len_time, nele,
                                                  dtype='int32')
        sweep_end_ray_index['data'] = np.arange(nele - 1, len_time, nele,
                                                dtype='int32')

    # azimuth, elevation
    azimuth = get_metadata('azimuth')
    elevation = get_metadata('elevation')
    _azimuth, _elevation = first_volume.get_azimuth_and_elev_array()
    azimuth['data'] = _azimuth.flatten()
    elevation['data'] = _elevation.flatten()

    # instrument_parameters
    prt = get_metadata('prt')
    prt_mode = get_metadata('prt_mode')
    nyquist_velocity = get_metadata('nyquist_velocity')
    unambiguous_range = get_metadata('unambiguous_range')

    pm_data, nv_data, pr_data, ur_data = first_volume.get_instr_params()
    prt['data'] = pr_data.flatten()
    prt_mode['data'] = pm_data
    nyquist_velocity['data'] = nv_data.flatten()
    unambiguous_range['data'] = ur_data.flatten()

    instrument_parameters = {'unambiguous_range': unambiguous_range,
                             'prt_mode': prt_mode, 'prt': prt,
                             'nyquist_velocity': nyquist_velocity}

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)


#####################
# private functions #
#####################

#   id  abbreviation    common_name     standard_name
#   0   DZ              DBZ_F           reflectivity_horizontal
#   1   VR              VEL_F           mean_doppler_velocity
#   2   SW              WIDTH
#   3   CZ              DBZ
#   4   ZT              DBZ             reflectivity_horizontal_filtered
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


VOLUMENUM2RSLNAME = {
    0: 'DZ',
    1: 'VR',
    2: 'SW',
    3: 'CZ',
    4: 'ZT',
    5: 'DR',
    6: 'LR',
    7: 'ZD',
    8: 'DM',
    9: 'RH',
    10: 'PH',
    11: 'XZ',
    12: 'CD',
    13: 'MZ',
    14: 'MD',
    15: 'ZE',
    16: 'VE',
    17: 'KD',
    18: 'TI',
    19: 'DX',
    20: 'CH',
    21: 'AH',
    22: 'CV',
    23: 'AV',
    24: 'SQ'}


RSLNAME2VOLUMENUM = dict([(v, k) for k, v in VOLUMENUM2RSLNAME.iteritems()])


VOLUMENUM2STANDARDNAME = {
    0: 'reflectivity_horizontal_filtered',
    1: 'mean_doppler_velocity',
    4: 'reflectivity_horizontal',
    7: 'diff_reflectivity',
    9: 'copol_coeff',
    10: 'dp_phase_shift',
    16: 'corrected_mean_doppler_velocity',
    17: 'diff_phase',
    24: 'norm_coherent_power'}
