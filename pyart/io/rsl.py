"""
pyart.io.rsl
============

Python wrapper around the RSL library.

.. autosummary::
    :toctree: generated/

    read_rsl
    VOLUMENUM2RSLNAME
    RSLNAME2VOLUMENUM

"""

# Nothing from this module is imported into pyart.io if RSL is not installed.
import numpy as np

from ..config import _FileMetadata, get_fillvalue
from . import _rsl_interface
from .radar import Radar
from .common import dms_to_d, make_time_unit_str


def read_rsl(filename, field_names=None, additional_metadata=None,
             file_field_names=False, exclude_fields=None,
             radar_format=None, callid=None):
    """
    Read a file supported by RSL

    Parameters
    ----------
    filename : str or RSL_radar
        Name of file whose format is supported by RSL.
    field_names : dict, optional
        Dictionary mapping RSL data type names to radar field names. If a
        data type found in the file does not appear in this dictionary or has
        a value of None it will not be placed in the radar.fields dictionary.
        A value of None, the default, will use the mapping defined in the
        Py-ART configuration file.
    additional_metadata : dict of dicts, optional
        Dictionary of dictionaries to retrieve metadata from during this read.
        This metadata is not used during any successive file reads unless
        explicitly included.  A value of None, the default, will not
        introduct any addition metadata and the file specific or default
        metadata as specified by the Py-ART configuration file will be used.
    file_field_names : bool, optional
        True to use the RSL data type names for the field names. If this
        case the field_names parameter is ignored. The field dictionary will
        likely only have a 'data' key, unless the fields are defined in
        `additional_metadata`.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters.
    radar_format : str or None
        Format of the radar file.  Must be 'wsr88d' or None.
    callid : str or None
        Four letter NEXRAD radar Call ID, only used when radar_format is
        'wsr88d'.

    Returns
    -------
    radar : Radar
        Radar object.

    """
    # create metadata retrieval object
    filemetadata = _FileMetadata('rsl', field_names, additional_metadata,
                                 file_field_names, exclude_fields)

    # read the file
    fillvalue = get_fillvalue()
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
    time = filemetadata('time')

    t_start = first_ray.get_datetime()

    last_sweep = first_volume.get_sweep(nsweeps - 1)
    last_ray = last_sweep.get_ray(nrays - 1)
    t_end = last_ray.get_datetime()

    t_span = (t_end - t_start).seconds
    time['data'] = np.linspace(0, t_span, nrays * nsweeps)
    time['units'] = make_time_unit_str(t_start)

    # range
    _range = filemetadata('range')
    gate0 = first_ray.range_bin1
    gate_size = first_ray.gate_size
    _range['data'] = gate0 + gate_size * np.arange(ngates, dtype='float32')
    _range['meters_to_center_of_first_gate'] = _range['data'][0]
    _range['meters_between_gates'] = np.array(gate_size, dtype='float32')

    # fields
    # transfer only those which are available and have a standard name
    fields = {}
    for volume_num in available_vols:

        rsl_field_name = VOLUMENUM2RSLNAME[volume_num]
        field_name = filemetadata.get_field_name(rsl_field_name)

        # extract the field, mask and reshape
        data = rslfile.get_volume_array(volume_num)
        data[np.where(np.isnan(data))] = fillvalue
        data[np.where(data == 131072)] = fillvalue
        data = np.ma.masked_equal(data, fillvalue)
        data.shape = (data.shape[0] * data.shape[1], data.shape[2])

        # create the field dictionary
        field_dic = filemetadata(field_name)
        field_dic['data'] = data
        field_dic['_FillValue'] = fillvalue
        fields[field_name] = field_dic

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

    # latitude
    latitude = filemetadata('latitude')
    lat = dms_to_d((rsl_dict['latd'], rsl_dict['latm'], rsl_dict['lats']))
    latitude['data'] = np.array([lat], dtype='float64')

    # longitude
    longitude = filemetadata('longitude')
    lon = dms_to_d((rsl_dict['lond'], rsl_dict['lonm'], rsl_dict['lons']))
    longitude['data'] = np.array([lon], dtype='float64')

    # altitude
    altitude = filemetadata('altitude')
    altitude['data'] = np.array([rsl_dict['height']], dtype='float64')

    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
    # sweep_end_ray_index
    sweep_number = filemetadata('sweep_number')
    sweep_mode = filemetadata('sweep_mode')
    fixed_angle = filemetadata('fixed_angle')
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')
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
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')
    _azimuth, _elevation = first_volume.get_azimuth_and_elev_array()
    azimuth['data'] = _azimuth.flatten()
    elevation['data'] = _elevation.flatten()

    # instrument_parameters
    prt = filemetadata('prt')
    prt_mode = filemetadata('prt_mode')
    nyquist_velocity = filemetadata('nyquist_velocity')
    unambiguous_range = filemetadata('unambiguous_range')

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
    24: 'SQ',
    25: 'VS',
    26: 'VL',
    27: 'VG',
    28: 'VT',
    29: 'NP',
    30: 'HC',
    31: 'VC',
    32: 'V2',
    33: 'S2',
    34: 'V3',
    35: 'S3',
}

RSLNAME2VOLUMENUM = dict([(v, k) for k, v in VOLUMENUM2RSLNAME.iteritems()])
