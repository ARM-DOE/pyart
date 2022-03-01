"""
Routines for reading ARM vertically-pointing radar ingest (e.g., a1) files.
These files are characterized by being NetCDF files that do not fully conform
to the CF/Radial convention. Nonetheless this module borrows heavily from the
existing CF/Radial module.

"""

import netCDF4
import numpy as np

from ..io import cfradial
from ..config import FileMetadata
from ..core.radar import Radar


def read_kazr(filename, field_names=None, additional_metadata=None,
              file_field_names=False, exclude_fields=None,
              include_fields=None):
    """
    Read K-band ARM Zenith Radar (KAZR) NetCDF ingest data.

    Parameters
    ----------
    filename : str
        Name of NetCDF file to read data from.
    field_names : dict, optional
        Dictionary mapping field names in the file names to radar field names.
        Unlike other read functions, fields not in this dictionary or having a
        value of None are still included in the radar.fields dictionary, to
        exclude them use the `exclude_fields` parameter. Fields which are
        mapped by this dictionary will be renamed from key to value.
    additional_metadata : dict of dicts, optional
        This parameter is not used, it is included for uniformity.
    file_field_names : bool, optional
        True to force the use of the field names from the file in which
        case the `field_names` parameter is ignored. False will use to
        `field_names` parameter to rename fields.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters. Set
        to None to include all fields specified by include_fields.
    include_fields : list or None, optional
        List of fields to include from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters. Set
        to None to include all fields not specified by exclude_fields.

    Returns
    -------
    radar : Radar
        Radar object.
    """

    # create metadata retrieval object
    filemetadata = FileMetadata(
        'cfradial', field_names, additional_metadata, file_field_names,
        exclude_fields, include_fields)

    # read the data
    ncobj = netCDF4.Dataset(filename)
    ncvars = ncobj.variables

    # 4.1 Global attribute -> move to metadata dictionary
    metadata = dict([(k, getattr(ncobj, k)) for k in ncobj.ncattrs()])
    metadata['n_gates_vary'] = 'false'

    # 4.2 Dimensions (do nothing)

    # 4.3 Global variable -> move to metadata dictionary
    if 'volume_number' in ncvars:
        metadata['volume_number'] = int(ncvars['volume_number'][:])
    else:
        metadata['volume_number'] = 0

    global_vars = {'platform_type': 'fixed', 'instrument_type': 'radar',
                   'primary_axis': 'axis_z'}
    # ignore time_* global variables, these are calculated from the time
    # variable when the file is written.
    for var, default_value in global_vars.items():
        if var in ncvars:
            metadata[var] = str(netCDF4.chartostring(ncvars[var][:]))
        else:
            metadata[var] = default_value

    # 4.4 coordinate variables -> create attribute dictionaries
    time = cfradial._ncvar_to_dict(ncvars['time'])
    _range = cfradial._ncvar_to_dict(ncvars['range'])

    # 4.5 Ray dimension variables

    # 4.6 Location variables -> create attribute dictionaries
    # the only difference in this section to cfradial.read_cfradial is the
    # minor variable name differences:
    # latitude -> lat
    # longitude -> lon
    # altitdue -> alt
    latitude = cfradial._ncvar_to_dict(ncvars['lat'])
    longitude = cfradial._ncvar_to_dict(ncvars['lon'])
    altitude = cfradial._ncvar_to_dict(ncvars['alt'])

    # 4.7 Sweep variables -> create atrribute dictionaries
    # this is the section that needed the most work since the initial NetCDF
    # file did not contain any sweep information
    sweep_number = filemetadata('sweep_number')
    sweep_number['data'] = np.array([0], dtype=np.int32)

    sweep_mode = filemetadata('sweep_mode')
    sweep_mode['data'] = np.array(['vertical_pointing'], dtype=np.str)

    fixed_angle = filemetadata('fixed_angle')
    fixed_angle['data'] = np.array([90.0], dtype=np.float32)

    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_start_ray_index['data'] = np.array([0], dtype=np.int32)

    sweep_end_ray_index = filemetadata('sweep_end_ray_index')
    sweep_end_ray_index['data'] = np.array(
        [ncvars['time'].size - 1], dtype=np.int32)

    # first sweep mode determines scan_type
    # this module is specific to vertically-pointing data
    scan_type = 'vpt'

    # 4.8 Sensor pointing variables -> create attribute dictionaries
    # this section also required some changes since the initial NetCDF did not
    # contain any sensor pointing variables
    azimuth = filemetadata('azimuth')
    azimuth['data'] = 0.0 * np.ones(ncvars['time'].size, dtype=np.float32)

    elevation = filemetadata('elevation')
    elevation['data'] = 90.0 * np.ones(ncvars['time'].size, dtype=np.float32)

    # 4.9 Moving platform geo-reference variables

    # 4.10 Moments field data variables -> field attribute dictionary
    # all variables with dimensions of 'time', 'range' are fields
    keys = [k for k, v in ncvars.items() if
            v.dimensions == ('time', 'range')]

    fields = {}
    for key in keys:
        field_name = filemetadata.get_field_name(key)
        if field_name is None:
            if exclude_fields is not None and key in exclude_fields:
                continue
            if include_fields is not None and not key in include_fields:
                continue
            field_name = key
        fields[field_name] = cfradial._ncvar_to_dict(ncvars[key])

    # 4.5 instrument_parameters sub-convention -> instrument_parameters dict
    # this section needed multiple changes and/or additions since the
    # instrument parameters were primarily located in the global attributes
    # this section is likely still incomplete
    omega = float(ncobj.radar_operating_frequency.split()[0])
    frequency = filemetadata('frequency')
    frequency['data'] = np.array([omega / 1e9], dtype=np.float32)

    prt_mode = filemetadata('prt_mode')
    prt_mode['data'] = np.array(['fixed'], dtype=np.str)

    prf = float(ncobj.pulse_repetition_frequency.split()[0])
    prt = filemetadata('prt')
    prt['data'] = (1.0 / prf) * np.ones(ncvars['time'].size, dtype=np.float32)

    v_nq = float(ncobj.nyquist_velocity.split()[0])
    nyquist_velocity = filemetadata('nyquist_velocity')
    nyquist_velocity['data'] = v_nq * np.ones(ncvars['time'].size,
                                              dtype=np.float32),
    samples = int(ncobj.num_spectral_averages)
    n_samples = filemetadata('n_samples')
    n_samples['data'] = samples * np.ones(ncvars['time'].size, dtype=np.int32)

    # 4.6 radar_parameters sub-convention -> instrument_parameters dict
    # this section needed multiple changes and/or additions since the
    # radar instrument parameters were primarily located in the global
    # attributes
    # this section is likely still incomplete
    instrument_parameters = {
        'frequency': frequency,
        'prt_mode': prt_mode,
        'prt': prt,
        'nyquist_velocity': nyquist_velocity,
        'n_samples': n_samples,
    }

    # 4.7 lidar_parameters sub-convention -> skip
    # 4.8 radar_calibration sub-convention -> skip

    # close NetCDF object
    ncobj.close()

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)
