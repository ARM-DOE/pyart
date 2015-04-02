"""
pyart.aux_io.arm_vpt
====================

Routines for reading ARM vertically-pointing radar ingest (e.g., a1) files.
These files are characterized by being NetCDF files that do not fully conform
to the CF/Radial convention. Nonetheless this module borrows heavily from the
existing CF/Radial module.

"""

import re
import netCDF4
import numpy as np

from ..io import cfradial
from ..config import FileMetadata
from ..core.radar import Radar


def read_kazr(filename, field_names=None, additional_metadata=None,
              file_field_names=False, exclude_fields=None):
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
        after the `file_field_names` and `field_names` parameters.

    Returns
    -------
    radar : Radar
        Radar object.
    """

    # create metadata retrieval object
    filemetadata = FileMetadata(
        'cfradial', field_names, additional_metadata, file_field_names,
        exclude_fields)

    # read the data
    ncobj = netCDF4.Dataset(filename)
    ncvars = ncobj.variables

    # 4.1 Global attribute -> move to metadata dictionary
    metadata = dict([(k, getattr(ncobj, k)) for k in ncobj.ncattrs()])
    if 'n_gates_vary' in metadata:
        metadata['n_gates_vary'] = 'false'  # corrected below

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
    for var, default_value in global_vars.iteritems():
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
    if 'altitude_agl' in ncvars:
        altitude_agl = cfradial._ncvar_to_dict(ncvars['altitude_agl'])
    else:
        altitude_agl = None

    # 4.7 Sweep variables -> create atrribute dictionaries
    # this is the section that needed the most work since the initial NetCDF
    # file did not contain any sweep information
    sweep_number = {
        'data': np.array([0], dtype=np.int32),
        'long_name': 'Sweep number',
        'standard_name': 'sweep_number',
        'units': 'count',
    }
    sweep_mode = {
        'data': np.array(['vertical_pointing'], dtype=np.str),
        'long_name': 'Sweep mode',
        'standard_name': 'sweep_mode',
        'units': 'unitless',
        'comment': ('Options are: "sector", "coplane", "rhi" '
                    '"vertical_pointing", "idle", "azimuth_surveillance" '
                    '"elevation_surveillance", "sunscan", "pointing" '
                    '"manual_ppi", "manual_rhi"')
    }
    fixed_angle = {
        'data': np.array([90.0], dtype=np.float32),
        'long_name': 'Target angle for sweep',
        'standard_name': 'target_fixed_angle',
        'units': 'degrees',
    }
    sweep_start_ray_index = {
        'data': np.array([0], dtype=np.int32),
        'long_name': 'Index of first ray in sweep, 0-based',
        'standard_name': 'index_of_first_ray_in_sweep',
        'units': 'count',
    }
    sweep_end_ray_index = {
        'data': np.array([ncvars['time'].size - 1], dtype=np.int32),
        'long_name': 'Index of last ray in sweep, 0-based',
        'standard_name': 'index_of_last_ray_in_sweep',
        'units': 'count',
    }

    #sweep_number = _ncvar_to_dict(ncvars['sweep_number'])
    #sweep_mode = _ncvar_to_dict(ncvars['sweep_mode'])
    #fixed_angle = _ncvar_to_dict(ncvars['fixed_angle'])
    #sweep_start_ray_index = _ncvar_to_dict(ncvars['sweep_start_ray_index'])
    #sweep_end_ray_index = _ncvar_to_dict(ncvars['sweep_end_ray_index'])

    if 'target_scan_rate' in ncvars:
        target_scan_rate = cfradial._ncvar_to_dict(ncvars['target_scan_rate'])
    else:
        target_scan_rate = None
    if 'rays_are_indexed' in ncvars:
        rays_are_indexed = cfradial._ncvar_to_dict(ncvars['rays_are_indexed'])
    else:
        rays_are_indexed = None
    if 'ray_angle_res' in ncvars:
        ray_angle_res = cfradial._ncvar_to_dict(ncvars['ray_angle_res'])
    else:
        ray_angle_res = None

    # first sweep mode determines scan_type
    # this module is specific to vertically-pointing data
    mode = 'vertical_pointing'
    scan_type = 'vpt'

    # 4.8 Sensor pointing variables -> create attribute dictionaries
    # this section also required some changes since the initial NetCDF did not
    # contain any sensor pointing variables
    azimuth = {
        'data': 0.0 * np.ones(ncvars['time'].size, dtype=np.float32),
        'long_name': 'Azimuth angle relative to true north',
        'standard_name': 'beam_azimuth_angle',
        'units': 'degrees',
        'axis': 'radial_azimuth_angle',
        'comment': 'Azimuth of antenna relative to true north',
    }
    elevation = {
        'data': 90.0 * np.ones(ncvars['time'].size, dtype=np.float32),
        'long_name': 'Elevation angle from horizontal plane',
        'standard_name': 'beam_elevation_angle',
        'units': 'degrees',
        'axis': 'radial_elevation_coordinate',
        'comment': 'Elevation of antenna relative to the horizontal plane',
    }

    #azimuth = _ncvar_to_dict(ncvars['azimuth'])
    #elevation = _ncvar_to_dict(ncvars['elevation'])

    if 'scan_rate' in ncvars:
        scan_rate = cfradial._ncvar_to_dict(ncvars['scan_rate'])
    else:
        scan_rate = None

    if 'antenna_transition' in ncvars:
        antenna_transition = cfradial._ncvar_to_dict(
            ncvars['antenna_transition'])
    else:
        antenna_transition = None

    # 4.9 Moving platform geo-reference variables
    # Aircraft specific varaibles
    if 'rotation' in ncvars:
        rotation = _ncvar_to_dict(ncvars['rotation'])
    else:
        rotation = None
    if 'tilt' in ncvars:
        tilt = _ncvar_to_dict(ncvars['tilt'])
    else:
        tilt = None
    if 'roll' in ncvars:
        roll = _ncvar_to_dict(ncvars['roll'])
    else:
        roll = None
    if 'drift' in ncvars:
        drift = _ncvar_to_dict(ncvars['drift'])
    else:
        drift = None
    if 'heading' in ncvars:
        heading = _ncvar_to_dict(ncvars['heading'])
    else:
        heading = None
    if 'pitch' in ncvars:
        pitch = _ncvar_to_dict(ncvars['pitch'])
    else:
        pitch = None
    if 'georefs_applied' in ncvars:
        georefs_applied = _ncvar_to_dict(ncvars['georefs_applied'])
    else:
        georefs_applied = None

    # 4.10 Moments field data variables -> field attribute dictionary
    if 'ray_n_gates' in ncvars:
        # all variables with dimensions of n_points are fields.
        keys = [k for k, v in ncvars.iteritems()
                if v.dimensions == ('n_points', )]
    else:
        # all variables with dimensions of 'time', 'range' are fields
        keys = [k for k, v in ncvars.iteritems()
                if v.dimensions == ('time', 'range')]

    fields = {}
    for key in keys:
        field_name = filemetadata.get_field_name(key)
        if field_name is None:
            if exclude_fields is not None and key in exclude_fields:
                continue
            field_name = key
        fields[field_name] = cfradial._ncvar_to_dict(ncvars[key])

    if 'ray_n_gates' in ncvars:
        shape = (len(ncvars['time']), len(ncvars['range']))
        ray_n_gates = ncvars['ray_n_gates'][:]
        ray_start_index = ncvars['ray_start_index'][:]
        for dic in fields.values():
            cfradial._unpack_variable_gate_field_dic(
                dic, shape, ray_n_gates, ray_start_index)

    # 4.5 instrument_parameters sub-convention -> instrument_parameters dict
    # this section needed multiple changes and/or additions since the
    # instrument parameters were primarily located in the global attributes
    # this section is likely still incomplete
    omega = float(ncobj.radar_operating_frequency.split()[0])
    frequency = {
        'data': np.array([omega / 1e9], dtype=np.float32),
        'long_name': 'Radar operating frequency',
        'standard_name': 'radiation_frequency',
        'units': 'Hz',
        'meta_group': 'instrument_parameters',
    }
    prt_mode = {
        'data': np.array(['fixed'], dtype=np.str),
        'long_name': 'Pulsing mode',
        'standard_name': 'transmit_pulse_mode',
        'units': 'unitless',
        'meta_group': 'instrument_parameters',
        'comment': ('Options are "fixed", "staggered", "dual". Assumed '
                    '"fixed" if missing'),
    }
    prf = float(ncobj.pulse_repetition_frequency.split()[0])
    prt = {
        'data': (1.0 / prf) * np.ones(ncvars['time'].size, dtype=np.float32),
        'long_name': 'Pulse repetition time',
        'standard_name': 'pulse_repetition_time',
        'units': 'seconds',
        'meta_group': 'instrument_parameters',
        'comments': 'For staggered PRT, also see prt_ratio',
    }
    v_nq = float(ncobj.nyquist_velocity.split()[0])
    nyquist_velocity = {
        'data': v_nq * np.ones(ncvars['time'].size, dtype=np.float32),
        'long_name': 'Nyquist velocity',
        'standard_name': 'unambiguous_doppler_velocity',
        'units': 'meters_per_second',
        'meta_group': 'instrument_parameters',
    }
    samples = int(ncobj.num_spectral_averages)
    n_samples = {
        'data': samples * np.ones(ncvars['time'].size, dtype=np.int32),
        'long_name': 'Number of samples used to compute moments',
        'standard_name': 'number_of_samples_used_to_compute_moments',
        'units': 'unitless',
        'meta_group': 'instrument_parameters',
    }

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

    #keys = [k for k in cfradial._INSTRUMENT_PARAMS_DIMS.keys() if k in ncvars]
    #instrument_parameters = dict(
    #    (k, cfradial._ncvar_to_dict(ncvars[k])) for k in keys)
    #if instrument_parameters == {}:  # if no parameters set to None
    #    instrument_parameters = None

    # 4.7 lidar_parameters sub-convention -> skip

    # 4.8 radar_calibration sub-convention -> radar_calibration
    keys = cfradial._find_all_meta_group_vars(ncvars, 'radar_calibration')
    radar_calibration = dict(
        (k, cfradial._ncvar_to_dict(ncvars[k])) for k in keys)
    if radar_calibration == {}:
        radar_calibration = None

    # close NetCDF object
    ncobj.close()

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters,
        radar_calibration=radar_calibration,
        altitude_agl=altitude_agl,
        scan_rate=scan_rate,
        antenna_transition=antenna_transition,
        target_scan_rate=target_scan_rate,
        rays_are_indexed=rays_are_indexed, ray_angle_res=ray_angle_res,
        rotation=rotation, tilt=tilt, roll=roll, drift=drift, heading=heading,
        pitch=pitch, georefs_applied=georefs_applied)
