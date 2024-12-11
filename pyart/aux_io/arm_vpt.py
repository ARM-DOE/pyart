"""
Routines for reading ARM vertically-pointing radar ingest (e.g., a1) files.
These files are characterized by being NetCDF files that do not fully conform
to the CF/Radial convention. Nonetheless this module borrows heavily from the
existing CF/Radial module.

"""

import netCDF4
import numpy as np

from ..config import FileMetadata
from ..core.radar import Radar
from ..io import cfradial


def read_kazr(
    filename,
    field_names=None,
    additional_metadata=None,
    file_field_names=False,
    exclude_fields=None,
    include_fields=None,
):
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
        "cfradial",
        field_names,
        additional_metadata,
        file_field_names,
        exclude_fields,
        include_fields,
    )

    # read the data
    ncobj = netCDF4.Dataset(filename)
    ncvars = ncobj.variables

    # 4.1 Global attribute -> move to metadata dictionary
    metadata = {k: getattr(ncobj, k) for k in ncobj.ncattrs()}
    metadata["n_gates_vary"] = "false"

    # 4.2 Dimensions (do nothing)

    # 4.3 Global variable -> move to metadata dictionary
    if "volume_number" in ncvars:
        metadata["volume_number"] = int(ncvars["volume_number"][:])
    else:
        metadata["volume_number"] = 0

    global_vars = {
        "platform_type": "fixed",
        "instrument_type": "radar",
        "primary_axis": "axis_z",
    }
    # ignore time_* global variables, these are calculated from the time
    # variable when the file is written.
    for var, default_value in global_vars.items():
        if var in ncvars:
            metadata[var] = str(netCDF4.chartostring(ncvars[var][:]))
        else:
            metadata[var] = default_value

    # 4.4 coordinate variables -> create attribute dictionaries
    time = cfradial._ncvar_to_dict(ncvars["time"])
    _range = cfradial._ncvar_to_dict(ncvars["range"])

    # 4.5 Ray dimension variables

    # 4.6 Location variables -> create attribute dictionaries
    # the only difference in this section to cfradial.read_cfradial is the
    # minor variable name differences:
    # latitude -> lat
    # longitude -> lon
    # altitdue -> alt
    latitude = cfradial._ncvar_to_dict(ncvars["lat"])
    longitude = cfradial._ncvar_to_dict(ncvars["lon"])
    altitude = cfradial._ncvar_to_dict(ncvars["alt"])

    # 4.7 Sweep variables -> create atrribute dictionaries
    # this is the section that needed the most work since the initial NetCDF
    # file did not contain any sweep information
    sweep_number = filemetadata("sweep_number")
    sweep_number["data"] = np.array([0], dtype=np.int32)

    sweep_mode = filemetadata("sweep_mode")
    sweep_mode["data"] = np.array(["vertical_pointing"], dtype=str)

    fixed_angle = filemetadata("fixed_angle")
    fixed_angle["data"] = np.array([90.0], dtype=np.float32)

    sweep_start_ray_index = filemetadata("sweep_start_ray_index")
    sweep_start_ray_index["data"] = np.array([0], dtype=np.int32)

    sweep_end_ray_index = filemetadata("sweep_end_ray_index")
    sweep_end_ray_index["data"] = np.array([ncvars["time"].size - 1], dtype=np.int32)

    # first sweep mode determines scan_type
    # this module is specific to vertically-pointing data
    scan_type = "vpt"

    # 4.8 Sensor pointing variables -> create attribute dictionaries
    # this section also required some changes since the initial NetCDF did not
    # contain any sensor pointing variables
    azimuth = filemetadata("azimuth")
    azimuth["data"] = 0.0 * np.ones(ncvars["time"].size, dtype=np.float32)

    elevation = filemetadata("elevation")
    elevation["data"] = 90.0 * np.ones(ncvars["time"].size, dtype=np.float32)

    # 4.9 Moving platform geo-reference variables

    # 4.10 Moments field data variables -> field attribute dictionary
    # all variables with dimensions of 'time', 'range' are fields
    keys = [k for k, v in ncvars.items() if v.dimensions == ("time", "range")]

    fields = {}
    for key in keys:
        field_name = filemetadata.get_field_name(key)
        if field_name is None:
            if exclude_fields is not None and key in exclude_fields:
                continue
            if include_fields is not None and key not in include_fields:
                continue
            field_name = key
        fields[field_name] = cfradial._ncvar_to_dict(ncvars[key])

    # 4.5 instrument_parameters sub-convention -> instrument_parameters dict
    # this section needed multiple changes and/or additions since the
    # instrument parameters were primarily located in the global attributes
    # this section is likely still incomplete
    omega = float(ncobj.radar_operating_frequency.split()[0])
    frequency = filemetadata("frequency")
    frequency["data"] = np.array([omega / 1e9], dtype=np.float32)

    prt_mode = filemetadata("prt_mode")
    prt_mode["data"] = np.array(["fixed"], dtype=str)

    prf = float(ncobj.pulse_repetition_frequency.split()[0])
    prt = filemetadata("prt")
    prt["data"] = (1.0 / prf) * np.ones(ncvars["time"].size, dtype=np.float32)

    v_nq = float(ncobj.nyquist_velocity.split()[0])
    nyquist_velocity = filemetadata("nyquist_velocity")
    nyquist_velocity["data"] = v_nq * np.ones(ncvars["time"].size, dtype=np.float32)
    samples = int(ncobj.num_spectral_averages)
    n_samples = filemetadata("n_samples")
    n_samples["data"] = samples * np.ones(ncvars["time"].size, dtype=np.int32)

    # 4.6 radar_parameters sub-convention -> instrument_parameters dict
    # this section needed multiple changes and/or additions since the
    # radar instrument parameters were primarily located in the global
    # attributes
    # this section is likely still incomplete
    instrument_parameters = {
        "frequency": frequency,
        "prt_mode": prt_mode,
        "prt": prt,
        "nyquist_velocity": nyquist_velocity,
        "n_samples": n_samples,
    }

    # 4.7 lidar_parameters sub-convention -> skip
    # 4.8 radar_calibration sub-convention -> skip

    # close NetCDF object
    ncobj.close()

    return Radar(
        time,
        _range,
        fields,
        metadata,
        scan_type,
        latitude,
        longitude,
        altitude,
        sweep_number,
        sweep_mode,
        fixed_angle,
        sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth,
        elevation,
        instrument_parameters=instrument_parameters,
    )


def read_mmcr(
    filename,
    field_names=None,
    additional_metadata=None,
    file_field_names=False,
    exclude_fields=None,
    include_fields=None,
    mode_names=None,
    mode_to_extract=None,
):
    """
    Read millimeter wavelength cloud radar (MMCR) b1 NetCDF ingest data.

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
    mode_names: dict
        Keys are mode dimension indices. Values are the corresponding mode
        abbreviations and prt. Using default values if None.
    mode_to_extract: str or None
        Abbreviation of mode to extract.
        Extracting all data for all modes if None.

    Returns
    -------
    radar_out : Radar or dict of Radar objects
        If `mode_to_extract` is None, returns a dict with keys representing operated modes and
        values the corresponding `Radar` object.
        If `mode_to_extract` is specified returning the `Radar` object corresponding to the that
        mode.

    """

    if mode_names is None:
        mode_names = {
            1: "bl",
            2: "ci",
            3: "ge",
            4: "pr",
            5: "dualpol_ch0",  # Dual-pol receiver 0
            6: "dualpol_ch1",  # Dual-pol receiver 1
        }

    # create metadata retrieval object
    filemetadata = FileMetadata(
        "cfradial",
        field_names,
        additional_metadata,
        file_field_names,
        exclude_fields,
        include_fields,
    )

    # read the data
    ncobj = netCDF4.Dataset(filename)
    ncvars = ncobj.variables

    # 4.1 Global attribute -> move to metadata dictionary
    metadata = {k: getattr(ncobj, k) for k in ncobj.ncattrs()}
    metadata["n_gates_vary"] = "false"

    # 4.2 Dimensions (do nothing)

    # 4.3 Global variable -> move to metadata dictionary
    if "volume_number" in ncvars:
        metadata["volume_number"] = int(ncvars["volume_number"][:])
    else:
        metadata["volume_number"] = 0

    global_vars = {
        "platform_type": "fixed",
        "instrument_type": "radar",
        "primary_axis": "axis_z",
    }
    # ignore time_* global variables, these are calculated from the time
    # variable when the file is written.
    for var, default_value in global_vars.items():
        if var in ncvars:
            metadata[var] = str(netCDF4.chartostring(ncvars[var][:]))
        else:
            metadata[var] = default_value

    # 4.4 coordinate variables -> create attribute dictionaries (`_range` is generated in loop below)

    # 4.5 Ray dimension variables

    # 4.6 Location variables -> create attribute dictionaries
    # the only difference in this section to cfradial.read_cfradial is the
    # minor variable name differences:
    # latitude -> lat
    # longitude -> lon
    # altitdue -> alt
    latitude = cfradial._ncvar_to_dict(ncvars["lat"])
    longitude = cfradial._ncvar_to_dict(ncvars["lon"])
    altitude = cfradial._ncvar_to_dict(ncvars["alt"])

    # 4.7 Sweep variables -> create atrribute dictionaries
    # this is the section that needed the most work since the initial NetCDF
    # file did not contain any sweep information
    sweep_number = filemetadata("sweep_number")
    sweep_number["data"] = np.array([0], dtype=np.int32)

    sweep_mode = filemetadata("sweep_mode")
    sweep_mode["data"] = np.array(["vertical_pointing"], dtype=str)

    fixed_angle = filemetadata("fixed_angle")
    fixed_angle["data"] = np.array([90.0], dtype=np.float32)

    sweep_start_ray_index = filemetadata("sweep_start_ray_index")
    sweep_start_ray_index["data"] = np.array([0], dtype=np.int32)

    # first sweep mode determines scan_type
    # this module is specific to vertically-pointing data
    scan_type = "vpt"

    # 4.8 Sensor pointing variables -> create attribute dictionaries
    # this section also required some changes since the initial NetCDF did not
    # contain any sensor pointing variables

    # 4.9 Moving platform geo-reference variables

    # 4.10 Moments field data variables -> field attribute dictionary
    # all variables with dimensions of 'time', 'heights' are fields
    keys = [k for k, v in ncvars.items() if v.dimensions == ("time", "heights")]

    # Determine if one or more modes are to be extracted and generate Radar object(s) output
    if mode_to_extract is None:
        modes = np.unique(ncvars["ModeNum"])
        radar_out = {}
    else:
        modes = [key for key in mode_names.keys() if mode_names[key] == mode_to_extract]
    for mode in modes:
        mode_sample_indices = ncvars["ModeNum"] == np.array(mode)
        active_range = ncvars["heights"][mode, :] > 0.0
        range_arr = ncvars["heights"][mode, active_range] - altitude["data"]
        time = cfradial._ncvar_to_dict(ncvars["time"])
        time["data"] = ncvars["time"][mode_sample_indices]
        _range = {
            "long_name": "Range (center of radar sample volume)",
            "units": "m AGL",
            "missing_value": -9999.0,
            "data": range_arr,
        }

        fields = {}
        for key in keys:
            field_name = filemetadata.get_field_name(key)
            if field_name is None:
                if exclude_fields is not None and key in exclude_fields:
                    continue
                if include_fields is not None and key not in include_fields:
                    continue
                field_name = key

            d = {
                k: getattr(ncvars[key], k)
                for k in ncvars[key].ncattrs()
                if k not in ["scale_factor", "add_offset"]
            }
            d["data"] = ncvars[key][mode_sample_indices, active_range]
            fields[field_name] = d

        # 4.5 instrument_parameters sub-convention -> instrument_parameters dict
        # this section needed multiple changes and/or additions since the
        # instrument parameters were primarily located in the global attributes
        # this section is likely still incomplete

        sweep_end_ray_index = filemetadata("sweep_end_ray_index")
        sweep_end_ray_index["data"] = np.array([time["data"].size - 1], dtype=np.int32)

        azimuth = filemetadata("azimuth")
        azimuth["data"] = 0.0 * np.ones(time["data"].size, dtype=np.float32)

        elevation = filemetadata("elevation")
        elevation["data"] = 90.0 * np.ones(time["data"].size, dtype=np.float32)

        omega = float(ncobj.radar_operating_frequency.split()[0])
        frequency = filemetadata("frequency")
        frequency["data"] = np.array([omega / 1e9], dtype=np.float32)

        prt_mode = filemetadata("prt_mode")
        prt_mode["data"] = np.array(["fixed"], dtype=str)

        prt = filemetadata("prt")
        prt["data"] = np.full(
            time["data"].size, ncvars["InterPulsePeriod"][mode] * 1e-9, dtype=np.float32
        )

        nyquist_velocity = filemetadata("nyquist_velocity")
        nyquist_velocity["data"] = np.full(
            time["data"].size, ncvars["NyquistVelocity"][mode], dtype=np.float32
        )

        n_samples = filemetadata("n_samples")
        n_samples["data"] = np.full(
            time["data"].size, ncvars["NumSpectralAverages"][mode], dtype=np.float32
        )

        # 4.6 radar_parameters sub-convention -> instrument_parameters dict
        # this section needed multiple changes and/or additions since the
        # radar instrument parameters were primarily located in the global
        # attributes
        # this section is likely still incomplete
        instrument_parameters = {
            "frequency": frequency,
            "prt_mode": prt_mode,
            "prt": prt,
            "nyquist_velocity": nyquist_velocity,
            "n_samples": n_samples,
        }

        # 4.7 lidar_parameters sub-convention -> skip
        # 4.8 radar_calibration sub-convention -> skip

        Radar_tmp = Radar(
            time,
            _range,
            fields,
            metadata,
            scan_type,
            latitude,
            longitude,
            altitude,
            sweep_number,
            sweep_mode,
            fixed_angle,
            sweep_start_ray_index,
            sweep_end_ray_index,
            azimuth,
            elevation,
            instrument_parameters=instrument_parameters,
        )
        if mode_to_extract is None:
            radar_out[mode_names[mode]] = Radar_tmp
        else:
            radar_out = Radar_tmp

    # close NetCDF object
    ncobj.close()

    return radar_out
