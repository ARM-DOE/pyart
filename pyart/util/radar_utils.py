"""
Functions for working radar instances.

"""

import copy

import numpy as np
import numpy.ma as ma
from netCDF4 import date2num

from ..config import get_fillvalue
from . import datetime_utils


def is_vpt(radar, offset=0.5):
    """
    Determine if a Radar appears to be a vertical pointing scan.

    This function only verifies that the object is a vertical pointing scan,
    use the :py:func:`to_vpt` function to convert the radar to a vpt scan
    if this function returns True.

    Parameters
    ----------
    radar : Radar
        Radar object to determine if.
    offset : float, optional
        Maximum offset of the elevation from 90 degrees to still consider
        to be vertically pointing.

    Returns
    -------
    flag : bool
        True if the radar appear to be verticle pointing, False if not.

    """
    # check that the elevation is within offset of 90 degrees.
    elev = radar.elevation["data"]
    return np.all((elev < 90.0 + offset) & (elev > 90.0 - offset))


def to_vpt(radar, single_scan=True):
    """
    Convert an existing Radar object to represent a vertical pointing scan.

    This function does not verify that the Radar object contains a vertical
    pointing scan. To perform such a check use :py:func:`is_vpt`.

    Parameters
    ----------
    radar : Radar
        Mislabeled vertical pointing scan Radar object to convert to be
        properly labeled. This object is converted in place, no copy of
        the existing data is made.
    single_scan : bool, optional
        True to convert the volume to a single scan, any azimuth angle data
        is lost. False will convert the scan to contain the same number of
        scans as rays, azimuth angles are retained.

    """
    if single_scan:
        nsweeps = 1
        radar.azimuth["data"][:] = 0.0
        seri = np.array([radar.nrays - 1], dtype="int32")
        radar.sweep_end_ray_index["data"] = seri
    else:
        nsweeps = radar.nrays
        # radar.azimuth not adjusted
        radar.sweep_end_ray_index["data"] = np.arange(nsweeps, dtype="int32")

    radar.scan_type = "vpt"
    radar.nsweeps = nsweeps
    radar.target_scan_rate = None  # no scanning
    radar.elevation["data"][:] = 90.0

    radar.sweep_number["data"] = np.arange(nsweeps, dtype="int32")
    radar.sweep_mode["data"] = np.array(["vertical_pointing"] * nsweeps)
    radar.fixed_angle["data"] = np.ones(nsweeps, dtype="float32") * 90.0
    radar.sweep_start_ray_index["data"] = np.arange(nsweeps, dtype="int32")

    if radar.instrument_parameters is not None:
        for key in ["prt_mode", "follow_mode", "polarization_mode"]:
            if key in radar.instrument_parameters:
                ip_dic = radar.instrument_parameters[key]
                ip_dic["data"] = np.array([ip_dic["data"][0]] * nsweeps)

    # Attributes that do not need any changes
    # radar.altitude
    # radar.altitude_agl
    # radar.latitude
    # radar.longitude

    # radar.range
    # radar.ngates
    # radar.nrays

    # radar.metadata
    # radar.radar_calibration

    # radar.time
    # radar.fields
    # radar.antenna_transition
    # radar.scan_rate
    return


def determine_sweeps(
    radar,
    max_offset=0.1,
    running_win_dt=5.0,
    deg_rng=(-5.0, 360.0),
    consider_2pi_jump=True,
):
    """
    Determine the number of sweeps using elevation data (PPI scans) or azimuth
    data (RHI scans) and update the input radar object

    Parameters
    ----------
    radar : Radar object
        The radar object containing the data.
    max_offset : float
        Maximum elevation offset (if is_ppi is True) or azimuth offset (if
        is_ppi is False) allowed to determine sweeps.
    running_win_dt: float
        running window period (in seconds) used to determine elevation or
        azimuth shifts.
        Note: set wisely: the method assumes that a single sweep is longer than this
        parameter.
    deg_rng: float
        angle range (azimuth or elevation) to consider for calculations.
        Assuming azimuths between 0 to 360, this should be equal to (0., 360.), but
        given that there could be ppi scan strategies at negative elevations,
        one might consider a negative values (current default), or , for example,
        -180 to 180 if the azimuth range goes from -180 to 180.
    consider_2pi_jump: bool
        if True and radar scan type is 'rhi', overwriting deg_rng to (0., 360.), and
        merging the first and last azimuth bins (to have shots just below 360 and
        just above 0 to be considered part of the same sweep).

    """
    # set fixed and variable coordinates depending on scan type
    # ======================
    if "rhi" in radar.scan_type.lower():
        var_array = radar.elevation["data"]
        fix_array = radar.azimuth["data"]
        if consider_2pi_jump:
            deg_rng = (0.0, 360.0)
    else:  # ppi or vpt
        var_array = radar.azimuth["data"]
        fix_array = radar.elevation["data"]

    # set bins and parameters and allocate lists
    # ======================
    angle_bins = np.arange(
        deg_rng[0] - max_offset, deg_rng[1] + max_offset + 1e-10, max_offset * 2.0
    )
    sample_dt = np.nanmean(np.diff(radar.time["data"]))
    win_size = int(np.ceil(running_win_dt / sample_dt))
    if win_size < 2:
        raise ValueError(
            "Window size <= 1; consider decreasing the value of running_win_dt"
        )
    sweep_start_index, sweep_end_index = [], []
    in_sweep = False  # determine if sweep is underway in current index

    # Loop through coordinate data and detect sweep edges
    # ======================
    t = 0
    while t < radar.time["data"].size - win_size + 1:
        var_win = var_array[t : t + win_size]
        fix_win = fix_array[t : t + win_size]
        idle_sweep = np.diff(var_win) == 0
        if idle_sweep[0]:  # sweep did not start
            t += 1
            continue
        bincounts, _ = np.histogram(fix_win, bins=angle_bins)
        if ("rhi" in radar.scan_type.lower()) & consider_2pi_jump:
            bincounts[0] += bincounts[-1]
            bincounts = bincounts[:-1]
        moving_radar = np.sum(bincounts > 0) > 1  # radar transition to a new sweep
        if in_sweep:
            if t == radar.time["data"].size - win_size:
                sweep_end_index.append(radar.time["data"].size - 1)
            elif moving_radar:
                in_sweep = False
                sweep_end_index.append(t + win_size - 2)
                t += win_size - 2
        elif np.all(~idle_sweep) & ~moving_radar:
            in_sweep = True
            sweep_start_index.append(t)
        t += 1
    sweep_number = np.arange(len(sweep_start_index))

    # Update radar object
    # ======================
    radar.sweep_start_ray_index["data"] = ma.array(sweep_start_index, dtype="int32")
    radar.sweep_end_ray_index["data"] = ma.array(sweep_end_index, dtype="int32")
    radar.sweep_number["data"] = ma.array(sweep_number, dtype="int32")
    fixed_angle = [
        np.mean(fix_array[si : ei + 1])
        for si, ei in zip(
            radar.sweep_start_ray_index["data"], radar.sweep_end_ray_index["data"]
        )
    ]
    radar.fixed_angle["data"] = ma.array(fixed_angle, dtype="float32")
    radar.nsweeps = len(sweep_number)
    transition = np.zeros(radar.nrays)
    for i in range(radar.nsweeps):
        if i == 0:
            transition[: sweep_start_index[i]] = 1
        else:
            transition[sweep_end_index[i - 1] : sweep_start_index[i]] = 1
    radar.antenna_transition["data"] = ma.array(transition, dtype="int32")
    bstr_entry = np.array([x for x in f"{radar.scan_type:<22}"], dtype="|S1")
    radar.sweep_mode = ma.array(np.tile(bstr_entry[np.newaxis, :], (radar.nsweeps, 1)))
    return


def subset_radar(
    radar,
    field_names,
    rng_min=None,
    rng_max=None,
    ele_min=None,
    ele_max=None,
    azi_min=None,
    azi_max=None,
):
    """
    Creates a subset of the radar object along new dimensions.

    Parameters
    ----------
    radar : radar object
        The radar object containing the data.
    field_names : str or None
        The fields to keep in the new radar.
    rng_min, rng_max : float
        The range limits [m]. If None the entire coverage of the radar is
        going to be used.
    ele_min, ele_max, azi_min, azi_max : float or None
        The limits of the grid [deg]. If None the limits will be the limits
        of the radar volume.

    Returns
    -------
    radar : radar object
        The radar object containing only the desired data.

    """
    radar_aux = copy.deepcopy(radar)

    if (
        rng_min is None
        and rng_max is None
        and ele_min is None
        and ele_max is None
        and azi_min is None
        and azi_max is None
    ):
        return radar_aux

    if rng_min is None:
        rng_min = 0.0
    if rng_max is None:
        rng_max = np.max(radar_aux.range["data"])

    ind_rng = np.where(
        np.logical_and(
            radar_aux.range["data"] >= rng_min, radar_aux.range["data"] <= rng_max
        )
    )[0]

    if ind_rng.size == 0:
        print("No range bins between " + str(rng_min) + " and " + str(rng_max) + " m")
        return None

    # Determine angle limits
    if radar_aux.scan_type == "ppi":
        if ele_min is None:
            ele_min = np.min(radar_aux.fixed_angle["data"])
        if ele_max is None:
            ele_max = np.max(radar_aux.fixed_angle["data"])
        if azi_min is None:
            azi_min = np.min(radar_aux.azimuth["data"])
        if azi_max is None:
            azi_max = np.max(radar_aux.azimuth["data"])
    else:
        if ele_min is None:
            ele_min = np.min(radar_aux.elevation["data"])
        if ele_max is None:
            ele_max = np.max(radar_aux.elevation["data"])
        if azi_min is None:
            azi_min = np.min(radar_aux.fixed_angle["data"])
        if azi_max is None:
            azi_max = np.max(radar_aux.fixed_angle["data"])

    if radar_aux.scan_type == "ppi":
        # Get radar elevation angles within limits
        ele_vec = np.sort(radar_aux.fixed_angle["data"])
        ele_vec = ele_vec[np.logical_and(ele_vec >= ele_min, ele_vec <= ele_max)]
        if ele_vec.size == 0:
            print(
                "No elevation angles between " + str(ele_min) + " and " + str(ele_max)
            )
            return None

        # get sweeps corresponding to the desired elevation angles
        ind_sweeps = []
        for ele in ele_vec:
            ind_sweeps.append(np.where(radar_aux.fixed_angle["data"] == ele)[0][0])
        radar_aux = radar_aux.extract_sweeps(ind_sweeps)

        # Get indices of rays within limits
        if azi_min < azi_max:
            ind_rays = np.where(
                np.logical_and(
                    radar_aux.azimuth["data"] >= azi_min,
                    radar_aux.azimuth["data"] <= azi_max,
                )
            )[0]
        else:
            ind_rays = np.where(
                np.logical_or(
                    radar_aux.azimuth["data"] >= azi_min,
                    radar_aux.azimuth["data"] <= azi_max,
                )
            )[0]

    else:
        # Get radar azimuth angles within limits
        azi_vec = radar_aux.fixed_angle["data"]
        if azi_min < azi_max:
            azi_vec = np.sort(
                azi_vec[np.logical_and(azi_vec >= azi_min, azi_vec <= azi_max)]
            )
        else:
            azi_vec = azi_vec[np.logical_or(azi_vec >= azi_min, azi_vec <= azi_max)]
        if azi_vec.size == 0:
            print("No azimuth angles between " + str(azi_min) + " and " + str(azi_max))
            return None

        # get sweeps corresponding to the desired azimuth angles
        ind_sweeps = []
        for azi in azi_vec:
            ind_sweeps.append(np.where(radar_aux.fixed_angle["data"] == azi)[0][0])
        radar_aux = radar_aux.extract_sweeps(ind_sweeps)

        # Get indices of rays within limits
        ind_rays = np.where(
            np.logical_and(
                radar_aux.elevation["data"] >= ele_min,
                radar_aux.elevation["data"] <= ele_max,
            )
        )[0]

    # get new sweep start index and stop index
    sweep_start_inds = copy.deepcopy(radar_aux.sweep_start_ray_index["data"])
    sweep_end_inds = copy.deepcopy(radar_aux.sweep_end_ray_index["data"])

    nrays = 0
    ind_rays_aux = []
    for j in range(radar_aux.nsweeps):
        # get ray indices for this sweep
        ind_rays_sweep = ind_rays[
            np.logical_and(
                ind_rays >= sweep_start_inds[j], ind_rays <= sweep_end_inds[j]
            )
        ]
        # order rays
        if radar_aux.scan_type == "ppi":
            ind = np.argsort(radar_aux.azimuth["data"][ind_rays_sweep])
            ind_rays_sweep = ind_rays_sweep[ind]
            # avoid large gaps in data
            azimuths = radar_aux.azimuth["data"][ind_rays_sweep]
            azi_steps = azimuths[1:] - azimuths[:-1]
            ind_gap = np.where(azi_steps > 2 * np.median(azi_steps))[0]
            if ind_gap.size > 0:
                ind_rays_sweep = np.append(
                    ind_rays_sweep[ind_gap[0] + 1 :], ind_rays_sweep[: ind_gap[0] + 1]
                )
            ind_rays_aux.extend(ind_rays_sweep)
        else:
            ind = np.argsort(radar_aux.elevation["data"][ind_rays_sweep])
            ind_rays_aux.extend(ind_rays_sweep[ind])

        rays_in_sweep = ind_rays_sweep.size
        radar_aux.rays_per_sweep["data"][j] = rays_in_sweep
        if j == 0:
            radar_aux.sweep_start_ray_index["data"][j] = 0
        else:
            radar_aux.sweep_start_ray_index["data"][j] = int(
                radar_aux.sweep_end_ray_index["data"][j - 1] + 1
            )
        radar_aux.sweep_end_ray_index["data"][j] = (
            radar_aux.sweep_start_ray_index["data"][j] + rays_in_sweep - 1
        )
        nrays += rays_in_sweep

    ind_rays = np.array(ind_rays_aux)

    # Update metadata
    radar_aux.range["data"] = radar_aux.range["data"][ind_rng]
    radar_aux.time["data"] = radar_aux.time["data"][ind_rays]
    radar_aux.azimuth["data"] = radar_aux.azimuth["data"][ind_rays]
    radar_aux.elevation["data"] = radar_aux.elevation["data"][ind_rays]
    radar_aux.init_gate_x_y_z()
    radar_aux.init_gate_longitude_latitude()
    radar_aux.init_gate_altitude()
    radar_aux.nrays = nrays
    radar_aux.ngates = ind_rng.size

    if radar_aux.instrument_parameters is not None:
        if "nyquist_velocity" in radar_aux.instrument_parameters:
            radar_aux.instrument_parameters["nyquist_velocity"][
                "data"
            ] = radar_aux.instrument_parameters["nyquist_velocity"]["data"][ind_rays]
        if "pulse_width" in radar_aux.instrument_parameters:
            radar_aux.instrument_parameters["pulse_width"][
                "data"
            ] = radar_aux.instrument_parameters["pulse_width"]["data"][ind_rays]
        if "number_of_pulses" in radar_aux.instrument_parameters:
            radar_aux.instrument_parameters["number_of_pulses"][
                "data"
            ] = radar_aux.instrument_parameters["number_of_pulses"]["data"][ind_rays]

    # Get new fields
    if field_names is None:
        radar_aux.fields = dict()
    else:
        fields_aux = copy.deepcopy(radar_aux.fields)
        radar_aux.fields = dict()
        for field_name in field_names:
            if field_name not in fields_aux:
                print("Field " + field_name + " not available")
                continue

            fields_aux[field_name]["data"] = fields_aux[field_name]["data"][:, ind_rng]
            fields_aux[field_name]["data"] = fields_aux[field_name]["data"][ind_rays, :]
            radar_aux.add_field(field_name, fields_aux[field_name])

    return radar_aux


def join_radar(radar1, radar2):
    """
    Combine two radar instances into one.

    Parameters
    ----------
    radar1 : Radar
        Radar object.
    radar2 : Radar
        Radar object.

    """
    # must have same gate spacing
    new_radar = copy.deepcopy(radar1)
    new_radar.azimuth["data"] = np.append(
        radar1.azimuth["data"], radar2.azimuth["data"]
    )
    new_radar.elevation["data"] = np.append(
        radar1.elevation["data"], radar2.elevation["data"]
    )
    new_radar.fixed_angle["data"] = np.append(
        radar1.fixed_angle["data"], radar2.fixed_angle["data"]
    )
    new_radar.sweep_number["data"] = np.append(
        radar1.sweep_number["data"], radar2.sweep_number["data"]
    )
    new_radar.sweep_start_ray_index["data"] = np.append(
        radar1.sweep_start_ray_index["data"],
        radar2.sweep_start_ray_index["data"] + radar1.nrays,
    )
    new_radar.sweep_end_ray_index["data"] = np.append(
        radar1.sweep_end_ray_index["data"],
        radar2.sweep_end_ray_index["data"] + radar1.nrays,
    )
    new_radar.nsweeps += radar2.nsweeps
    new_radar.sweep_mode["data"] = np.append(
        radar1.sweep_mode["data"], radar2.sweep_mode["data"]
    )
    if (radar1.rays_are_indexed is not None) and (radar2.rays_are_indexed is not None):
        new_radar.rays_are_indexed["data"] = np.append(
            radar1.rays_are_indexed["data"], radar2.rays_are_indexed["data"]
        )
    else:
        new_radar.rays_are_indexed = None

    if new_radar.instrument_parameters is not None:
        if "nyquist_velocity" in new_radar.instrument_parameters:
            new_radar.instrument_parameters["nyquist_velocity"]["data"] = np.append(
                radar1.instrument_parameters["nyquist_velocity"]["data"],
                radar2.instrument_parameters["nyquist_velocity"]["data"],
            )
        if "pulse_width" in new_radar.instrument_parameters:
            new_radar.instrument_parameters["pulse_width"]["data"] = np.append(
                radar1.instrument_parameters["pulse_width"]["data"],
                radar2.instrument_parameters["pulse_width"]["data"],
            )
        if "number_of_pulses" in new_radar.instrument_parameters:
            new_radar.instrument_parameters["number_of_pulses"]["data"] = np.append(
                radar1.instrument_parameters["number_of_pulses"]["data"],
                radar2.instrument_parameters["number_of_pulses"]["data"],
            )
        if "prt" in new_radar.instrument_parameters:
            new_radar.instrument_parameters["prt"]["data"] = np.append(
                radar1.instrument_parameters["prt"]["data"],
                radar2.instrument_parameters["prt"]["data"],
            )

    if (radar1.ray_angle_res is not None) and (radar2.ray_angle_res is not None):
        new_radar.ray_angle_res["data"] = np.append(
            radar1.ray_angle_res["data"], radar2.ray_angle_res["data"]
        )
    else:
        new_radar.ray_angle_res = None

    if len(radar1.range["data"]) >= len(radar2.range["data"]):
        new_radar.range["data"] = radar1.range["data"]
    else:
        new_radar.range["data"] = radar2.range["data"]
    new_radar.ngates = len(new_radar.range["data"])

    if (radar1.target_scan_rate is not None) and (radar2.target_scan_rate is not None):
        new_radar.target_scan_rate["data"] = np.append(
            radar1.target_scan_rate["data"], radar2.target_scan_rate["data"]
        )
    else:
        new_radar.target_scan_rate = None

    # to combine times we need to reference them to a standard
    # for this we'll use epoch time
    r1num = datetime_utils.datetimes_from_radar(radar1, epoch=True)
    r2num = datetime_utils.datetimes_from_radar(radar2, epoch=True)
    new_radar.time["data"] = date2num(
        np.append(r1num, r2num), datetime_utils.EPOCH_UNITS
    )
    new_radar.time["units"] = datetime_utils.EPOCH_UNITS
    new_radar.nrays = len(new_radar.time["data"])

    fields_to_remove = []
    for var in new_radar.fields.keys():
        # if the field is present in both radars combine both fields
        # otherwise remove it from new radar
        if var in radar1.fields and var in radar2.fields:
            sh1 = radar1.fields[var]["data"].shape
            sh2 = radar2.fields[var]["data"].shape
            new_field_shape = (sh1[0] + sh2[0], max(sh1[1], sh2[1]))
            new_field = np.ma.masked_all(new_field_shape)
            new_field.set_fill_value(get_fillvalue())
            new_field[0 : sh1[0], 0 : sh1[1]] = radar1.fields[var]["data"]
            new_field[sh1[0] :, 0 : sh2[1]] = radar2.fields[var]["data"]
            new_radar.fields[var]["data"] = new_field
        else:
            print("Field " + var + " not present in both radars")
            fields_to_remove.append(var)

    if fields_to_remove:
        for field_name in fields_to_remove:
            new_radar.fields.pop(field_name, None)

    # radar locations
    # TODO moving platforms - any more?
    if (
        len(radar1.latitude["data"])
        == 1 & len(radar2.latitude["data"])
        == 1 & len(radar1.longitude["data"])
        == 1 & len(radar2.longitude["data"])
        == 1 & len(radar1.altitude["data"])
        == 1 & len(radar2.altitude["data"])
        == 1
    ):
        lat1 = float(radar1.latitude["data"])
        lon1 = float(radar1.longitude["data"])
        alt1 = float(radar1.altitude["data"])
        lat2 = float(radar2.latitude["data"])
        lon2 = float(radar2.longitude["data"])
        alt2 = float(radar2.altitude["data"])

        if (lat1 != lat2) or (lon1 != lon2) or (alt1 != alt2):
            ones1 = np.ones(len(radar1.time["data"]), dtype="float32")
            ones2 = np.ones(len(radar2.time["data"]), dtype="float32")
            new_radar.latitude["data"] = np.append(ones1 * lat1, ones2 * lat2)
            new_radar.longitude["data"] = np.append(ones1 * lon1, ones2 * lon2)
            new_radar.latitude["data"] = np.append(ones1 * alt1, ones2 * alt2)
        else:
            new_radar.latitude["data"] = radar1.latitude["data"]
            new_radar.longitude["data"] = radar1.longitude["data"]
            new_radar.altitude["data"] = radar1.altitude["data"]

    else:
        new_radar.latitude["data"] = np.append(
            radar1.latitude["data"], radar2.latitude["data"]
        )
        new_radar.longitude["data"] = np.append(
            radar1.longitude["data"], radar2.longitude["data"]
        )
        new_radar.altitude["data"] = np.append(
            radar1.altitude["data"], radar2.altitude["data"]
        )

    return new_radar


def image_mute_radar(radar, field, mute_field, mute_threshold, field_threshold=None):
    """
    This function will split a field based on thresholds from another field.

    Specifically, it was designed to separate areas of reflectivity where
    the correlation coefficient is less than a certain threshold to discern
    melting precipitation.

    Author: Laura Tomkins (@lauratomkins)

    Parameters
    ----------
    radar : Radar
        Radar instance which provides the fields for muting.
    field : str
        Name of field to image mute.
    mute_field : str
        Name of field to image mute by.
    mute_threshold : float
        Threshold value to mute by.
    field_threshold : float
        Additional threshold to mask.

    Returns
    -------
    radar : Radar
        Radar object with 2 new fields from input field, one muted and one not muted.

    References
    ----------
    Tomkins, L. M., Yuter, S. E., Miller, M. A., and Allen, L. R., 2022:
    Image muting of mixed precipitation to improve identification of regions
    of heavy snow in radar data. Atmos. Meas. Tech., 15, 5515–5525,
    https://doi.org/10.5194/amt-15-5515-2022

    """
    # add checks for field availability
    if field not in radar.fields.keys():
        raise KeyError("Failed - ", field, " field to mute not found in Radar object.")

    if mute_field not in radar.fields.keys():
        raise KeyError(
            "Failed - ", mute_field, " field to mute by not found in Radar object."
        )

    # get data from fields
    data_to_mute = radar.fields[field]["data"]
    data_mute_by = radar.fields[mute_field]["data"]

    # create filters
    # field_filter is used if user wants to use additional criteria in the
    # original field
    if field_threshold is not None:
        field_filter = data_to_mute >= field_threshold
    else:
        field_filter = None

    # mute_filter will be the primary filter for determining muted regions
    mute_filter = data_mute_by <= mute_threshold

    # mute_mask is the combined filter
    if field_filter is None:
        mute_mask = mute_filter
    else:
        mute_mask = mute_filter & field_filter

    # break up the field into muted regions and non muted regions
    non_muted_field = np.ma.masked_where(mute_mask, data_to_mute)
    non_muted_field = np.ma.masked_invalid(non_muted_field)

    muted_field = np.ma.masked_where(~mute_mask, data_to_mute)
    muted_field = np.ma.masked_invalid(muted_field)

    # add fields to a dictionary and save to radar object
    non_muted_dict = radar.fields[field].copy()
    non_muted_dict["data"] = non_muted_field
    non_muted_dict["long_name"] = "Non-muted " + field
    radar.add_field("nonmuted_" + field, non_muted_dict)

    muted_dict = radar.fields[field].copy()
    muted_dict["data"] = muted_field
    muted_dict["long_name"] = "Muted " + field
    radar.add_field("muted_" + field, muted_dict)
    return radar
