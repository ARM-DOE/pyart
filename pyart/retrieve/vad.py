"""
Retrieval of VADs from a radar object.

"""

import numpy as np

from ..config import get_field_name
from ..core import HorizontalWindProfile


def vad_michelson(
    radar,
    vel_field=None,
    z_want=None,
    gatefilter=None,
    valid_ray_min=16,
    max_speed_error=2.0,
):
    """
    Velocity azimuth display.

    Creates a VAD object containing U Wind, V Wind and height that
    can then be used to plot and produce the velocity azimuth display.

    Parameters
    ----------
    radar : Radar
        Radar object used.
    vel_field : string, optional
        Velocity field to use for VAD calculation.
    z_want : array, optional
        Heights for where to sample vads from.
        None will result in np.linespace(0, 10000, 100).
    gatefilter : GateFilter, optional
        A GateFilter indicating radar gates that should be excluded
        from the import vad calculation.
    valid_ray_min : int, optional
        Minimum number of rays with valid velocity a gate needs to be used.
        Gates with fewer rays give noisy fits and are left out. Raising it
        smooths the profile but leaves more heights empty. Default is 16,
        the same as vad_browning.
    max_speed_error : float, optional
        Largest standard error of the fitted wind speed, in m/s, for a gate
        to be used. The error grows when valid rays are few, bunched in one
        part of the circle, or noisy, so this rejects unreliable fits.
        Gates are then averaged into height bins weighted by 1 / error**2.
        None turns the check off. Default is 2.0.

    Returns
    -------
    vad : HorizontalWindProfile
        A velocity azimuth display object containing height, speed, direction,
        u_wind, v_wind from a radar object.

    References
    ----------
    Michelson, D. B., Andersson, T., Koistinen, J., Collier, C. G., Riedl, J.,
    Szturc, J., Gjertsen, U., Nielsen, A. and Overgaard, S. (2000) BALTEX Radar
    Data Centre Products and their Methodologies. In SMHI Reports. Meteorology
    and Climatology. Swedish Meteorological and Hydrological Institute,
    Norrkoping.

    """
    speeds = []
    angles = []
    errors = []
    heights = []

    # Pulling z data from radar
    z_gate_data = radar.gate_z["data"]

    # Setting parameters
    if z_want is None:
        z_want = np.linspace(0, 1000, 100)

    # Parse field parameters
    if vel_field is None:
        radar.check_field_exists("velocity")
        vel_field = get_field_name("velocity")

    # Selecting what velocity data to use based on gatefilter
    if gatefilter is not None:
        velocities = np.ma.masked_where(
            gatefilter.gate_excluded, radar.fields[vel_field]["data"]
        )
    else:
        velocities = radar.fields[vel_field]["data"]

    # Getting radar sweep index values
    for i in range(len(radar.sweep_start_ray_index["data"])):
        index_start = radar.sweep_start_ray_index["data"][i]
        index_end = radar.sweep_end_ray_index["data"][i]
        if not (index_end - index_start) % 2 == 0:
            index_end = index_end - 1

        used_velocities = velocities[index_start:index_end]
        azimuth = radar.azimuth["data"][index_start:index_end]
        elevation = radar.fixed_angle["data"][i]

        # Calculating speed and angle
        speed, angle, speed_error = _vad_calculation_m(
            used_velocities, azimuth, elevation, valid_ray_min, max_speed_error
        )

        print("max height", z_gate_data[index_start, :].max(), "meters")

        # Filling empty arrays with data
        speeds.append(speed)
        angles.append(angle)
        errors.append(speed_error)
        heights.append(z_gate_data[index_start, :])

    # Combining arrays and sorting
    speed_array = np.concatenate(speeds)
    angle_array = np.concatenate(angles)
    height_array = np.concatenate(heights)
    arg_order = height_array.argsort()
    speed_ordered = speed_array[arg_order]
    height_ordered = height_array[arg_order]
    angle_ordered = angle_array[arg_order]
    # Weight each gate by its confidence; a floor keeps exact fits finite
    error_ordered = np.concatenate(errors)[arg_order]
    weight_ordered = 1.0 / np.maximum(error_ordered, 0.01) ** 2

    # Calculating U and V wind
    u_ordered, v_ordered = _sd_to_uv(speed_ordered, angle_ordered)
    u_mean = _interval_mean(u_ordered, height_ordered, z_want, weight_ordered)
    v_mean = _interval_mean(v_ordered, height_ordered, z_want, weight_ordered)
    vad = HorizontalWindProfile.from_u_and_v(z_want, u_mean, v_mean)
    return vad


def _vad_calculation_m(
    velocity_field, azimuth, elevation, valid_ray_min=0, max_speed_error=None
):
    """Fit v = c + a sin(az) + b cos(az) at each gate by least squares.

    Returns the horizontal speed, the angle and the standard error of the
    speed. Gates with fewer than valid_ray_min valid rays, or a speed error
    above max_speed_error, are set to nan.
    """
    # Masked gates become nan so they are left out of the fit
    velocity_field = np.ma.filled(
        np.ma.asarray(velocity_field, dtype=np.float64), np.nan
    )
    valid = ~np.isnan(velocity_field)
    velocity = np.where(valid, velocity_field, 0.0)
    n_valid = valid.sum(axis=0)

    # Design matrix, one row per ray: [1, sin(az), cos(az)]
    az = np.deg2rad(azimuth)
    design = np.stack([np.ones_like(az), np.sin(az), np.cos(az)], axis=1)

    # Normal equations for every gate at once, using valid rays only
    normal = np.einsum("ij,ik,il->jkl", valid.astype(np.float64), design, design)
    rhs = np.einsum("ij,ik->jk", velocity, design)
    solvable = (n_valid > 3) & (np.abs(np.linalg.det(normal)) > 1e-9)
    normal[~solvable] = np.eye(3)
    normal_inv = np.linalg.inv(normal)
    params = np.einsum("jkl,jl->jk", normal_inv, rhs)

    # Standard errors from the scatter left after the fit
    residual = np.where(valid, velocity - design @ params.T, 0.0)
    variance = (residual**2).sum(axis=0) / np.maximum(n_valid - 3, 1)
    cov = variance[:, np.newaxis, np.newaxis] * normal_inv

    a_value = params[:, 1]
    b_value = params[:, 2]
    amplitude = np.hypot(a_value, b_value)
    amplitude_var = (
        a_value**2 * cov[:, 1, 1]
        + b_value**2 * cov[:, 2, 2]
        + 2 * a_value * b_value * cov[:, 1, 2]
    ) / np.maximum(amplitude**2, np.finfo(np.float64).tiny)

    cos_elevation = np.cos(np.deg2rad(elevation))
    speed = amplitude / cos_elevation
    speed_error = np.sqrt(amplitude_var) / cos_elevation
    angle = np.arctan2(a_value, b_value)

    # Leave out gates without enough valid rays or with an unreliable fit
    rejected = ~solvable | (n_valid < valid_ray_min)
    if max_speed_error is not None:
        rejected |= speed_error > max_speed_error
    speed[rejected] = np.nan
    angle[rejected] = np.nan
    speed_error[rejected] = np.nan
    return speed, angle, speed_error


def _interval_mean(data, current_z, wanted_z, weights=None):
    """Find the mean of data indexed by current_z
    at wanted_z on intervals wanted_z+/- delta
    wanted_z. Nan values are skipped, and weights, when given,
    make it a weighted mean."""
    delta = wanted_z[1] - wanted_z[0]
    if weights is None:
        weights = np.ones_like(data)
    # argmin returns the first of tied gate heights on every platform,
    # unlike argsort, whose tie order depends on the sort implementation
    pos_lower = [
        np.argmin((current_z - (wanted_z[i] - delta / 2.0)) ** 2)
        for i in range(len(wanted_z))
    ]
    pos_upper = [
        np.argmin((current_z - (wanted_z[i] + delta / 2.0)) ** 2)
        for i in range(len(wanted_z))
    ]
    mean_values = np.full(len(wanted_z), np.nan)
    for i in range(len(wanted_z)):
        values = data[pos_lower[i] : pos_upper[i]]
        value_weights = weights[pos_lower[i] : pos_upper[i]]
        use = np.isfinite(values) & np.isfinite(value_weights)
        if use.any():
            mean_values[i] = np.average(values[use], weights=value_weights[use])
    return mean_values


def _sd_to_uv(speed, direction):
    """Takes speed and direction to create u_mean and v_mean."""
    return (np.sin(direction) * speed), (np.cos(direction) * speed)


def vad_browning(
    radar,
    velocity,
    z_want=None,
    valid_ray_min=16,
    gatefilter=None,
    window=2,
    weight="equal",
):
    """
    Velocity azimuth display.
    Note: This code uses only one sweep. Before using the
    velocity_azimuth_display function, use, for example:
    one_sweep_radar = radar.extract_sweeps([0])

    Parameters
    ----------
    radar : Radar
        Radar object used.
    velocity : string
        Velocity field to use for VAD calculation.

    Other Parameters
    ----------------
    z_want : array
        Array of desired heights to be sampled for the vad
        calculation.
    valid_ray_min : int
        Amount of rays required to include that level in
        the VAD calculation.
    gatefilter : GateFilter
        A GateFilter indicating radar gates that should be excluded when
        from the import vad calculation.
    window : int
        Value to use for window when determining new values in the
        _Averag1D function.
    weight : string
        A string to indicate weighting method to use. 'equal' for
        equal weighting when interpolating or 'idw' for inverse
        distribution squared weighting for interpolating.
        Default is 'equal'.

    Returns
    -------
    height : array
        Heights in meters above sea level at which horizontal winds were
        sampled.
    speed : array
        Horizontal wind speed in meters per second at each height.
    direction : array
        Horizontal wind direction in degrees at each height.
    u_wind : array
        U-wind mean in meters per second.
    v_wind : array
        V-wind mean in meters per second.

    Reference
    ----------
    K. A. Browning and R. Wexler, 1968: The Determination
    of Kinematic Properties of a Wind Field Using Doppler
    Radar. J. Appl. Meteor., 7, 105–113

    """
    velocities = radar.fields[velocity]["data"]
    if gatefilter is not None:
        velocities = np.ma.masked_where(gatefilter.gate_excluded, velocities)
    azimuths = radar.azimuth["data"][:]
    elevation = radar.fixed_angle["data"][0]

    u_wind, v_wind = _vad_calculation_b(velocities, azimuths, elevation, valid_ray_min)
    bad = np.logical_or(np.isnan(u_wind), np.isnan(v_wind))
    good_u_wind = u_wind[~bad]
    good_v_wind = v_wind[~bad]
    radar_height = radar.gate_z["data"][0]
    good_height = radar_height[~bad]
    if z_want is None:
        z_want = np.linspace(0, 1000, 100)[:50]
    try:
        print("max height", np.max(good_height), " meters")
        print("min height", np.min(good_height), " meters")
    except ValueError:
        raise ValueError(
            "Not enough data in this radar sweep " "for a vad calculation."
        )

    u_interp = _Average1D(
        good_height, good_u_wind, z_want[1] - z_want[0] / window, weight
    )
    v_interp = _Average1D(
        good_height, good_v_wind, z_want[1] - z_want[0] / window, weight
    )

    u_wanted = u_interp(z_want)
    v_wanted = v_interp(z_want)
    u_wanted = np.ma.masked_equal(u_wanted, 99999.0)
    v_wanted = np.ma.masked_equal(v_wanted, 99999.0)

    vad = HorizontalWindProfile.from_u_and_v(z_want, u_wanted, v_wanted)
    return vad


def _vad_calculation_b(velocities, azimuths, elevation, valid_ray_min):
    """Calculates VAD for a scan and returns u_mean and
    v_mean. velocities is a 2D array, azimuths is a 1D
    array, elevation is a number.
    Note:
    We need to solve: Ax = b
    where:
    A = [sum_sin_squared_az, sum_sin_cos_az    ] = [a, b]
        [sum_sin_cos_az,     sum_cos_squared_az]   [c, d]
    b = [sum_sin_vel_dev] = [b_1]
        [sum_cos_vel_dev]   [b_2]
    The solution to this is:
    x = A-1 * b
    A-1 is:
     1    [ d,  -b ]
    --- * [ -c,  a ]
    |A|
    and the determinate, det is: det = a*d - b*c
    Therefore the elements of x are:
    x_1 = (d* b_1  + -b * b_2) / det = (d*b_1 - b*b_2) / det
    x_2 = (-c * b_1 +  a * b_2) / det = (a*b_2 - c*b_1) / det
    """
    velocities = velocities.filled(np.nan)
    shape = velocities.shape
    _, nbins = velocities.shape

    invalid = np.isnan(velocities)
    valid_rays_per_gate = np.sum(~np.isnan(velocities), axis=0)
    too_few_valid_rays = valid_rays_per_gate < valid_ray_min
    invalid[:, too_few_valid_rays] = True

    sin_az = np.sin(np.deg2rad(azimuths))
    cos_az = np.cos(np.deg2rad(azimuths))
    sin_az = np.repeat(sin_az, nbins).reshape(shape)
    cos_az = np.repeat(cos_az, nbins).reshape(shape)
    sin_az[invalid] = np.nan
    cos_az[invalid] = np.nan

    mean_velocity_per_gate = np.nanmean(velocities, axis=0).reshape(1, -1)
    velocity_deviation = velocities - mean_velocity_per_gate

    sum_cos_vel_dev = np.nansum(cos_az * velocity_deviation, axis=0)
    sum_sin_vel_dev = np.nansum(sin_az * velocity_deviation, axis=0)

    sum_sin_cos_az = np.nansum(sin_az * cos_az, axis=0)
    sum_sin_squared_az = np.nansum(sin_az**2, axis=0)
    sum_cos_squared_az = np.nansum(cos_az**2, axis=0)

    # The A matrix
    a = sum_sin_squared_az
    b = sum_sin_cos_az
    c = sum_sin_cos_az
    d = sum_cos_squared_az

    # The b vector
    b_1 = sum_sin_vel_dev
    b_2 = sum_cos_vel_dev

    # solve for the x vector
    determinant = a * d - b * c
    x_1 = (d * b_1 - b * b_2) / determinant
    x_2 = (a * b_2 - c * b_1) / determinant

    # calculate horizontal components of winds
    elevation_scale = 1 / np.cos(np.deg2rad(elevation))
    u_mean = x_1 * elevation_scale
    v_mean = x_2 * elevation_scale
    return u_mean, v_mean


def _inverse_dist_squared(dist):
    """Obtaining distance weights by using distance weighting
    interpolation, using the inverse distance-squared relationship.
    """
    weights = 1 / (dist * dist)
    weights[np.isnan(weights)] = 99999.0
    return weights


class _Average1D:
    """Used to find the nearest gate height and horizontal wind
    value with respect to the user's desired height."""

    def __init__(self, x, y, window, weight, fill_value=99999.0):
        sort_idx = np.argsort(x)
        self.x_sorted = x[sort_idx]
        self.y_sorted = y[sort_idx]
        self.window = window
        self.fill_value = fill_value

        if weight == "equal":
            self.weight_func = lambda x: None
        elif weight == "idw":
            self.weight_func = _inverse_dist_squared
        elif callable(weight):
            self.weight_func = weight
        else:
            raise ValueError("Invalid weight argument:", weight)

    def __call__(self, x_new, window=None):
        if window is None:
            window = self.window

        y_new = np.zeros_like(x_new, dtype=self.y_sorted.dtype)
        for i, center in enumerate(x_new):
            bottom = center - window
            top = center + window
            start = np.searchsorted(self.x_sorted, bottom)
            stop = np.searchsorted(self.x_sorted, top)

            x_in_window = self.x_sorted[start:stop]
            y_in_window = self.y_sorted[start:stop]
            if len(x_in_window) == 0:
                y_new[i] = self.fill_value
            else:
                distances = x_in_window - center
                weights = self.weight_func(distances)
                y_new[i] = np.average(y_in_window, weights=weights)
        return y_new
