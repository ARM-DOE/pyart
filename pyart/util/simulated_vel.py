"""
pyart.util.simulated_vel
========================

Function for creating simulated velocity fields.

.. autosummary::
    :toctree: generated/

    simulated_vel_from_profile

"""

import numpy as np
from scipy.interpolate import interp1d

from ..config import get_metadata, get_field_name


def simulated_vel_from_profile(
        radar, profile, interp_kind='linear', sim_vel_field=None):
    """
    Create simulated radial velocities from a profile of horizontal winds.

    Parameters
    ----------
    radar : Radar
        Radar instance which provides the scanning parameters for the
        simulated radial velocities.
    profile : HorizontalWindProfile
        Profile of horizontal winds.
    interp_kind : str, optional
        Specifies the kind of interpolation used to determine the winds at a
        given height. Must be one of 'linear', 'nearest', 'zero', 'slinear',
        'quadratic', or 'cubic'. The the documentation for the SciPy
        scipy.interpolate.interp1d function for descriptions.
    sim_vel_field : str, optional
        Name to use for the simulated velocity field metadata. None will use
        the default field name from the Py-ART configuration file.

    Returns
    -------
    sim_vel : dict
        Dictionary containing a radar field of simulated radial velocities.

    """
    # parse parameters
    if sim_vel_field is None:
        sim_vel_field = get_field_name('simulated_velocity')

    # radar parameters
    azimuths = np.deg2rad(radar.azimuth['data']).reshape(-1, 1)
    elevations = np.deg2rad(radar.elevation['data']).reshape(-1, 1)
    gate_altitudes = radar.gate_altitude['data']

    if isinstance(gate_altitudes, np.ma.MaskedArray):
        gate_altitudes = gate_altitudes.filled(np.nan)

    # prepare wind profile for interpolation
    if isinstance(profile.height, np.ma.MaskedArray):
        height = profile.height.filled(np.nan)
    else:
        height = profile.height

    height_is_not_nan = ~np.isnan(height)
    winds = np.empty((2, len(height)), dtype=np.float64)
    if isinstance(profile.u_wind, np.ma.MaskedArray):
        winds[0] = profile.u_wind.filled(np.nan)
    else:
        winds[0] = profile.u_wind

    if isinstance(profile.v_wind, np.ma.MaskedArray):
        winds[1] = profile.v_wind.filled(np.nan)
    else:
        winds[1] = profile.v_wind

    wind_is_not_nan = np.logical_and(~np.isnan(winds[0]), ~np.isnan(winds[1]))
    no_nans = np.logical_and(height_is_not_nan, wind_is_not_nan)
    height = height[no_nans]
    winds[0] = winds[0][no_nans]
    winds[1] = winds[1][no_nans]
    wind_interp = interp1d(
        height, winds, kind=interp_kind, bounds_error=False)

    # interpolated wind speeds at all gates altitudes
    gate_winds = wind_interp(gate_altitudes)
    gate_u = np.ma.masked_invalid(gate_winds[0])
    gate_v = np.ma.masked_invalid(gate_winds[1])

    # calculate the radial velocity for all gates
    radial_vel = (gate_u * np.sin(azimuths) * np.cos(elevations) +
                  gate_v * np.cos(azimuths) * np.cos(elevations))

    sim_vel = get_metadata(sim_vel_field)
    sim_vel['data'] = radial_vel
    return sim_vel
