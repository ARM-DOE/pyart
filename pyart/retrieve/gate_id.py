"""
pyart.retrieve.gate_id
======================

.. autosummary::
    :toctree: generated/

    map_profile_to_gates
    fetch_radar_time_profile

"""

try:
    from netCDF4 import num2date, datetime
except ImportError:
    from cftime import num2date, datetime

import numpy as np
from scipy import interpolate

from ..config import get_fillvalue, get_metadata, get_field_name
from ..core.transforms import antenna_to_cartesian


def map_profile_to_gates(profile, heights, radar, toa=None,
                         profile_field=None, height_field=None):
    """
    Given a profile of a variable map it to the gates of radar assuming 4/3Re.

    Parameters
    ----------
    profile : array
        Profile array to map.
    heights : array
        Monotonically increasing heights in meters with same shape as profile.
    radar : Radar
        Radar to map to.
    toa : float, optional
        Top of atmosphere, where to use profile up to. If None check for
        mask and use lowest element, if no mask uses whole profile.
    height_field : str, optional
        Name to use for height field metadata. None will use the default field
        name from the Py-ART configuration file.
    profile_field : str, optional
        Name to use for interpolate profile field metadata. None will use the
        default field name from the Py-ART configuration file.

    Returns
    -------
    height_dict, profile_dict : dict
        Field dictionaries containing the height of the gates and the profile
        interpolated onto the radar gates.

    """
    # retrieve the Z coordinates of the radar gates
    rg, azg = np.meshgrid(radar.range['data'], radar.azimuth['data'])
    rg, eleg = np.meshgrid(radar.range['data'], radar.elevation['data'])
    _, _, z = antenna_to_cartesian(rg / 1000.0, azg, eleg)

    # Check that z is not a MaskedArray
    if isinstance(z, np.ma.MaskedArray):
        z = z.filled(np.NaN)

    # find toa is not provided
    if toa is None:
        ismasked = np.where(np.ma.getmaskarray(profile))[0]
        if len(ismasked) == 0:
            toa = None
        else:
            toa = ismasked.min()

    # interpolate
    f_interp = interpolate.interp1d(
        heights[:toa], profile[:toa], bounds_error=False,
        fill_value=get_fillvalue())
    fld = np.ma.masked_equal(f_interp(z + radar.altitude['data'][0]),
                             get_fillvalue())

    # populate field dictionaries
    if height_field is None:
        height_field = get_field_name('height')
    height_dict = get_metadata(height_field)
    height_dict['data'] = z + radar.altitude['data'][0]

    if profile_field is None:
        profile_field = get_field_name('interpolated_profile')
    profile_dict = get_metadata(profile_field)
    profile_dict['data'] = fld

    return height_dict, profile_dict


def fetch_radar_time_profile(sonde_dset, radar, time_key='time',
                             height_key='height', nvars=None):
    """
    Extract the correct profile from a interpolated sonde.

    This is an ARM specific method which extract the correct profile out of
    netCDF Variables from a Interpolated Sonde VAP for the volume start time
    of a radar object.

    Parameters
    ----------
    sonde_dset : Dataset
        Interpolate sonde Dataset.
    radar : Radar
        Radar object from which the nearest profile will be found.
    time_key : string, optional
        Key to find a CF startard time variable.
    height_key : string, optional
        Key to find profile height data.
    nvars : list, optional
        NetCDF variable to generated profiles for. If None (the default) all
        variables with dimension of time, height will be found in ncvars.

    Returns
    -------
    return_dic : dict
        Profiles at the start time of the radar.

    """
    ncvars = sonde_dset.variables
    if nvars is None:
        time_height_shape = (len(ncvars[time_key]), len(ncvars[height_key]))
        nvars = [k for k, v in ncvars.items() if v.shape == time_height_shape]

    radar_start = num2date(radar.time['data'][0], radar.time['units'])
    radar_day_start = datetime(radar_start.year, radar_start.month,
                               radar_start.day)
    seconds_since_start_of_day = (radar_start - radar_day_start).seconds
    time_index = abs(ncvars[time_key][:] - seconds_since_start_of_day).argmin()

    return_dic = dict([(key, ncvars[key][time_index, :]) for key in nvars])
    return_dic[height_key] = ncvars[height_key][:]
    return return_dic
