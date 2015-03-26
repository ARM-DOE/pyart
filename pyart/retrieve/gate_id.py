"""
pyart.retrieve.gate_id
=========================

gate_id

.. autosummary::
    :toctree: generated/

    gate_id

"""

import numpy as np
import netCDF4
from scipy import interpolate, ndimage


from ..config import get_fillvalue, get_field_name, get_metadata
from ..io.common import radar_coords_to_cart

def map_profile_to_gates(profile, heights,
                         radar, info_dict = None, toa= None):
    """
    Given a profile of a variable map it to the gates of radar assuming 4/3Re
    ----------
    profile : numpy array
        Array to map
    heights : numpy array
        monotonically increasing heights in meters of shape(profile)
    radar : Radar
        Radar to map to

    Returns
    -------
    new_fields : (height_dict,profile_dict)
        twople of field dictionaries containing the height of the gates and
        the profile
        interpolated onto the radar gates

    Other Parameters
    ----------------
    info_dict: dictionary
        Dictionary of metadata for field, Defaults to None
    toa: float
        Top of atmosphere, where to use profile up to. If
        None check for mask and use lowest
        element, if no mask uses whole profile.

    References
    ----------

    """
    #first retrieve the Z coordinates of the radar gates
    rg, azg = np.meshgrid(radar.range['data'], radar.azimuth['data'] )
    rg, eleg = np.meshgrid(radar.range['data'], radar.elevation['data'])
    x,y,z = radar_coords_to_cart(rg /1000.0, azg, eleg)

    #set up the interpolation
    if toa == None:
        if 'mask' in dir(my_profile['temp']):
            toa = np.where(my_profile['temp'].mask)[0].min()
        else:
            toa = -1

    f_interp = interpolate.interp1d(heights[0:toa], profile[0:toa],
            bounds_error = False,
            fill_value = pyart.config.get_fillvalue())

    #interpolate
    fld = np.ma.masked_equal(f_interp(z + radar.altitude['data'][0]),
                             get_fillvalue())

    #populate field dictionaries

    height_dict = get_metadata(radar.fields.keys()[0])
    height_dict['data'] = z + radar.altitude['data'][0]
    height_dict['long_name'] = 'Height of radar beam agl'
    height_dict['standard_name'] = 'height'
    height_dict['valid_min'] = -100.0
    height_dict['valid_max'] = 400000.0
    height_dict['units'] = 'meters'

    profile_dict = get_metadata(radar.fields.keys()[0])
    profile_dict['data'] = fld
    if info_dict != None:
        for key in info_dict.keys():
            profile_dict[key] = info_dict[key]
    else:
        profile_dict['long_name'] = 'Interpolated profile'
        profile_dict['standard_name'] = 'interp_profile'
        profile_dict['valid_min'] = -9999
        profile_dict['valid_max'] = 9999
        profile_dict['units'] = 'unknown'

    return height_dict, profile_dict

