"""
pyart.retrieve.simple_moment_calculations
=========================

steiner_conv_strat

.. autosummary::
    :toctree: generated/

    calculate_snr_from_reflectivity
"""

import numpy as np

from ..config import get_metadata
from ..io.common import radar_coords_to_cart

def calculate_snr_from_reflectivity(radar,
        reflectivity_key = None,
        toa = None):
    """
    Back out the Signal to noise ratio, in dB, from the reflectivity field
    ----------
    radar : Radar
        Radar object

    Returns
    -------
    snr : field dictionary
        Field dictionary containing the signal to noise ratio

    Other Parameters
    ----------------
    reflectivity_key : string
        In case reflectivity is not in reflectivity
    toa : float
        height above which to take noise floor measurements

    References
    ----------

    """
    if reflectivity_key == None:
        reflectivity_key = 'reflectivity'

    if toa == None:
        toa = 25000.


    range_grid = np.meshgrid(radar.range['data'],
                             np.ma.ones(radar.time['data'].shape))[0] +1.0

    # remove range scale.. This is basically the radar constant scaled dBm
    pseudo_power = radar.fields[reflectivity_key]['data']\
            -20.0*np.log10(range_grid/1000.0)

    #Noise floor estimate
    #25km.. should be no scatterers, not even planes, this high
    #we could get undone by AP though.. also sun

    rg, azg = np.meshgrid(radar.range['data'], radar.azimuth['data'] )
    rg, eleg = np.meshgrid(radar.range['data'], radar.elevation['data'])
    x,y,z = radar_coords_to_cart(rg /1000.0, azg, eleg) #FLAG: need to fix

    points_above = np.where(z > toa)
    noise_floor_estimate = pseudo_power[points_above].mean()
    snr_dict = get_metadata(radar.fields.keys()[0])
    snr_dict['data'] = pseudo_power - noise_floor_estimate
    snr_dict['units'] = 'dB'
    snr_dict['standard_name'] = 'SNR'
    snr_dict['long_name'] = 'Signal to Noise Ratio'
    return snr_dict
