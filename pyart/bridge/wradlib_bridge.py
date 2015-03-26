"""
pyart.bridge.wradlib
=================

Pyart methods linking to wradlib functions http://wradlib.bitbucket.org/

.. autosummary::
    :toctree: generated/

     texture_of_complex_phase

"""
import wradlib
import numpy as np
from ..config import get_metadata

def texture_of_complex_phase(radar, phidp_key = 'differential_phase'):
    """
    Calculate the texture of the real part of the complex differential
    phase field
    ----------
    radar : Radar
        Radar object to use.

    Returns
    -------
    texture : dict
        Field dictionary containing the texture of the real part
        of the complex differential
    phase field.

    Other Parameters
    ----------------
    phidp_key : string
        key for the differential phase field


    References
    ----------
    Gourley, J. J., P. Tabary, and J. Parent du Chatelet, 2007
    sing Polarimetric Radar Observations.
    """
    #Grab the phase data
    phidp_array = radar.fields[phidp_key]['data']

    #convert to complex number
    complex_phase = np.exp(1j*(phidp_array*np.pi/180.0))

    #calculate texture using wradlib
    w_texture_complex = wradlib.dp.texture((np.real(complex_phase) +1.0)*180)
    texture_field = get_metadata(phidp_key)
    texture_field['data'] = w_texture_complex
    texture_field['long_name'] = 'Texture of Differential phase (PhiDP)'
    texture_field['standard_name'] = 'differential_phase_hv_texture'
    texture_field['valid_min'] = 0.0
    texture_field['valid_max'] = 400.0
    return texture_field
