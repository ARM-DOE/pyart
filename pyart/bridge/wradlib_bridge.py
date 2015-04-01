"""
pyart.bridge.wradlib
====================

Py-ART methods linking to wradlib functions, http://wradlib.bitbucket.org/

.. autosummary::
    :toctree: generated/

     texture_of_complex_phase

"""


import wradlib
import numpy as np

from ..config import get_metadata, get_field_name


def texture_of_complex_phase(radar, phidp_field=None,
                             phidp_texture_field=None):
    """
    Calculate the texture of the differential phase field.

    Calculate the texture of the real part of the complex differential
    phase field

    Parameters
    ----------
    radar : Radar
        Radar object from which to .
    phidp_field : str, optional
        Name of field in radar which contains the differential phase shift.
        None will use the default field name in the Py-ART configuration file.
    phidp_texture_field : str, optional
        Name to use for the differential phase texture field metadata.
        None will use the default field name in the Py-ART configuration file.

    Returns
    -------
    texture_field : dict
        Field dictionary containing the texture of the real part
        of the complex differential phase.

    References
    ----------
    Gourley, J. J., P. Tabary, and J. Parent du Chatelet,
    A fuzzy logic algorithm for the separation of precipitating from
    nonprecipitating echoes using polarimetric radar observations,
    Journal of Atmospheric and Oceanic Technology 24 (8), 1439-1451

    """
    # parse field names
    if phidp_field is None:
        phidp_field = get_field_name('differential_phase')
    if phidp_texture_field is None:
        phidp_field = get_field_name('differential_phase')

    # Grab the phase data
    phidp = radar.fields[phidp_field]['data']

    # convert to complex number
    complex_phase = np.exp(1j*(phidp*np.pi/180.0))

    # calculate texture using wradlib
    w_texture_complex = wradlib.dp.texture(
        (np.real(complex_phase) + 1.0) * 180)

    texture_field = get_metadata(phidp_texture_field)
    texture_field['data'] = w_texture_complex
    return texture_field
