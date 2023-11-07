"""
Rescale reflectivity to precipitation rate

"""

import numpy as np

from pyart.core import Radar


def ZtoR(radar, ref_field="reflectivity", a=300, b=1.4, save_name='NWS_primary_prate'):
    """
    Convert reflectivity (dBZ) to precipitation rate (mm/hr)

    Author: Laura Tomkins

    Parameters
    ----------
    radar : Radar
        Radar object used.
    ref_field : str
        Reflectivity field name to use to look up reflectivity data. In the
        radar object. Default field name is 'reflectivity'. Units are expected
        to be dBZ.
    a : float
        a value (coefficient) in the Z-R relationship
    b: float
        b value (exponent) in the Z-R relationship

    Returns
    -------
    radar : Radar
        The radar object containing the precipitation rate field

    References
    ----------
    American Meteorological Society, 2022: "Z-R relation". Glossary of Meteorology,
    https://glossary.ametsoc.org/wiki/Z-r_relation

    """

    # get reflectivity data
    ref_data = radar.fields[ref_field]['data']
    ref_data = np.ma.masked_invalid(ref_data)

    # convert to linear reflectivity
    ref_linear = 10 ** (ref_data / 10)
    precip_rate = (ref_linear / a) ** (1 / b)

    # create dictionary
    prate_dict = {
        'data': precip_rate,
        'standard_name': save_name,
        'long_name': "{} rescaled from linear reflectivity".format(save_name),
        'units': 'mm/hr',
        'valid_min': 0,
        'valid_max': 10000
    }

    # add field to radar object
    radar.add_field(save_name, prate_dict, replace_existing=True)

    return radar
