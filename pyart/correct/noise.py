"""
pyart.correct.noise
===================

Corrects polarimetric variables for noise

.. autosummary::
    :toctree: generated/

    correct_noise_rhohv

"""

import numpy as np

from ..config import get_metadata, get_field_name, get_fillvalue


def correct_noise_rhohv(radar, urhohv_field=None, snr_field=None,
                        zdr_field=None, nh_field=None, nv_field=None,
                        rhohv_field=None):
    """
    Corrects RhoHV for noise

    Parameters
    ----------
    radar : Radar
        radar object

    urhohv_field : str
        name of the RhoHV uncorrected for noise field

    snr_field, zdr_field, nh_field, nv_field: str
        names of the SNR, ZDR, horizontal channel noise in dBZ and vertical
        channel noise in dBZ used to correct RhoHV

    rhohv_field: str
        name of the rhohv field to output

    Returns
    -------
    rhohv : dictionary of dictionaries
        noise corrected RhoHV field

    """
    # parse the field parameters
    if urhohv_field is None:
        rhohv_field = get_field_name('uncorrected_cross_correlation_ratio')
    if snr_field is None:
        snr_field = get_field_name('signal_to_noise_ratio')
    if zdr_field is None:
        zdr_field = get_field_name('differential_reflectivity')
    if nh_field is None:
        nh_field = get_field_name('noisedBZ_hh')
    if nv_field is None:
        nv_field = get_field_name('noisedBZ_vv')
    if rhohv_field is None:
        rhohv_field = get_field_name('cross_correlation_ratio')

    # extract fields from radar
    if urhohv_field in radar.fields:
        urhohv = radar.fields[urhohv_field]['data']
    else:
        raise KeyError('Field not available: ' + urhohv_field)
    if snr_field in radar.fields:
        snrdB_h = radar.fields[snr_field]['data']
    else:
        raise KeyError('Field not available: ' + snr_field)
    if zdr_field in radar.fields:
        zdrdB = radar.fields[zdr_field]['data']
    else:
        raise KeyError('Field not available: ' + zdr_field)
    if nh_field in radar.fields:
        nh = radar.fields[nh_field]['data']
    else:
        raise KeyError('Field not available: ' + nh_field)
    if nv_field in radar.fields:
        nv = radar.fields[nv_field]['data']
    else:
        raise KeyError('Field not available: ' + nv_field)

    snr_h = np.ma.power(10., 0.1*snrdB_h)
    zdr = np.ma.power(10., 0.1*zdrdB)
    alpha = np.ma.power(10., 0.1*(nh-nv))

    mask = np.ma.getmaskarray(urhohv)
    fill_value = urhohv.get_fill_value()

    rhohv_data = np.ma.masked_where(
        mask, urhohv*np.ma.sqrt((1.+1./snr_h)*(1.+zdr/(alpha*snr_h))))
    rhohv_data.set_fill_value(fill_value)
    rhohv_data.data[mask.nonzero()] = fill_value
    is_above1 = rhohv_data > 1.
    rhohv_data[is_above1.nonzero()] = 1.

    rhohv = get_metadata(rhohv_field)
    rhohv['data'] = rhohv_data

    return rhohv
