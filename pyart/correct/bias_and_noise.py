"""
Corrects polarimetric variables for noise

"""

import numpy as np

from ..config import get_field_name, get_metadata


def correct_noise_rhohv(
    radar,
    urhohv_field=None,
    snr_field=None,
    zdr_field=None,
    nh_field=None,
    nv_field=None,
    rhohv_field=None,
):
    """
    Corrects RhoHV for noise according to eq. 6 in Gourley et al. 2006.
    This correction should only be performed if noise has not been subtracted
    from the signal during the moments computation.

    Parameters
    ----------
    radar : Radar
        Radar object.
    urhohv_field : str, optional
        Name of the RhoHV uncorrected for noise field.
    snr_field, zdr_field, nh_field, nv_field : str, optional
        Names of the SNR, ZDR, horizontal channel noise in dBZ and vertical
        channel noise in dBZ used to correct RhoHV.
    rhohv_field : str, optional
        Name of the rhohv field to output.

    Returns
    -------
    rhohv : dict
        Noise corrected RhoHV field.

    References
    ----------
    Gourley et al. Data Quality of the Meteo-France C-Band Polarimetric
    Radar, JAOT, 23, 1340-1356

    """
    # parse the field parameters
    if urhohv_field is None:
        urhohv_field = get_field_name("uncorrected_cross_correlation_ratio")
    if snr_field is None:
        snr_field = get_field_name("signal_to_noise_ratio")
    if zdr_field is None:
        zdr_field = get_field_name("differential_reflectivity")
    if nh_field is None:
        nh_field = get_field_name("noisedBZ_hh")
    if nv_field is None:
        nv_field = get_field_name("noisedBZ_vv")
    if rhohv_field is None:
        rhohv_field = get_field_name("cross_correlation_ratio")

    # extract fields from radar
    if urhohv_field in radar.fields:
        urhohv = radar.fields[urhohv_field]["data"]
    else:
        raise KeyError("Field not available: " + urhohv_field)
    if snr_field in radar.fields:
        snrdB_h = radar.fields[snr_field]["data"]
    else:
        raise KeyError("Field not available: " + snr_field)
    if zdr_field in radar.fields:
        zdrdB = radar.fields[zdr_field]["data"]
    else:
        raise KeyError("Field not available: " + zdr_field)
    if nh_field in radar.fields:
        nh = radar.fields[nh_field]["data"]
    else:
        raise KeyError("Field not available: " + nh_field)
    if nv_field in radar.fields:
        nv = radar.fields[nv_field]["data"]
    else:
        raise KeyError("Field not available: " + nv_field)

    snr_h = np.ma.power(10.0, 0.1 * snrdB_h)
    zdr = np.ma.power(10.0, 0.1 * zdrdB)
    alpha = np.ma.power(10.0, 0.1 * (nh - nv))

    rhohv_data = urhohv * np.ma.sqrt(
        (1.0 + 1.0 / snr_h) * (1.0 + zdr / (alpha * snr_h))
    )
    rhohv_data[rhohv_data > 1.0] = 1.0

    rhohv = get_metadata(rhohv_field)
    rhohv["data"] = rhohv_data

    return rhohv


def correct_bias(radar, bias=0.0, field_name=None):
    """
    Corrects a radar data bias. If field name is none the correction is
    applied to horizontal reflectivity by default.

    Parameters
    ----------
    radar : Radar
        Radar object.
    bias : float, optional
        The bias magnitude.
    field_name: str, optional
        Names of the field to be corrected.

    Returns
    -------
    corrected_field : dict
        The corrected field

    """
    # parse the field parameters
    if field_name is None:
        field_name = get_field_name("reflectivity")

    # extract fields from radar
    if field_name in radar.fields:
        field_data = radar.fields[field_name]["data"]
    else:
        raise KeyError("Field not available: " + field_name)

    corr_field_data = field_data - bias

    if field_name.startswith("corrected_"):
        corr_field_name = field_name
    else:
        corr_field_name = "corrected_" + field_name

    corr_field = get_metadata(corr_field_name)
    corr_field["data"] = corr_field_data

    return corr_field

def calc_zdr_offset(radar, gatefilter=None, height_range=None, zdr_var=None):
    """
    Function for calculating the ZDR bias from a VPT scan.

    Parameters
    ----------
    radar : PyART radar object
        Radar object with radar data
    gatefilter: PyART GateFilter
        Gatefilter for filtering out data for calculating ZDR bias
    height_range: tuple 
        The minimum and maximum elevations for the scan. 
    zdr_var: string or None
        The name of the ZDR variable. Set to None to have PyART try to determine this automatically.
    
    Returns
    -------
    obj : dict
        The mean vertical profiles of each radar moment are extracted along with the ZDR bias.

    """
    if height_range is None:
        height_range = (0, 100000.)

    if zdr_var is None:
        zdr_var = get_field_name("differential_reflectivity")

    height_mask = np.logical_and(radar.range["data"] >= height_range[0],
                                 radar.range["data"] <= height_range[1])
    
    mask = gatefilter.gate_included
    Zdr = radar.fields[zdr_var]["data"]
    if isinstance(Zdr, np.ma.MaskedArray):
        Zdr = Zdr.filled(np.nan)
    Zdr = np.where(mask == True, Zdr, np.nan)
    bias = np.nanmean(Zdr[:, height_mask])
    results = {
        'bias': bias,
        'profile_zdr': np.nanmean(Zdr[:, height_mask], axis=0),
        'range': radar.range["data"][height_mask],
    }
    for k in radar.fields.keys():
        if k != "range":
            field_data = radar.fields[k]["data"].astype(float)
            if isinstance(field_data, np.ma.MaskedArray):
                field_data = field_data.filled(np.nan)
            
            field_data = np.where(mask == True, field_data, np.nan)
            results['profile_' + k] = np.nanmean(field_data[:, :], axis=0)
    return results