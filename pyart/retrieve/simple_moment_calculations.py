"""
pyart.retrieve.simple_moment_calculations
=========================================

Simple moment calculations.

.. autosummary::
    :toctree: generated/

    calculate_snr_from_reflectivity
    compute_noisedBZ
    compute_snr
    compute_l
    compute_cdr

"""

import numpy as np

from ..config import get_metadata, get_field_name, get_fillvalue
from ..core.transforms import antenna_to_cartesian


def calculate_snr_from_reflectivity(
        radar, refl_field=None, snr_field=None, toa=25000.):
    """
    Calculate the signal to noise ratio, in dB, from the reflectivity field.

    Parameters
    ----------
    radar : Radar
        Radar object from which to retrieve reflectivity field.
    refl_field : str, optional
        Name of field in radar which contains the reflectivity.
        None will use the default field name in the Py-ART configuration file.
    snr_field : str, optional
        Name to use for snr metadata.
        None will use the default field name in the Py-ART configuration file.
    toa : float, optional
        Height above which to take noise floor measurements, in meters.

    Returns
    -------
    snr : field dictionary
        Field dictionary containing the signal to noise ratio.

    """
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if snr_field is None:
        snr_field = get_field_name('signal_to_noise_ratio')

    range_grid = np.meshgrid(radar.range['data'],
                             np.ma.ones(radar.time['data'].shape))[0] + 1.0

    # remove range scale.. This is basically the radar constant scaled dBm
    pseudo_power = (radar.fields[refl_field]['data'] -
                    20.0*np.log10(range_grid / 1000.0))

    # Noise floor estimate
    # 25km.. should be no scatterers, not even planes, this high
    # we could get undone by AP though.. also sun
    rg, azg = np.meshgrid(radar.range['data'], radar.azimuth['data'])
    rg, eleg = np.meshgrid(radar.range['data'], radar.elevation['data'])
    x, y, z = antenna_to_cartesian(rg / 1000.0, azg, eleg)  # XXX: need to fix

    points_above = np.where(z > toa)
    noise_floor_estimate = pseudo_power[points_above].mean()

    snr_dict = get_metadata(snr_field)
    snr_dict['data'] = pseudo_power - noise_floor_estimate
    return snr_dict


def compute_noisedBZ(nrays, noisedBZ_val, range, ref_dist,
                     noise_field=None):
    """
    Computes noise in dBZ from reference noise value.

    Parameters
    ----------
    nrays: int
        number of rays in the reflectivity field

    noisedBZ_val: float
        Estimated noise value in dBZ at reference distance

    range: np array of floats
        range vector in m

    ref_dist: float
        reference distance in Km

    noise_field: str
        name of the noise field to use

    Returns
    -------
    noisedBZ : dict
        the noise field

    """
    # parse the field parameters
    if noise_field is None:
        noise_field = get_field_name('noisedBZ_hh')

    noisedBZ_vec = noisedBZ_val+20.*np.ma.log10(1e-3*range/ref_dist)

    noisedBZ = get_metadata(noise_field)
    noisedBZ['data'] = np.tile(noisedBZ_vec, (nrays, 1))
#    noisedBZ['data'].set_fill_value(get_fillvalue())

    return noisedBZ


def compute_snr(radar, refl_field=None, noise_field=None, snr_field=None):
    """
    Computes SNR from a reflectivity field and the noise in dBZ.

    Parameters
    ----------
    radar : Radar
        radar object

    refl_field, noise_field : str
        name of the reflectivity and noise field used for the calculations

    snr_field : str
        name of the SNR field

    Returns
    -------
    snr : dict
        the SNR field

    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if noise_field is None:
        noise_field = get_field_name('noisedBZ_hh')
    if snr_field is None:
        snr_field = get_field_name('signal_to_noise_ratio')

    # extract fields from radar
    if refl_field in radar.fields:
        refl = radar.fields[refl_field]['data']
    else:
        raise KeyError('Field not available: ' + refl_field)
    if noise_field in radar.fields:
        noisedBZ = radar.fields[noise_field]['data']
    else:
        raise KeyError('Field not available: ' + noise_field)

    snr_data = refl-noisedBZ
#    snr_data.set_fill_value(get_fillvalue())

    snr = get_metadata(snr_field)
    snr['data'] = snr_data

    return snr


def compute_l(radar, rhohv_field=None, l_field=None):
    """
    Computes Rhohv in logarithmic scale according to L=-log10(1-RhoHV)

    Parameters
    ----------
    radar : Radar
        radar object

    rhohv_field : str
        name of the RhoHV field used for the calculation

    l_field : str
        name of the L field

    Returns
    -------
    l : dict
        L field

    """
    # parse the field parameters
    if rhohv_field is None:
        rhohv_field = get_field_name('cross_correlation_ratio')
    if l_field is None:
        l_field = get_field_name('logarithmic_cross_correlation_ratio')

    # extract rhohv field from radar
    if rhohv_field in radar.fields:
        rhohv = radar.fields[rhohv_field]['data']
    else:
        raise KeyError('Field not available: ' + rhohv_field)

#    mask = np.ma.getmaskarray(rhohv)

    rhohv[rhohv >= 1.] = 0.9999
    l_data = -np.ma.log10(1.-rhohv)
#    l_data.data[mask] = np.ma.masked
#    l_data.set_fill_value(get_fillvalue())

    l = get_metadata(l_field)
    l['data'] = l_data

    return l


def compute_cdr(radar, rhohv_field=None, zdr_field=None, cdr_field=None):
    """
    Computes the Circular Depolarization Ratio

    Parameters
    ----------
    radar : Radar
        radar object

    rhohv_field, zdr_field : str
        name of the input RhoHV and ZDR fields

    cdr_field : str
        name of the CDR field

    Returns
    -------
    cdr : dict
        CDR field

    """
    # parse the field parameters
    if rhohv_field is None:
        rhohv_field = get_field_name('cross_correlation_ratio')
    if zdr_field is None:
        zdr_field = get_field_name('differential_reflectivity')
    if cdr_field is None:
        cdr_field = get_field_name('circular_depolarization_ratio')

    # extract fields from radar
    if rhohv_field in radar.fields:
        rhohv = radar.fields[rhohv_field]['data']
    else:
        raise KeyError('Field not available: ' + rhohv_field)
    if zdr_field in radar.fields:
        zdrdB = radar.fields[zdr_field]['data']
    else:
        raise KeyError('Field not available: ' + zdr_field)

    zdr = np.ma.power(10., 0.1*zdrdB)

    cdr_data = (
        10.*np.ma.log10(
            (1.+1./zdr-2.*rhohv*np.ma.sqrt(1./zdr)) /
            (1.+1./zdr+2.*rhohv*np.ma.sqrt(1./zdr))))

    cdr = get_metadata(cdr_field)
    cdr['data'] = cdr_data

    return cdr
