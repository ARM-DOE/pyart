"""
pyart.retrieve.qpe
=========================================

Functions for rainfall rate estimation

.. autosummary::
    :toctree: generated/

    get_lowest_elevation
    rr_zpoly
    rr_z
    rr_kdp
    rr_a
    rr_zkdp
    rr_za
    rr_hydro

"""


import os
import sys
from time import time
from warnings import warn

import numpy as np
import netCDF4
import matplotlib as mpl
import pylab as pl

from ..util import met
from ..config import get_metadata, get_field_name, get_fillvalue
# import ballsy_masked as ballsy


def get_lowest_elevation(mmcg_ncf, req_moments, corrections):
    pass
    # for each of the moments in req_moments fetch the lowest valid return


def rr_zpoly(radar, refl_field=None, rr_field=None):
    """
    Estimates rainfall rate from reflectivity using a polynomial Z-R relation
    developed at McGill University

    Parameters
    ----------
    radar : Radar
        Radar object

    refl_field : str
        name of the reflectivity field to use

    rr_field : str
        name of the rainfall rate field

    Returns
    -------
    rain : dict
        Field dictionary containing the rainfall rate.

    """

    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if rr_field is None:
        rr_field = get_field_name('radar_estimated_rain_rate')

    if refl_field in radar.fields:
        refl = radar.fields[refl_field]['data']
    else:
        raise KeyError('Field not available: ' + refl_field)

    refl2 = refl*refl
    refl3 = refl*refl2
    refl4 = refl*refl3

    rr_data = np.ma.power(
        10., -2.3+0.17*refl-5.1e-3*refl2+9.8e-5*refl3-6e-7*refl4)

    rain = get_metadata(rr_field)
    rain['data'] = rr_data

    return rain


def rr_z(radar, alpha=0.0376, beta=0.6112, refl_field=None, rr_field=None):
    """
    Estimates rainfall rate from reflectivity using a power law

    Parameters
    ----------
    radar : Radar
        Radar object

    alpha,beta : floats
        factor (alpha) and exponent (beta) of the power law

    refl_field : str
        name of the reflectivity field to use

    rr_field : str
        name of the rainfall rate field

    Returns
    -------
    rain : dict
        Field dictionary containing the rainfall rate.

    """

    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if rr_field is None:
        rr_field = get_field_name('radar_estimated_rain_rate')

    if refl_field in radar.fields:
        refl = radar.fields[refl_field]['data']
    else:
        raise KeyError('Field not available: ' + refl_field)

    rr_data = alpha*np.ma.power(np.ma.power(10., 0.1*refl), beta)

    rain = get_metadata(rr_field)
    rain['data'] = rr_data

    return rain


def rr_kdp(radar, alpha=None, beta=None, kdp_field=None, rr_field=None):
    """
    Estimates rainfall rate from kdp using alpha power law

    Parameters
    ----------
    radar : Radar
        Radar object

    alpha,beta : floats
        Optional. factor (alpha) and exponent (beta) of the power law.
        If not set the factors are going to be determined according
        to the radar frequency

    kdp_field : str
        name of the specific differential phase field to use

    rr_field : str
        name of the rainfall rate field

    Returns
    -------
    rain : dict
        Field dictionary containing the rainfall rate.

    """
    # select the coefficients as alpha function of frequency band
    if alpha is None or beta is None:
        # assign coefficients according to radar frequency
        if 'frequency' in radar.instrument_parameters:
            freq = radar.instrument_parameters['frequency']['data'][0]
            # S band: Beard and Chuang coefficients
            if freq >= 2e9 and freq < 4e9:
                freq_band = 'S'
                alpha = 50.7
                beta = 0.85
            # C band: Beard and Chuang coefficients
            elif freq >= 4e9 and freq < 8e9:
                freq_band = 'C'
                alpha = 29.7
                beta = 0.85
            # X band: Brandes coefficients
            elif freq >= 8e9 and freq <= 12e9:
                freq_band = 'X'
                alpha = 15.810
                beta = 0.7992
            else:
                if freq < 2e9:
                    freq_band = 'S'
                    alpha = 50.7
                    beta = 0.85
                else:
                    freq_band = 'X'
                    alpha = 15.810
                    beta = 0.7992
                warn('Radar frequency out of range. ' +
                     'Coefficients only applied to S, C or X band. ' +
                     freq_band + ' band coefficients will be used')
        else:
            freq_band = 'C'
            alpha = 29.7
            beta = 0.85
            warn('Radar frequency unknown. ' +
                 'Default coefficients for C band will be applied')

    # parse the field parameters
    if kdp_field is None:
        kdp_field = get_field_name('specific_differential_phase')
    if rr_field is None:
        rr_field = get_field_name('radar_estimated_rain_rate')

    if kdp_field in radar.fields:
        kdp = radar.fields[kdp_field]['data']
    else:
        raise KeyError('Field not available: ' + kdp_field)

    kdp[kdp < 0] = 0.
    rr_data = alpha*np.ma.power(kdp, beta)

    rain = get_metadata(rr_field)
    rain['data'] = rr_data

    return rain


def rr_a(radar, alpha=None, beta=None, a_field=None, rr_field=None):
    """
    Estimates rainfall rate from specific attenuation using alpha power law

    Parameters
    ----------
    radar : Radar
        Radar object

    alpha,beta : floats
        Optional. factor (alpha) and exponent (beta) of the power law.
        If not set the factors are going to be determined according
        to the radar frequency

    a_field : str
        name of the specific attenuation field to use

    rr_field : str
        name of the rainfall rate field

    Returns
    -------
    rain : dict
        Field dictionary containing the rainfall rate.

    References
    ----------
    Diederich M., Ryzhkov A., Simmer C., Zhang P. and Tromel S., 2015: Use of
    Specific Attenuation for Rainfall Measurement at X-Band Radar Wavelenghts.
    Part I: Radar Calibration and Partial Beam Blockage Estimation. Journal of
    Hydrometeorology, 16, 487-502.

    Ryzhkov A., Diederich M., Zhang P. and Simmer C., 2014: Potential
    Utilization of Specific Attenuation for Rainfall Estimation, Mitigation of
    Partial Beam Blockage, and Radar Networking. Journal of Atmospheric and
    Oceanic Technology, 31, 599-619.

    """
    # select the coefficients as alpha function of frequency band
    if alpha is None or beta is None:
        # assign coefficients according to radar frequency
        if 'frequency' in radar.instrument_parameters:
            freq = radar.instrument_parameters['frequency']['data'][0]
            # S band: at 10°C according to tables from
            # Ryzhkov et al. 2014
            if freq >= 2e9 and freq < 4e9:
                freq_band = 'S'
                alpha = 3100.
                beta = 1.03
            # C band: at 10°C according to tables from
            # Diederich et al. 2015
            elif freq >= 4e9 and freq < 8e9:
                freq_band = 'C'
                alpha = 250.
                beta = 0.91
            # X band: at 10°C according to tables from
            # Diederich et al. 2015
            elif freq >= 8e9 and freq <= 12e9:
                freq_band = 'X'
                alpha = 45.5
                beta = 0.83
            else:
                if freq < 2e9:
                    freq_band = 'S'
                    alpha = 3100.
                    beta = 1.03
                else:
                    freq_band = 'X'
                    alpha = 45.5
                    beta = 0.83
                warn('Radar frequency out of range. ' +
                     'Coefficients only applied to S, C or X band. ' +
                      freq_band + ' band coefficients will be used')
        else:
            freq_band = 'C'
            alpha = 250.
            beta = 0.91
            warn('Radar frequency unknown. ' +
                 'Default coefficients for C band will be applied')

    # parse the field parameters
    if a_field is None:
        a_field = get_field_name('specific_attenuation')
    if rr_field is None:
        rr_field = get_field_name('radar_estimated_rain_rate')

    if a_field in radar.fields:
        att = radar.fields[a_field]['data']
    else:
        raise KeyError('Field not available: ' + a_field)

    rr_data = alpha*np.ma.power(att, beta)

    rain = get_metadata(rr_field)
    rain['data'] = rr_data

    return rain


def rr_zkdp(radar, alphaz=0.0376, betaz=0.6112, alphakdp=None, betakdp=None,
            refl_field=None, kdp_field=None, rr_field=None,
            master_field=None, thresh=None, thresh_max=True):
    """
    Estimates rainfall rate from a blending of power law r-kdp and r-z
    relations.

    Parameters
    ----------
    radar : Radar
        Radar object

    alphaz,betaz : floats
        factor (alpha) and exponent (beta) of the z-r power law.

    alphakdp, betakdp : floats
        Optional. factor (alpha) and exponent (beta) of the kdp-r power law.
        If not set the factors are going to be determined according
        to the radar frequency

    refl_field : str
        name of the reflectivity field to use

    kdp_field : str
        name of the specific differential phase field to use

    rr_field : str
        name of the rainfall rate field

    master_field : str
        name of the field that is going to act as master. Has to be
        either refl_field or kdp_field. Default is refl_field
    thresh : float
        value of the threshold that determines when to use the slave
        field.
    thresh_max : Boolean
        If true the master field is used up to the thresh value maximum.
        Otherwise the master field is not used below thresh value.

    Returns
    -------
    rain_master : dict
        Field dictionary containing the rainfall rate.

    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if kdp_field is None:
        kdp_field = get_field_name('specific_differential_phase')
    if rr_field is None:
        rr_field = get_field_name('radar_estimated_rain_rate')

    rain_z = rr_z(radar, alpha=alphaz, beta=betaz, refl_field=refl_field,
                  rr_field=rr_field)
    rain_kdp = rr_kdp(radar, alpha=alphakdp, beta=betakdp,
                      kdp_field=kdp_field, rr_field=rr_field)

    if master_field == refl_field:
        slave_field = kdp_field
        rain_master = rain_z
        rain_slave = rain_kdp
    elif master_field == kdp_field:
        slave_field = refl_field
        rain_master = rain_kdp
        rain_slave = rain_z
    elif master_field is None:
        master_field = refl_field
        slave_field = kdp_field
        rain_master = rain_z
        rain_slave = rain_kdp
    else:
        master_field = refl_field
        slave_field = kdp_field
        rain_master = rain_z
        rain_slave = rain_kdp
        thresh = 40.
        thresh_max = True
        warn('Unknown master field. Using ' + refl_field +
             ' with threshold '+str(thresh))

    if thresh_max:
        is_slave = rain_master['data'] > thresh
    else:
        is_slave = rain_master['data'] < thresh
    rain_master['data'][is_slave] = (
        rain_slave['data'][is_slave()])

    return rain_master


def rr_za(radar, alphaz=0.0376, betaz=0.6112, alphaa=None, betaa=None,
          refl_field=None, a_field=None, rr_field=None,
          master_field=None, thresh=None, thresh_max=False):
    """
    Estimates rainfall rate from a blending of power law r-alpha and r-z
    relations.

    Parameters
    ----------
    radar : Radar
        Radar object

    alphaz,betaz : floats
        factor (alpha) and exponent (beta) of the z-r power law.

    alphaa,betaa : floats
        Optional. factor (alpha) and exponent (beta) of the a-r power law.
        If not set the factors are going to be determined according
        to the radar frequency

    refl_field : str
        name of the reflectivity field to use

    a_field : str
        name of the specific attenuation field to use

    rr_field : str
        name of the rainfall rate field

    master_field : str
        name of the field that is going to act as master. Has to be
        either refl_field or kdp_field. Default is refl_field
    thresh : float
        value of the threshold that determines when to use the slave
        field.
    thresh_max : Boolean
        If true the master field is used up to the thresh value maximum.
        Otherwise the master field is not used below thresh value.

    Returns
    -------
    rain_master : dict
        Field dictionary containing the rainfall rate.

    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if a_field is None:
        a_field = get_field_name('specific_attenuation')
    if rr_field is None:
        rr_field = get_field_name('radar_estimated_rain_rate')

    rain_z = rr_z(radar, alpha=alphaz, beta=betaz, refl_field=refl_field,
                  rr_field=rr_field)
    rain_a = rr_a(radar, alpha=alphaa, beta=betaa, a_field=a_field,
                  rr_field=rr_field)

    if master_field == refl_field:
        slave_field = a_field
        rain_master = rain_z
        rain_slave = rain_a
    elif master_field == a_field:
        slave_field = refl_field
        rain_master = rain_a
        rain_slave = rain_z
    elif master_field is None:
        master_field = a_field
        slave_field = refl_field
        rain_master = rain_a
        rain_slave = rain_z
    else:
        master_field = a_field
        slave_field = refl_field
        rain_master = rain_a
        rain_slave = rain_z
        thresh = 0.04
        thresh_max = False
        warn('Unknown master field. Using ' + a_field + ' with threshold ' +
             str(thresh))

    if thresh_max:
        is_slave = rain_master['data'] > thresh
    else:
        is_slave = rain_master['data'] < thresh

    rain_master['data'][is_slave] = (
        rain_slave['data'][is_slave])

    return rain_master


def rr_hydro(radar, alphazr=0.0376, betazr=0.6112, alphazs=0.1, betazs=0.5,
             alphaa=None, betaa=None, mp_factor=0.6, refl_field=None,
             a_field=None, hydro_field=None, rr_field=None, master_field=None,
             thresh=None, thresh_max=False):
    """
    Estimates rainfall rate using different relations between R and the
    polarimetric variables depending on the hydrometeor type

    Parameters
    ----------
    radar : Radar
        Radar object

    alphazr,betazr : floats
        factor (alpha) and exponent (beta) of the z-r power law for rain.

    alphazs,betazs : floats
        factor (alpha) and exponent (beta) of the z-s power law for snow.

    alphaa,betaa : floats
        Optional. factor (alpha) and exponent (beta) of the a-r power law.
        If not set the factors are going to be determined according
        to the radar frequency

    mp_factor : float
        factor applied to z-r relation in the melting layer

    refl_field : str
        name of the reflectivity field to use

    a_field : str
        name of the specific attenuation field to use

    hydro_field : str
        name of the hydrometeor classification field to use

    rr_field : str
        name of the rainfall rate field

    master_field : str
        name of the field that is going to act as master. Has to be
        either refl_field or kdp_field. Default is refl_field
    thresh : float
        value of the threshold that determines when to use the slave
        field.
    thresh_max : Boolean
        If true the master field is used up to the thresh value maximum.
        Otherwise the master field is not used below thresh value.

    Returns
    -------
    rain : dict
        Field dictionary containing the rainfall rate.

    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if a_field is None:
        a_field = get_field_name('specific_attenuation')
    if hydro_field is None:
        hydro_field = get_field_name('radar_echo_classification')
    if rr_field is None:
        rr_field = get_field_name('radar_estimated_rain_rate')

    # extract fields and parameters from radar
    if hydro_field in radar.fields:
        hydroclass = radar.fields[hydro_field]['data']
    else:
        raise KeyError('Field not available: ' + hydro_field)

    # get the location of each hydrometeor class
    is_ds = hydroclass == 1
    is_cr = hydroclass == 2
    is_lr = hydroclass == 3
    is_gr = hydroclass == 4
    is_rn = hydroclass == 5
    is_vi = hydroclass == 6
    is_ws = hydroclass == 7
    is_mh = hydroclass == 8
    is_ih = hydroclass == 9

    # compute z-r (in rain) z-r in snow and z-a relations
    rain_z = rr_z(radar, alpha=alphazr, beta=betazr,
                  refl_field=refl_field, rr_field=rr_field)
    snow_z = rr_z(radar, alpha=alphazs, beta=betazs,
                  refl_field=refl_field, rr_field=rr_field)
    rain_a = rr_a(radar, alpha=alphaa, beta=betaa,
                  a_field=a_field, rr_field=rr_field)

    # initialize rainfall rate field
    rr_data = np.ma.zeros(hydroclass.shape, dtype='float32')
    rr_data[:] = np.ma.masked
    rr_data.set_fill_value(get_fillvalue())

    # apply the relations for each type
    # solid phase
    rr_data[is_ds] = snow_z['data'][is_ds]
    rr_data[is_cr] = snow_z['data'][is_cr]
    rr_data[is_vi] = snow_z['data'][is_vi]
    rr_data[is_gr] = snow_z['data'][is_gr]
    rr_data[is_ih] = snow_z['data'][is_ih]

    # rain
    if master_field == refl_field:
        slave_field = a_field
        rain_master = rain_z
        rain_slave = rain_a
    elif master_field == a_field:
        slave_field = refl_field
        rain_master = rain_a
        rain_slave = rain_z
    elif master_field is None:
        master_field = a_field
        slave_field = refl_field
        rain_master = rain_a
        rain_slave = rain_z
    else:
        master_field = a_field
        slave_field = refl_field
        rain_master = rain_a
        rain_slave = rain_z
        thresh = 0.04
        thresh_max = False
        warn('Unknown master field. Using ' + a_field + ' with threshold ' +
             str(thresh))

    if thresh_max:
        is_slave = rain_master['data'] > thresh
    else:
        is_slave = rain_master['data'] < thresh

    rain_master['data'][is_slave] = rain_slave['data'][is_slave]

    rr_data[is_lr] = rain_master['data'][is_lr]
    rr_data[is_rn] = rain_master['data'][is_rn]

    # mixed phase
    rr_data[is_ws] = mp_factor*rain_z['data'][is_ws]
    rr_data[is_mh] = mp_factor*rain_z['data'][is_mh]

    rain = get_metadata(rr_field)
    rain['data'] = rr_data

    return rain
