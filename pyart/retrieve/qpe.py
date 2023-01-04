"""
Functions for rainfall rate estimation.

"""

from warnings import warn

import numpy as np

from ..config import get_field_name, get_fillvalue, get_metadata
from .echo_class import get_freq_band


def est_rain_rate_zpoly(radar, refl_field=None, rr_field=None):
    """
    Estimates rainfall rate from reflectivity using a polynomial Z-R relation
    developed at McGill University.

    Parameters
    ----------
    radar : Radar
        Radar object.
    refl_field : str, optional
        Name of the reflectivity field to use.
    rr_field : str, optional
        Name of the rainfall rate field.

    Returns
    -------
    rain : dict
        Field dictionary containing the rainfall rate.

    References
    ----------
    Doelling et al. Systematic variations of Z–R-relationships from drop size
    distributions measured in northern Germany during seven years. 1998. Atmos.
    Ocean. Technol, 21, 1545-1556.

    Joss et al. Operational Use of Radar for Precipitation Measurements
    in Switzerland. 1998. Vdf Hochschulverlag AG ETH Zurich: 134.

    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name("reflectivity")
    if rr_field is None:
        rr_field = get_field_name("radar_estimated_rain_rate")

    radar.check_field_exists(refl_field)
    refl = radar.fields[refl_field]["data"]

    refl2 = refl * refl
    refl3 = refl * refl2
    refl4 = refl * refl3

    rr_data = np.ma.power(
        10.0, -2.3 + 0.17 * refl - 5.1e-3 * refl2 + 9.8e-5 * refl3 - 6e-7 * refl4
    )

    rain = get_metadata(rr_field)
    rain["data"] = rr_data

    return rain


def est_rain_rate_z(radar, alpha=0.0376, beta=0.6112, refl_field=None, rr_field=None):
    """
    Estimates rainfall rate from reflectivity using a power law.

    Parameters
    ----------
    radar : Radar
        Radar object.
    alpha, beta : floats, optional
        Factor (alpha) and exponent (beta) of the power law.
    refl_field : str, optional
        Name of the reflectivity field to use.
    rr_field : str, optional
        Name of the rainfall rate field.

    Returns
    -------
    rain : dict
        Field dictionary containing the rainfall rate.

    Reference
    ---------
    Fabry, Frédéric. Radar Meterology. 2015. Ch 9. pg 148-165.
    https://doi.org/10.1017/CBO9781107707405

    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name("reflectivity")
    if rr_field is None:
        rr_field = get_field_name("radar_estimated_rain_rate")

    radar.check_field_exists(refl_field)
    refl = radar.fields[refl_field]["data"]

    rr_data = alpha * np.ma.power(np.ma.power(10.0, 0.1 * refl), beta)

    rain = get_metadata(rr_field)
    rain["data"] = rr_data

    return rain


def est_rain_rate_kdp(radar, alpha=None, beta=None, kdp_field=None, rr_field=None):
    """
    Estimates rainfall rate from kdp using alpha power law.

    Parameters
    ----------
    radar : Radar
        Radar object.
    alpha, beta : floats, optional
        Factor (alpha) and exponent (beta) of the power law. If not set the
        factors are going to be determined according to the radar frequency.
    kdp_field : str, optional
        Name of the specific differential phase field to use.
    rr_field : str, optional
        Name of the rainfall rate field.

    Returns
    -------
    rain : dict
        Field dictionary containing the rainfall rate.

    Reference
    ---------
    Figueras et al. Long-term monitoring of French polarimetric radar data
    quality and evaluation of several polarimetric quantitative precipitation
    estimators in ideal conditions for operational implementation at C-band.
    Quarterly Journal of the Royal Meteorological Society. 2012.
    https://doi.org/10.1002/qj.1934

    """
    # select the coefficients as alpha function of frequency band
    if alpha is None or beta is None:
        # assign coefficients according to radar frequency
        if "frequency" in radar.instrument_parameters:
            alpha, beta = _get_coeff_rkdp(
                radar.instrument_parameters["frequency"]["data"][0]
            )
        else:
            alpha, beta = _coeff_rkdp_table()["C"]
            warn(
                "Radar frequency unknown. "
                "Default coefficients for C band will be applied."
            )

    # parse the field parameters
    if kdp_field is None:
        kdp_field = get_field_name("specific_differential_phase")
    if rr_field is None:
        rr_field = get_field_name("radar_estimated_rain_rate")

    radar.check_field_exists(kdp_field)
    kdp = radar.fields[kdp_field]["data"]

    kdp[kdp < 0] = 0.0
    rr_data = alpha * np.ma.power(kdp, beta)

    rain = get_metadata(rr_field)
    rain["data"] = rr_data

    return rain


def est_rain_rate_a(radar, alpha=None, beta=None, a_field=None, rr_field=None):
    """
    Estimates rainfall rate from specific attenuation using alpha power law.

    Parameters
    ----------
    radar : Radar
        Radar object.
    alpha, beta : floats, optional
        Factor (alpha) and exponent (beta) of the power law. If not set the
        factors are going to be determined according to the radar frequency.
    a_field : str, optional
        Name of the specific attenuation field to use.
    rr_field : str, optional
        Name of the rainfall rate field.

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
        if "frequency" in radar.instrument_parameters:
            alpha, beta = _get_coeff_ra(
                radar.instrument_parameters["frequency"]["data"][0]
            )
        else:
            alpha, beta = _coeff_ra_table()["C"]
            warn(
                "Radar frequency unknown. "
                "Default coefficients for C band will be applied."
            )

    # parse the field parameters
    if a_field is None:
        a_field = get_field_name("specific_attenuation")
    if rr_field is None:
        rr_field = get_field_name("radar_estimated_rain_rate")

    radar.check_field_exists(a_field)
    att = radar.fields[a_field]["data"]

    rr_data = alpha * np.ma.power(att, beta)

    rain = get_metadata(rr_field)
    rain["data"] = rr_data

    return rain


def est_rain_rate_zkdp(
    radar,
    alphaz=0.0376,
    betaz=0.6112,
    alphakdp=None,
    betakdp=None,
    refl_field=None,
    kdp_field=None,
    rr_field=None,
    main_field=None,
    thresh=None,
    thresh_max=True,
):
    """
    Estimates rainfall rate from a blending of power law r-kdp and r-z
    relations.

    Parameters
    ----------
    radar : Radar
        Radar object.
    alphaz, betaz : floats, optional
        Factor (alpha) and exponent (beta) of the z-r power law.
    alphakdp, betakdp : floats, optional
        Factor (alpha) and exponent (beta) of the kdp-r power law.
        If not set the factors are going to be determined according
        to the radar frequency.
    refl_field : str, optional
        Name of the reflectivity field to use.
    kdp_field : str, optional
        Name of the specific differential phase field to use.
    rr_field : str, optional
        Name of the rainfall rate field.
    main_field : str, optional
        Name of the field that is going to act as main. Has to be
        either refl_field or kdp_field. Default is refl_field.
    thresh : float, optional
        Value of the threshold that determines when to use the secondary
        field.
    thresh_max : Bool, optional
        If true the main field is used up to the thresh value maximum.
        Otherwise the main field is not used below thresh value.

    Returns
    -------
    rain_main : dict
        Field dictionary containing the rainfall rate.

    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name("reflectivity")
    if kdp_field is None:
        kdp_field = get_field_name("specific_differential_phase")
    if rr_field is None:
        rr_field = get_field_name("radar_estimated_rain_rate")

    rain_z = est_rain_rate_z(
        radar, alpha=alphaz, beta=betaz, refl_field=refl_field, rr_field=rr_field
    )
    rain_kdp = est_rain_rate_kdp(
        radar, alpha=alphakdp, beta=betakdp, kdp_field=kdp_field, rr_field=rr_field
    )

    if main_field == refl_field:
        rain_main = rain_z
        rain_secondary = rain_kdp
    elif main_field == kdp_field:
        rain_main = rain_kdp
        rain_secondary = rain_z
    elif main_field is None:
        main_field = refl_field
        rain_main = rain_z
        rain_secondary = rain_kdp
    else:
        main_field = refl_field
        rain_main = rain_z
        rain_secondary = rain_kdp
        thresh = 40.0
        thresh_max = True
        warn(
            "Unknown main field. Using " + refl_field + " with threshold " + str(thresh)
        )

    if thresh_max:
        is_secondary = rain_main["data"] > thresh
    else:
        is_secondary = rain_main["data"] < thresh
    rain_main["data"][is_secondary] = rain_secondary["data"][is_secondary]

    return rain_main


def est_rain_rate_za(
    radar,
    alphaz=0.0376,
    betaz=0.6112,
    alphaa=None,
    betaa=None,
    refl_field=None,
    a_field=None,
    rr_field=None,
    main_field=None,
    thresh=None,
    thresh_max=False,
):
    """
    Estimates rainfall rate from a blending of power law r-alpha and r-z
    relations.

    Parameters
    ----------
    radar : Radar
        Radar object
    alphaz, betaz : floats, optional
        Factor (alpha) and exponent (beta) of the z-r power law.
    alphaa,betaa : floats, optional
        Factor (alpha) and exponent (beta) of the a-r power law. If not set
        the factors are going to be determined according to the radar frequency.
    refl_field : str, optional
        Name of the reflectivity field to use.
    a_field : str, optional
        Name of the specific attenuation field to use.
    rr_field : str, optional
        Name of the rainfall rate field.
    main_field : str, optional
        Name of the field that is going to act as main. Has to be
        either refl_field or kdp_field. Default is refl_field.
    thresh : float, optional
        Value of the threshold that determines when to use the secondary
        field.
    thresh_max : Bool, optional
        If true the main field is used up to the thresh value maximum.
        Otherwise the main field is not used below thresh value.

    Returns
    -------
    rain_main : dict
        Field dictionary containing the rainfall rate.

    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name("reflectivity")
    if a_field is None:
        a_field = get_field_name("specific_attenuation")
    if rr_field is None:
        rr_field = get_field_name("radar_estimated_rain_rate")

    rain_z = est_rain_rate_z(
        radar, alpha=alphaz, beta=betaz, refl_field=refl_field, rr_field=rr_field
    )
    rain_a = est_rain_rate_a(
        radar, alpha=alphaa, beta=betaa, a_field=a_field, rr_field=rr_field
    )

    if main_field == refl_field:
        rain_main = rain_z
        rain_secondary = rain_a
    elif main_field == a_field:
        rain_main = rain_a
        rain_secondary = rain_z
    elif main_field is None:
        main_field = a_field
        rain_main = rain_a
        rain_secondary = rain_z
    else:
        main_field = a_field
        rain_main = rain_a
        rain_secondary = rain_z
        thresh = 0.04
        thresh_max = False
        warn("Unknown main field. Using " + a_field + " with threshold " + str(thresh))

    if thresh_max:
        is_secondary = rain_main["data"] > thresh
    else:
        is_secondary = rain_main["data"] < thresh

    rain_main["data"][is_secondary] = rain_secondary["data"][is_secondary]

    return rain_main


def est_rain_rate_hydro(
    radar,
    alphazr=0.0376,
    betazr=0.6112,
    alphazs=0.1,
    betazs=0.5,
    alphaa=None,
    betaa=None,
    mp_factor=0.6,
    refl_field=None,
    a_field=None,
    hydro_field=None,
    rr_field=None,
    main_field=None,
    thresh=None,
    thresh_max=False,
):
    """
    Estimates rainfall rate using different relations between R and the
    polarimetric variables depending on the hydrometeor type.

    Parameters
    ----------
    radar : Radar
        Radar object.
    alphazr, betazr : floats, optional
        Factor (alpha) and exponent (beta) of the z-r power law for rain.
    alphazs, betazs : floats, optional
        Factor (alpha) and exponent (beta) of the z-s power law for snow.
    alphaa, betaa : floats, optional
        Factor (alpha) and exponent (beta) of the a-r power law.
        If not set the factors are going to be determined according
        to the radar frequency.
    mp_factor : float, optional
        Factor applied to z-r relation in the melting layer.
    refl_field : str, optional
        Name of the reflectivity field to use.
    a_field : str, optional
        Name of the specific attenuation field to use.
    hydro_field : str, optional
        Name of the hydrometeor classification field to use.
    rr_field : str, optional
        Name of the rainfall rate field.
    main_field : str, optional
        Name of the field that is going to act as main. Has to be
        either refl_field or kdp_field. Default is refl_field.
    thresh : float, optional
        Value of the threshold that determines when to use the secondary
        field.
    thresh_max : Bool, optional
        If true the main field is used up to the thresh value maximum.
        Otherwise the main field is not used below thresh value.

    Returns
    -------
    rain : dict
        Field dictionary containing the rainfall rate.

    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name("reflectivity")
    if a_field is None:
        a_field = get_field_name("specific_attenuation")
    if hydro_field is None:
        hydro_field = get_field_name("radar_echo_classification")
    if rr_field is None:
        rr_field = get_field_name("radar_estimated_rain_rate")

    # extract fields and parameters from radar
    if hydro_field in radar.fields:
        hydroclass = radar.fields[hydro_field]["data"]
    else:
        raise KeyError("Field not available: " + hydro_field)

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
    rain_z = est_rain_rate_z(
        radar, alpha=alphazr, beta=betazr, refl_field=refl_field, rr_field=rr_field
    )
    snow_z = est_rain_rate_z(
        radar, alpha=alphazs, beta=betazs, refl_field=refl_field, rr_field=rr_field
    )
    rain_a = est_rain_rate_a(
        radar, alpha=alphaa, beta=betaa, a_field=a_field, rr_field=rr_field
    )

    # initialize rainfall rate field
    rr_data = np.ma.zeros(hydroclass.shape, dtype="float32")
    rr_data[:] = np.ma.masked
    rr_data.set_fill_value(get_fillvalue())

    # apply the relations for each type
    # solid phase
    rr_data[is_ds] = snow_z["data"][is_ds]
    rr_data[is_cr] = snow_z["data"][is_cr]
    rr_data[is_vi] = snow_z["data"][is_vi]
    rr_data[is_gr] = snow_z["data"][is_gr]
    rr_data[is_ih] = snow_z["data"][is_ih]

    # rain
    if main_field == refl_field:
        rain_main = rain_z
        rain_secondary = rain_a
    elif main_field == a_field:
        rain_main = rain_a
        rain_secondary = rain_z
    elif main_field is None:
        main_field = a_field
        rain_main = rain_a
        rain_secondary = rain_z
    else:
        main_field = a_field
        rain_main = rain_a
        rain_secondary = rain_z
        thresh = 0.04
        thresh_max = False
        warn("Unknown main field. Using " + a_field + " with threshold " + str(thresh))

    if thresh_max:
        is_secondary = rain_main["data"] > thresh
    else:
        is_secondary = rain_main["data"] < thresh

    rain_main["data"][is_secondary] = rain_secondary["data"][is_secondary]

    rr_data[is_lr] = rain_main["data"][is_lr]
    rr_data[is_rn] = rain_main["data"][is_rn]

    # mixed phase
    rr_data[is_ws] = mp_factor * rain_z["data"][is_ws]
    rr_data[is_mh] = mp_factor * rain_z["data"][is_mh]

    rain = get_metadata(rr_field)
    rain["data"] = rr_data

    return rain


def _get_coeff_rkdp(freq):
    """
    Get the R(kdp) power law coefficients for a particular frequency.

    Parameters
    ----------
    freq : float
        Radar frequency [Hz].

    Returns
    -------
    alpha, beta : floats
        The coefficient and exponent of the power law.

    """
    coeff_rkdp_dict = _coeff_rkdp_table()

    freq_band = get_freq_band(freq)
    if (freq_band is not None) and (freq_band in coeff_rkdp_dict):
        return coeff_rkdp_dict[freq_band]

    if freq < 2e9:
        freq_band_aux = "S"
    elif freq > 12e9:
        freq_band_aux = "X"

    warn(
        "Radar frequency out of range. "
        + "Coefficients only applied to S, C or X band. "
        + freq_band_aux
        + " band coefficients will be used."
    )

    return coeff_rkdp_dict[freq_band_aux]


def _coeff_rkdp_table():
    """
    Defines the R(kdp) power law coefficients for each frequency band.

    Returns
    -------
    coeff_rkdp_dict : dict
        A dictionary with the coefficients at each band.

    """
    coeff_rkdp_dict = dict()

    # S band: Beard and Chuang coefficients
    coeff_rkdp_dict.update({"S": (50.70, 0.8500)})

    # C band: Beard and Chuang coefficients
    coeff_rkdp_dict.update({"C": (29.70, 0.8500)})

    # X band: Brandes coefficients
    coeff_rkdp_dict.update({"X": (15.81, 0.7992)})

    return coeff_rkdp_dict


def _get_coeff_ra(freq):
    """
    Get the R(A) power law coefficients for a particular frequency.

    Parameters
    ----------
    freq : float
        Radar frequency [Hz].

    Returns
    -------
    alpha, beta : floats
        The coefficient and exponent of the power law.

    """
    coeff_ra_dict = _coeff_ra_table()

    freq_band = get_freq_band(freq)
    if (freq_band is not None) and (freq_band in coeff_ra_dict):
        return coeff_ra_dict[freq_band]

    if freq < 2e9:
        freq_band_aux = "S"
    elif freq > 12e9:
        freq_band_aux = "X"

    warn(
        "Radar frequency out of range. "
        + "Coefficients only applied to S, C or X band. "
        + freq_band_aux
        + " band coefficients will be used."
    )

    return coeff_ra_dict[freq_band_aux]


def _coeff_ra_table():
    """
    Defines the R(A) power law coefficients for each frequency band.

    Returns
    -------
    coeff_ra_dict : dict
        A dictionary with the coefficients at each band.

    """
    coeff_ra_dict = dict()

    # S band: at 10 deg C according to tables from Ryzhkov et al. 2014
    coeff_ra_dict.update({"S": (3100.0, 1.03)})

    # C band: at 10 deg C according to tables from Diederich et al. 2015
    coeff_ra_dict.update({"C": (250.0, 0.91)})

    # X band: at 10 deg C according to tables from Diederich et al. 2015
    coeff_ra_dict.update({"X": (45.5, 0.83)})

    return coeff_ra_dict
