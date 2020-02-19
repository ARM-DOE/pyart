"""
Attenuation correction from polarimetric radars. Code adapted from method in
Gu et al, JAMC 2011, 50, 39. Adapted by Scott Collis and Scott Giangrande,
refactored by Jonathan Helmus. New code added by Meteo Swiss and inserted into
Py-ART by Robert Jackson.

"""

from copy import deepcopy
from warnings import warn

import numpy as np
from scipy.integrate import cumtrapz

from ..config import get_metadata, get_field_name, get_fillvalue
from .phase_proc import smooth_masked, det_process_range, smooth_and_trim
from ..filters import temp_based_gate_filter, iso0_based_gate_filter
from ..retrieve import get_freq_band


def calculate_attenuation_zphi(radar, doc=None, fzl=None, smooth_window_len=5,
                               gatefilter=None, a_coef=None, beta=None, c=None,
                               d=None, refl_field=None, phidp_field=None,
                               zdr_field=None, temp_field=None,
                               iso0_field=None, spec_at_field=None,
                               pia_field=None, corr_refl_field=None,
                               spec_diff_at_field=None, pida_field=None,
                               corr_zdr_field=None, temp_ref='temperature'):
    """
    Calculate the attenuation and the differential attenuation from a
    polarimetric radar using Z-PHI method..
    The attenuation is computed up to a user defined freezing level height
    or up to where temperatures in a temperature field are positive.
    The coefficients are either user-defined or radar frequency dependent.

    Parameters
    ----------
    radar : Radar
        Radar object to use for attenuation calculations. Must have
        phidp and refl fields.
    doc : float, optional
        Number of gates at the end of each ray to to remove from the
        calculation.
    fzl : float, optional
        Freezing layer, gates above this point are not included in the
        correction.
    gatefilter : GateFilter, optional
        The gates to exclude from the calculation. This, combined with
        the gates above fzl, will be excluded from the correction. Set to
        None to not use a gatefilter.
    smooth_window_len : int, optional
        Size, in range bins, of the smoothing window
    a_coef : float, optional
        A coefficient in attenuation calculation.
    beta : float, optional
        Beta parameter in attenuation calculation.
    c, d : float, optional
        coefficient and exponent of the power law that relates attenuation
        with differential attenuation
    refl_field : str, optional
        Name of the reflectivity field used for the attenuation correction.
        A value of None for any of these parameters will use the default
        field name as defined in the Py-ART configuration file.
    phidp_field : str, optional
        Name of the differential phase field used for the attenuation
        correction. A value of None for any of these parameters will use the
        default field name as defined in the Py-ART configuration file.
    zdr_field : str, optional
        Name of the differential reflectivity field used for the attenuation
        correction. A value of None for any of these parameters will use the
        default field name as defined in the Py-ART configuration file. This
        will only be used if it is available.
    temp_field : str, optional
        Name of the temperature field used for the attenuation
        correction. A value of None for any of these parameters will use the
        default field name as defined in the Py-ART configuration file.
    iso0_field : str, optional
        Name of the field for the height above the 0C isotherm for the
        attenuation correction. A value of None for any of these parameters
        will use the default field name as defined in the Py-ART configuration
        file. This will only be used if it is available.
    spec_at_field : str, optional
        Name of the specific attenuation field that will be used to fill in
        the metadata for the returned fields. A value of None for any of these
        parameters will use the default field names as defined in the Py-ART
        configuration file.
    pia_field : str, optional
        Name of the path integrated attenuation field that will be used to fill
        in the metadata for the returned fields. A value of None for any of
        these parameters will use the default field names as defined in the
        Py-ART configuration file.
    corr_refl_field : str, optional
        Name of the corrected reflectivity field that will be used to fill in
        the metadata for the returned fields. A value of None for any of these
        parameters will use the default field names as defined in the Py-ART
        configuration file.
    spec_diff_at_field : str, optional
        Name of the specific differential attenuation field that will be used
        to fill in the metadata for the returned fields. A value of None for
        any of these parameters will use the default field names as defined
        in the Py-ART configuration file. This will only be calculated if ZDR
        is available.
    pida_field : str, optional
        Name of the path integrated differential attenuation field that will
        be used to fill in the metadata for the returned fields. A value of
        None for any of these parameters will use the default field names as
        defined in the Py-ART configuration file. This will only be calculated
        if ZDR is available.
    corr_zdr_field : str, optional
        Name of the corrected differential reflectivity field that will
        be used to fill in the metadata for the returned fields. A value of
        None for any of these parameters will use the default field names as
        defined in the Py-ART configuration file. This will only be calculated
        if ZDR is available.
    temp_ref : str, optional
        the field use as reference for temperature. Can be either temperature,
        height_over_iso0 or fixed_fzl

    Returns
    -------
    spec_at : dict
        Field dictionary containing the specific attenuation.
    pia_dict : dict
        Field dictionary containing the path integrated attenuation.
    cor_z : dict
        Field dictionary containing the corrected reflectivity.
    spec_diff_at : dict
        Field dictionary containing the specific differential attenuation.
    pida_dict : dict
        Field dictionary containing the path integrated differential
        attenuation.
    cor_zdr : dict
        Field dictionary containing the corrected differential reflectivity.

    References
    ----------
    Gu et al. Polarimetric Attenuation Correction in Heavy Rain at C Band,
    JAMC, 2011, 50, 39-58.

    Ryzhkov et al. Potential Utilization of Specific Attenuation for Rainfall
    Estimation, Mitigation of Partial Beam Blockage, and Radar Networking,
    JAOT, 2014, 31, 599-619.

    """
    # select the coefficients as a function of frequency band
    if (a_coef is None) or (beta is None) or (c is None) or (d is None):
        if 'frequency' in radar.instrument_parameters:
            a_coef, beta, c, d = _get_param_attzphi(
                radar.instrument_parameters['frequency']['data'][0])
        else:
            a_coef, beta, c, d = _param_attzphi_table()['C']
            warn("Radar frequency unknown. Default coefficients "
                 + "for C band will be applied.")

    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if zdr_field is None:
        zdr_field = get_field_name('differential_reflectivity')
    if phidp_field is None:
        # use corrrected_differential_phase or unfolded_differential_phase
        # fields if they are available, if not use differential_phase field
        phidp_field = get_field_name('corrected_differential_phase')
        if phidp_field not in radar.fields:
            phidp_field = get_field_name('unfolded_differential_phase')
        if phidp_field not in radar.fields:
            phidp_field = get_field_name('differential_phase')
    if spec_at_field is None:
        spec_at_field = get_field_name('specific_attenuation')
    if pia_field is None:
        pia_field = get_field_name('path_integrated_attenuation')
    if corr_refl_field is None:
        corr_refl_field = get_field_name('corrected_reflectivity')
    if spec_diff_at_field is None:
        spec_diff_at_field = get_field_name(
            'specific_differential_attenuation')
    if pida_field is None:
        pida_field = get_field_name('path_integrated_differential_attenuation')
    if corr_zdr_field is None:
        corr_zdr_field = get_field_name(
            'corrected_differential_reflectivity')

    if temp_ref == 'temperature':
        if temp_field is None:
            temp_field = get_field_name('temperature')
    elif temp_ref == 'height_over_iso0':
        if iso0_field is None:
            iso0_field = get_field_name('height_over_iso0')

    # extract fields and parameters from radar if they exist
    # reflectivity and differential phase must exist
    # create arrays to hold the output data
    radar.check_field_exists(refl_field)
    refl = radar.fields[refl_field]['data']

    radar.check_field_exists(phidp_field)
    phidp = deepcopy(radar.fields[phidp_field]['data'])

    ah = np.ma.zeros(refl.shape, dtype='float64')
    pia = np.ma.zeros(refl.shape, dtype='float64')

    try:
        radar.check_field_exists(zdr_field)
        zdr = radar.fields[zdr_field]['data']

        adiff = np.ma.zeros(zdr.shape, dtype='float64')
        pida = np.ma.zeros(zdr.shape, dtype='float64')
    except KeyError:
        zdr = None

    # determine the valid data (i.e. data below freezing level)
    mask_fzl, end_gate_arr = get_mask_fzl(
        radar, fzl=fzl, doc=doc, min_temp=0, max_h_iso0=0., thickness=None,
        beamwidth=None, temp_field=temp_field, iso0_field=iso0_field,
        temp_ref=temp_ref)

    if gatefilter is None:
        mask = np.ma.getmaskarray(refl)
    else:
        mask = gatefilter.gate_excluded
        mask_fzl = np.logical_or(mask, mask_fzl)

    # prepare phidp: filter out values above freezing level and negative
    # makes sure phidp is monotonously increasing
    corr_phidp = _prepare_phidp(phidp, mask_fzl)

    # calculate initial reflectivity correction and gate spacing (in km)
    init_refl_correct = refl + corr_phidp * a_coef
    dr = (radar.range['data'][1] - radar.range['data'][0]) / 1000.0

    if smooth_window_len > 0:
        sm_refl = smooth_masked(init_refl_correct, wind_len=smooth_window_len,
                                min_valid=1, wind_type='mean')
    else:
        sm_refl = init_refl_correct
    refl_linear = np.ma.power(10.0, 0.1 * beta * sm_refl).filled(fill_value=0)

    for ray in range(radar.nrays):
        # perform attenuation calculation on a single ray
        # if number of valid range bins larger than smoothing window
        if end_gate_arr[ray] < 0:
            continue

        if end_gate_arr[ray] > smooth_window_len:
            # extract the ray's phase shift,
            # init. refl. correction and mask
            ray_phase_shift = corr_phidp[ray, 0:end_gate_arr[ray]]
            ray_mask = mask[ray, 0:end_gate_arr[ray]]
            ray_refl_linear = refl_linear[ray, 0:end_gate_arr[ray]]

            # perform calculation if there is valid data
            last_six_good = np.where(
                np.ndarray.flatten(ray_mask) == 0)[0][-6:]
            if(len(last_six_good)) == 6:
                phidp_max = np.median(ray_phase_shift[last_six_good])
                self_cons_number = (
                    10.0 ** (0.1 * beta * a_coef * phidp_max) - 1.0)
                I_indef = cumtrapz(0.46 * beta * dr * ray_refl_linear[::-1])
                I_indef = np.append(I_indef, I_indef[-1])[::-1]

                # set the specific attenutation and attenuation
                ah[ray, 0:end_gate_arr[ray]] = (
                    ray_refl_linear * self_cons_number /
                    (I_indef[0] + self_cons_number * I_indef))

                pia[ray, :-1] = cumtrapz(ah[ray, :]) * dr * 2.0
                pia[ray, -1] = pia[ray, -2]

                # if ZDR exists, set the specific differential attenuation
                # and differential attenuation
                if zdr is not None:
                    adiff[ray, 0:end_gate_arr[ray]] = (
                        c * np.ma.power(ah[ray, 0:end_gate_arr[ray]], d))

                    pida[ray, :-1] = cumtrapz(adiff[ray, :]) * dr * 2.0
                    pida[ray, -1] = pida[ray, -2]

    # prepare output field dictionaries
    # for specific attenuation and corrected reflectivity
    spec_at = get_metadata(spec_at_field)
    temp_array = np.ma.masked_where(mask, ah)
    spec_at['data'] = temp_array
    spec_at['_FillValue'] = temp_array.fill_value

    pia_array = np.ma.masked_where(mask, pia)
    pia_dict = get_metadata(pia_field)
    pia_dict['data'] = pia_array
    pia_dict['_FillValue'] = pia_array.fill_value

    cor_z = get_metadata(corr_refl_field)
    cor_z_array = np.ma.masked_where(mask, pia + refl)
    cor_z['data'] = cor_z_array
    cor_z['_FillValue'] = cor_z_array.fill_value

    # prepare output field dictionaries
    # for specific diff attenuation and corrected ZDR
    if zdr is not None:
        sda = np.ma.masked_where(mask, adiff)
        spec_diff_at = get_metadata(spec_diff_at_field)
        spec_diff_at['data'] = sda
        spec_diff_at['_FillValue'] = sda.fill_value

        pida_array = np.ma.masked_where(mask, pida)
        pida_dict = get_metadata(pida_field)
        pida_dict['data'] = pida_array
        pida_dict['_FillValue'] = pida_array.fill_value

        cor_zdr = get_metadata(corr_zdr_field)
        czdr = np.ma.masked_where(mask, pida + zdr)
        cor_zdr['data'] = czdr
        cor_zdr['_FillValue'] = czdr.fill_value
    else:
        spec_diff_at = None
        cor_zdr = None
        pida_dict = None

    return spec_at, pia_dict, cor_z, spec_diff_at, pida_dict, cor_zdr


def calculate_attenuation_philinear(
        radar, doc=None, fzl=None, pia_coef=None, gatefilter=None,
        pida_coef=None, refl_field=None, phidp_field=None, zdr_field=None,
        temp_field=None, iso0_field=None, spec_at_field=None,
        pia_field=None, corr_refl_field=None, spec_diff_at_field=None,
        pida_field=None, corr_zdr_field=None, temp_ref='temperature'):
    """
    Calculate the attenuation and the differential attenuation from a
    polarimetric radar using linear dependece with PhiDP.
    The attenuation is computed up to a user defined freezing level height,
    where temperatures in a temperature field are positive or where the height
    relative to the iso0 is 0.
    The coefficients are either user-defined or radar frequency dependent.

    Parameters
    ----------
    radar : Radar
        Radar object to use for attenuation calculations. Must have
        phidp and refl fields.
    doc : float, optional
        Number of gates at the end of each ray to to remove from the
        calculation.
    fzl : float, optional
        Freezing layer, gates above this point are not included in the
        correction.
    gatefilter : GateFilter, optional
        The gates to exclude from the calculation. This, combined with
        the gates above fzl, will be excluded from the correction. Set to
        None to not use a gatefilter.
    pia_coef : float, optional
        Coefficient in path integrated attenuation calculation
    pida_coeff : float, optional
        Coefficient in path integrated differential attenuation calculation
    refl_field : str, optional
        Name of the reflectivity field used for the attenuation correction.
        A value of None for any of these parameters will use the default
        field name as defined in the Py-ART configuration file.
    phidp_field : str, optional
        Name of the differential phase field used for the attenuation
        correction. A value of None for any of these parameters will use the
        default field name as defined in the Py-ART configuration file.
    zdr_field : str, optional
        Name of the differential reflectivity field used for the attenuation
        correction. A value of None for any of these parameters will use the
        default field name as defined in the Py-ART configuration file. This
        will only be used if it is available.
    temp_field : str, optional
        Name of the temperature field used for the attenuation
        correction. A value of None for any of these parameters will use the
        default field name as defined in the Py-ART configuration file.
    iso0_field : str, optional
        Name of the field for the height above the 0C isotherm for the
        attenuation correction. A value of None for any of these parameters
        will use the default field name as defined in the Py-ART configuration
        file. This will only be used if it is available.
    spec_at_field : str, optional
        Name of the specific attenuation field that will be used to fill in
        the metadata for the returned fields. A value of None for any of these
        parameters will use the default field names as defined in the Py-ART
        configuration file.
    pia_field : str, optional
        Name of the path integrated attenuation field that will be used to fill
        in the metadata for the returned fields. A value of None for any of
        these parameters will use the default field names as defined in the
        Py-ART configuration file.
    corr_refl_field : str, optional
        Name of the corrected reflectivity field that will be used to fill in
        the metadata for the returned fields. A value of None for any of these
        parameters will use the default field names as defined in the Py-ART
        configuration file.
    spec_diff_at_field : str, optional
        Name of the specific differential attenuation field that will be used
        to fill in the metadata for the returned fields. A value of None for
        any of these parameters will use the default field names as defined
        in the Py-ART configuration file. This will only be calculated if ZDR
        is available.
    corr_zdr_field : str, optional
        Name of the corrected differential reflectivity field that will
        be used to fill in the metadata for the returned fields. A value of
        None for any of these parameters will use the default field names as
        defined in the Py-ART configuration file. This will only be calculated
        if ZDR is available.
    temp_ref : str, optional
        The field use as reference for temperature. Can be either temperature,
        height_over_iso0 or fixed_fzl.

    Returns
    -------
    spec_at : dict
        Field dictionary containing the specific attenuation.
    pia_dict : dict
        Field dictionary containing the path integrated attenuation.
    cor_z : dict
        Field dictionary containing the corrected reflectivity.
    spec_diff_at : dict
        Field dictionary containing the specific differential attenuation.
    pida_dict : dict
        Field dictionary containing the path integrated differential
        attenuation.
    cor_zdr : dict
        Field dictionary containing the corrected differential reflectivity.

    """
    # select the coefficients as a function of frequency band
    if (pia_coef is None) or (pida_coef is None):
        if 'frequency' in radar.instrument_parameters:
            pia_coef, pida_coef = _get_param_attphilinear(
                radar.instrument_parameters['frequency']['data'][0])
        else:
            pia_coef, pida_coef = _param_attphilinear_table()['C']
            warn("Radar frequency unknown. Default "
                 + "coefficients for C band will be applied.")

    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if zdr_field is None:
        zdr_field = get_field_name('differential_reflectivity')
    if phidp_field is None:
        # use corrrected_differential_phase or unfolded_differential_phase
        # fields if they are available, if not use differential_phase field
        phidp_field = get_field_name('corrected_differential_phase')
        if phidp_field not in radar.fields:
            phidp_field = get_field_name('unfolded_differential_phase')
        if phidp_field not in radar.fields:
            phidp_field = get_field_name('differential_phase')
    if spec_at_field is None:
        spec_at_field = get_field_name('specific_attenuation')
    if pia_field is None:
        pia_field = get_field_name('path_integrated_attenuation')
    if corr_refl_field is None:
        corr_refl_field = get_field_name('corrected_reflectivity')
    if spec_diff_at_field is None:
        spec_diff_at_field = get_field_name(
            'specific_differential_attenuation')
    if pida_field is None:
        pida_field = get_field_name(
            'path_integrated_differential_attenuation')
    if corr_zdr_field is None:
        corr_zdr_field = get_field_name(
            'corrected_differential_reflectivity')

    if temp_ref == 'temperature':
        if temp_field is None:
            temp_field = get_field_name('temperature')
    elif temp_ref == 'height_over_iso0':
        if iso0_field is None:
            iso0_field = get_field_name('height_over_iso0')

    # extract fields and parameters from radar if they exist
    # reflectivity and differential phase must exist
    # create arrays to hold the output data
    radar.check_field_exists(refl_field)
    refl = radar.fields[refl_field]['data']

    radar.check_field_exists(phidp_field)
    phidp = deepcopy(radar.fields[phidp_field]['data'])

    try:
        radar.check_field_exists(zdr_field)
        zdr = radar.fields[zdr_field]['data']
    except KeyError:
        zdr = None

    # determine the valid data (i.e. data below freezing level)
    mask_fzl, _ = get_mask_fzl(
        radar, fzl=fzl, doc=doc, min_temp=0, max_h_iso0=0., thickness=None,
        beamwidth=None, temp_field=temp_field, iso0_field=iso0_field,
        temp_ref=temp_ref)

    if gatefilter is None:
        mask = np.ma.getmaskarray(refl)
    else:
        mask = gatefilter.gate_excluded
        mask_fzl = np.logical_or(mask, mask_fzl)

    # prepare phidp: filter out values above freezing level and negative
    # makes sure phidp is monotonously increasing
    corr_phidp = _prepare_phidp(phidp, mask_fzl)
    dr = (radar.range['data'][1] - radar.range['data'][0]) / 1000.0

    pia = pia_coef * corr_phidp
    ah = 0.5 * np.gradient(pia, dr, axis=1)

    # prepare output field dictionaries
    # for specific attenuation and corrected reflectivity
    spec_at = get_metadata(spec_at_field)
    spec_at['data'] = np.ma.masked_where(mask, np.ma.array(ah))

    pia_dict = get_metadata(pia_field)
    pia_dict['data'] = np.ma.masked_where(mask, np.ma.array(pia))

    cor_z = get_metadata(corr_refl_field)
    cor_z['data'] = np.ma.masked_where(mask, np.ma.array(pia + refl))

    # prepare output field dictionaries
    # for specific diff attenuation and corrected ZDR
    if zdr is not None:
        pida = pida_coef * corr_phidp
        adiff = 0.5 * np.gradient(pida, dr, axis=1)

        spec_diff_at = get_metadata(spec_diff_at_field)
        spec_diff_at['data'] = np.ma.masked_where(mask, adiff)

        pida_dict = get_metadata(pida_field)
        pida_dict['data'] = np.ma.masked_where(mask, pida)

        cor_zdr = get_metadata(corr_zdr_field)
        cor_zdr['data'] = np.ma.masked_where(mask, pida + zdr)

    return spec_at, pia_dict, cor_z, spec_diff_at, pida_dict, cor_zdr


def get_mask_fzl(radar, fzl=None, doc=None, min_temp=0., max_h_iso0=0.,
                 thickness=None, beamwidth=None, temp_field=None,
                 iso0_field=None, temp_ref='temperature'):
    """
    Constructs a mask to mask data placed thickness m below data at min_temp
    and beyond.

    Parameters
    ----------
    radar : Radar
        The radar object.
    fzl : float, optional
        Freezing layer, gates above this point are not included in the
        correction.
    doc : float, optional
        Number of gates at the end of each ray to to remove from the
        calculation.
    min_temp : float, optional
        Minimum temperature below which the data is mask in degrees.
    max_h_iso0 : float, optional
        Maximum height relative to the iso0 below which the data is mask in
        meters.
    thickness : float, optional
        Extent of the layer below the first gate where min_temp is reached
        that is going to be masked.
    beamwidth : float, optional
        The radar antenna 3 dB beamwidth.
    temp_field: str, optional
        The temperature field. A value of None will use the default
        field name as defined in the Py-ART configuration file. It is going
        to be used only if available.
    iso0_field: str, optional
        The field containing the height over the 0C isotherm. A value of None
        will use the default field name as defined in the Py-ART
        configuration file. It is going to be used only if available.
    temp_ref : str, optional
        The field use as reference for temperature. Can be either temperature,
        height_over_iso0 or fixed_fzl.

    Returns
    -------
    mask_fzl : 2D array
        The values that should be masked.
    end_gate_arr : 1D array
        The index of the last valid gate in the ray.

    """
    if temp_ref == 'temperature':
        if temp_field is None:
            temp_field = get_field_name('temperature')
    elif temp_ref == 'height_over_iso0':
        if iso0_field is None:
            iso0_field = get_field_name('height_over_iso0')

    if temp_ref == 'fixed_fzl':
        if fzl is None:
            fzl = 4000.
            doc = 15
            warn('Freezing level height not specified. ' +
                 'Using default '+str(fzl)+' [m]')
        end_gate_arr = np.zeros(radar.nrays, dtype='int32')
        mask_fzl = np.zeros((radar.nrays, radar.ngates), dtype=np.bool)
        for sweep in range(radar.nsweeps):
            end_gate, start_ray, end_ray = (
                det_process_range(radar, sweep, fzl, doc=doc))
            end_gate_arr[start_ray:end_ray] = end_gate
            mask_fzl[start_ray:end_ray, end_gate+1:] = True

    elif temp_ref == 'temperature':
        if temp_field in radar.fields:
            gatefilter = temp_based_gate_filter(
                radar, temp_field=temp_field, min_temp=min_temp,
                thickness=thickness, beamwidth=beamwidth)
            end_gate_arr = np.zeros(radar.nrays, dtype='int32')
            for ray in range(radar.nrays):
                ind_rng = np.where(gatefilter.gate_excluded[ray, :] == 1)[0]
                if len(ind_rng) > 0:
                    # there are filtered gates: The last valid gate is one
                    # before the first filter gate
                    if ind_rng[0] > 0:
                        end_gate_arr[ray] = ind_rng[0]-1
                    else:
                        end_gate_arr[ray] = 0
                else:
                    # there are no filter gates: all gates are valid
                    end_gate_arr[ray] = radar.ngates-1
            mask_fzl = gatefilter.gate_excluded == 1
        else:
            fzl = 4000.
            doc = 15
            warn('Temperature field not available.' +
                 'Using default freezing level height ' +
                 str(fzl) + ' [m].')
    else:
        if iso0_field in radar.fields:
            gatefilter = iso0_based_gate_filter(
                radar, iso0_field=iso0_field, max_h_iso0=max_h_iso0,
                thickness=thickness, beamwidth=beamwidth)
            end_gate_arr = np.zeros(radar.nrays, dtype='int32')
            for ray in range(radar.nrays):
                ind_rng = np.where(gatefilter.gate_excluded[ray, :] == 1)[0]
                if len(ind_rng) > 0:
                    # there are filtered gates: The last valid gate is one
                    # before the first filter gate
                    if ind_rng[0] > 0:
                        end_gate_arr[ray] = ind_rng[0]-1
                    else:
                        end_gate_arr[ray] = 0
                else:
                    # there are no filter gates: all gates are valid
                    end_gate_arr[ray] = radar.ngates-1
            mask_fzl = gatefilter.gate_excluded == 1
        else:
            fzl = 4000.
            doc = 15
            warn("Height over iso0 field not available."
                 + "Using default freezing level height "
                 + str(fzl) + " [m].")

    return mask_fzl, end_gate_arr


def _prepare_phidp(phidp, mask_fzl):
    """
    Prepares phidp to be used in attenuation correction by masking values
    above freezing level setting negative values to 0 and make sure it is
    monotously increasing.

    Parameters
    ----------
    phidp : ndarray 2D
        The phidp field.
    mask_fzl : ndarray 2D
        A mask of the data above freezing level height.

    Returns
    -------
    corr_phidp: ndarray 2D
        The corrected PhiDP field.

    """
    mask_phidp = np.ma.getmaskarray(phidp)
    mask_phidp = np.logical_or(mask_phidp, mask_fzl)
    mask_phidp = np.logical_or(mask_phidp, phidp < 0.)
    corr_phidp = np.ma.masked_where(mask_phidp, phidp)

    return np.maximum.accumulate(corr_phidp.filled(fill_value=0.), axis=1)


def _get_param_attzphi(freq):
    """
    Get the parameters of Z-Phi attenuation estimation for a particular
    frequency.

    Parameters
    ----------
    freq : float
        Radar frequency [Hz].

    Returns
    -------
    a_coeff, beta, c, d : floats
        The coefficient and exponent of the power law.

    """
    param_att_dict = _param_attzphi_table()

    freq_band = get_freq_band(freq)
    if (freq_band is not None) and (freq_band in param_att_dict):
        return param_att_dict[freq_band]

    if freq < 2e9:
        freq_band_aux = 'S'
    elif freq > 12e9:
        freq_band_aux = 'X'

    warn("Radar frequency out of range. "
         + "Coefficients only applied to S, C or X band. "
         + freq_band + " band coefficients will be used.")

    return param_att_dict[freq_band_aux]


def _param_attzphi_table():
    """
    Defines the parameters of Z-Phi attenuation estimation at each frequency
    band.

    Returns
    -------
    param_att_dict : dict
        A dictionary with the coefficients at each band.

    """
    param_att_dict = dict()

    # S band:
    param_att_dict.update({'S': (0.02, 0.64884, 0.15917, 1.0804)})

    # C band:
    param_att_dict.update({'C': (0.08, 0.64884, 0.3, 1.0804)})

    # X band:
    param_att_dict.update({'X': (0.31916, 0.64884, 0.15917, 1.0804)})

    return param_att_dict


def _get_param_attphilinear(freq):
    """
    Get the parameters of attenuation estimation based on phidp for a
    particular frequency.

    Parameters
    ----------
    freq : float
        Radar frequency [Hz].

    Returns
    -------
    a_coeff, beta, c, d : floats
        The coefficient and exponent of the power law.

    """
    param_att_dict = _param_attphilinear_table()

    freq_band = get_freq_band(freq)
    if (freq_band is not None) and (freq_band in param_att_dict):
        return param_att_dict[freq_band]

    if freq < 2e9:
        freq_band_aux = 'S'
    elif freq > 12e9:
        freq_band_aux = 'X'

    warn("Radar frequency out of range. "
         + "Coefficients only applied to S, C or X band. "
         + freq_band_aux + " band coefficients will be used.")

    return param_att_dict[freq_band_aux]


def _param_attphilinear_table():
    """
    Defines the parameters of attenuation estimation based on phidp at each
    frequency band.

    Returns
    -------
    param_att_dict : dict
        A dictionary with the coefficients at each band.

    """
    param_att_dict = dict()

    # S band:
    param_att_dict.update({'S': (0.04, 0.004)})

    # C band:
    param_att_dict.update({'C': (0.08, 0.03)})

    # X band:
    param_att_dict.update({'X': (0.28, 0.04)})

    return param_att_dict


def calculate_attenuation(radar, z_offset, debug=False, doc=15, fzl=4000.0,
                          rhv_min=0.8, ncp_min=0.5, a_coef=0.06, beta=0.8,
                          refl_field=None, ncp_field=None, rhv_field=None,
                          phidp_field=None, spec_at_field=None,
                          corr_refl_field=None):
    """
    Calculate the attenuation from a polarimetric radar using Z-PHI method.

    Parameters
    ----------
    radar : Radar
        Radar object to use for attenuation calculations. Must have
        copol_coeff, norm_coherent_power, proc_dp_phase_shift,
        reflectivity_horizontal fields.
    z_offset : float
        Horizontal reflectivity offset in dBZ.
    debug : bool, optional
        True to print debugging information, False supressed this printing.
    doc : float, optional
        Number of gates at the end of each ray to to remove from the
        calculation.
    fzl : float, optional
        Freezing layer, gates above this point are not included in the
        correction.
    rhv_min : float, optional
        Minimum copol_coeff value to consider valid.
    ncp_min : float, optional
        Minimum norm_coherent_power to consider valid.
    a_coef : float, optional
        A coefficient in attenuation calculation.
    beta : float, optional
        Beta parameter in attenuation calculation.
    refl_field : str, optional
        Name of the reflectivity field used for the attenuation correction.
        A value of None for any of these parameters will use the default
        field name as defined in the Py-ART configuration file.
    phidp_field : str, optional
        Name of the differential phase field used for the attenuation
        correction. A value of None for any of these parameters will use the
        default field name as defined in the Py-ART configuration file.
    ncp_field : str, optional
        Name of the normalized coherent power field used for the attenuation
        correction. A value of None for any of these parameters will use the
        default field name as defined in the Py-ART configuration file.
    zdr_field : str, optional
        Name of the differential reflectivity field used for the attenuation
        correction. A value of None for any of these parameters will use the
        default field name as defined in the Py-ART configuration file. This
        will only be used if it is available.
    spec_at_field : str, optional
        Name of the specific attenuation field that will be used to fill in
        the metadata for the returned fields. A value of None for any of these
        parameters will use the default field names as defined in the Py-ART
        configuration file.
    corr_refl_field : str, optional
        Name of the corrected reflectivity field that will be used to fill in
        the metadata for the returned fields. A value of None for any of these
        parameters will use the default field names as defined in the Py-ART
        configuration file.

    Returns
    -------
    spec_at : dict
        Field dictionary containing the specific attenuation.
    cor_z : dict
        Field dictionary containing the corrected reflectivity.

    References
    ----------
    Gu et al. Polarimetric Attenuation Correction in Heavy Rain at C Band,
    JAMC, 2011, 50, 39-58.

    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')
    if rhv_field is None:
        rhv_field = get_field_name('cross_correlation_ratio')
    if phidp_field is None:
        # use corrrected_differential_phae or unfolded_differential_phase
        # fields if they are available, if not use differential_phase field
        phidp_field = get_field_name('corrected_differential_phase')
        if phidp_field not in radar.fields:
            phidp_field = get_field_name('unfolded_differential_phase')
        if phidp_field not in radar.fields:
            phidp_field = get_field_name('differential_phase')
    if spec_at_field is None:
        spec_at_field = get_field_name('specific_attenuation')
    if corr_refl_field is None:
        corr_refl_field = get_field_name('corrected_reflectivity')

    # extract fields and parameters from radar
    norm_coherent_power = radar.fields[ncp_field]['data']
    copol_coeff = radar.fields[rhv_field]['data']
    reflectivity_horizontal = radar.fields[refl_field]['data']
    proc_dp_phase_shift = radar.fields[phidp_field]['data']
    nsweeps = int(radar.nsweeps)

    # determine where the reflectivity is valid, mask out bad locations.
    is_cor = copol_coeff > rhv_min
    is_coh = norm_coherent_power > ncp_min
    is_good = np.logical_and(is_cor, is_coh)
    mask = np.logical_not(is_good)
    refl = np.ma.masked_where(mask, reflectivity_horizontal + z_offset)

    # calculate initial reflectivity correction and gate spacing (in km)
    init_refl_correct = refl + proc_dp_phase_shift * a_coef
    dr = (radar.range['data'][1] - radar.range['data'][0]) / 1000.0

    # create array to hold specific attenuation and attenuation
    specific_atten = np.zeros(reflectivity_horizontal.shape, dtype='float32')
    atten = np.zeros(reflectivity_horizontal.shape, dtype='float32')

    for sweep in range(nsweeps):
        # loop over the sweeps
        if debug:
            print("Doing ", sweep)
        end_gate, start_ray, end_ray = det_process_range(
            radar, sweep, fzl, doc=doc)

        for i in range(start_ray, end_ray):
            # perform attenuation calculation on a single ray

            # extract the ray's phase shift and init. refl. correction
            ray_phase_shift = proc_dp_phase_shift[i, 0:end_gate]
            ray_init_refl = init_refl_correct[i, 0:end_gate]

            # perform calculation
            last_six_good = np.where(is_good[i, 0:end_gate])[0][-6:]
            phidp_max = np.median(ray_phase_shift[last_six_good])
            sm_refl = smooth_and_trim(ray_init_refl, window_len=5)
            reflectivity_linear = 10.0 ** (0.1 * beta * sm_refl)
            self_cons_number = 10.0 ** (0.1 * beta * a_coef * phidp_max) - 1.0
            I_indef = cumtrapz(0.46 * beta * dr * reflectivity_linear[::-1])
            I_indef = np.append(I_indef, I_indef[-1])[::-1]

            # set the specific attenutation and attenuation
            specific_atten[i, 0:end_gate] = (
                reflectivity_linear * self_cons_number /
                (I_indef[0] + self_cons_number * I_indef))

            atten[i, :-1] = cumtrapz(specific_atten[i, :]) * dr * 2.0
            atten[i, -1] = atten[i, -2]

    # prepare output field dictionaries
    spec_at = get_metadata(spec_at_field)
    spec_at['data'] = specific_atten
    spec_at['_FillValue'] = get_fillvalue()

    cor_z = get_metadata(corr_refl_field)
    cor_z['data'] = atten + reflectivity_horizontal + z_offset
    cor_z['data'].mask = init_refl_correct.mask
    cor_z['_FillValue'] = get_fillvalue()

    return spec_at, cor_z
