"""
Calculating attenuation from polarimetric radars

Code adapted from PAPER by Scott Giangrande et al

Adapted by Scott Collis and Scott Giangrande,

"""
import copy

import numpy as np
import scipy

from pyart.correct import phase_proc


def calculate_attenuation(radar, z_offset, **kwargs):
    doc = kwargs.get('doc', 15)
    fzl = kwargs.get('fzl', 4000.0)
    rhv_min = kwargs.get('rhv_min', 0.8)
    ncp_min = kwargs.get('ncp_min', 0.5)
    a_coef = kwargs.get('a_coef', 0.06)
    beta = kwargs.get('beta', 0.8)
    debug = kwargs.get('debug', False)
    is_cor = radar.fields['copol_coeff']['data'] > rhv_min
    is_coh = radar.fields['norm_coherent_power']['data'] > ncp_min
    is_good = np.logical_and(is_cor, is_cor)
    good_indeces = np.where(is_good)
    refl = np.ma.masked_where(np.logical_or(
        (radar.fields['norm_coherent_power']['data'] < ncp_min),
        (radar.fields['copol_coeff']['data'] < rhv_min)),
        copy.deepcopy(radar.fields['reflectivity_horizontal']['data']) +
        z_offset)
    ref_init_correct = (refl + radar.fields['proc_dp_phase_shift']['data'] *
                        a_coef)
    npts_good = len(good_indeces[0])
    dr = (radar.range['data'][1] - radar.range['data'][0]) / 1000.0
    specific_atten = np.zeros(
        radar.fields['reflectivity_horizontal']['data'].shape)
    atten = np.zeros(radar.fields['reflectivity_horizontal']['data'].shape)
    for sweep in range(len(radar.sweep_info['sweep_start_ray_index']['data'])):
        if debug:
            print "Doing ", sweep
        end_gate, start_ray, end_ray = phase_proc.det_process_range(
            radar, sweep, fzl, doc=doc)
        for i in range(start_ray, end_ray + 1):
            is_good = np.logical_and(is_cor[i, 0:end_gate],
                                     is_cor[i, 0:end_gate])
            good_indeces = np.where(is_good)
            maximum_phidp = np.median(
                radar.fields['proc_dp_phase_shift']['data'][i, 0:end_gate]
                [good_indeces[0][-6:]])
            smoothed_reflectivity = phase_proc.smooth_and_trim(
                ref_init_correct[i, 0:end_gate], window_len=5)
            reflectivity_in_linear_units = 10.0 ** (0.1 * beta *
                                                    smoothed_reflectivity)
            I_indef = scipy.integrate.cumtrapz(
                0.46 * beta * dr * reflectivity_in_linear_units[::-1])
            I_indef = np.append(I_indef, I_indef[-1])[::-1]
            self_cons_number = 10.0 ** (0.1 * beta * a_coef *
                                        maximum_phidp) - 1.0
            specific_atten[i, 0:end_gate] = (
                reflectivity_in_linear_units * self_cons_number /
                (I_indef[0] + self_cons_number*(I_indef)))
            atten[i, :-1] = scipy.integrate.cumtrapz(
                specific_atten[i, :]) * dr * 2.0
            atten[i, -1] = atten[i, -2]

    spec_at = copy.deepcopy(radar.fields[radar.fields.keys()[0]])
    spec_at['data'] = specific_atten
    spec_at['valid_min'] = 0.0
    spec_at['valid_max'] = 1.0
    spec_at['standard_name'] = 'specific_attenuation'
    spec_at['long_name'] = 'specific_attenuation'
    spec_at['least_significant_digit'] = 4
    spec_at['units'] = 'dB/km'
    cor_z = copy.deepcopy(radar.fields['reflectivity_horizontal'])
    cor_z['data'] = ref_init_correct  # atten+cor_z['data']+z_offset
    cor_z['least_significant_digit'] = 2
    cor_z['valid_max'] = 80.0
    return spec_at, cor_z
