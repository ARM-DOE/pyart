"""
pyart.correct.attenuation
=========================

Attenuation correction from polarimetric radars.

Code adapted from method in Gu et al, JAMC 2011, 50, 39.

Adapted by Scott Collis and Scott Giangrande, refactored by Jonathan Helmus.

.. autosummary::
    :toctree: generated/

    calculate_attenuation


"""
import copy

import numpy as np
from scipy.integrate import cumtrapz

from . import phase_proc


def calculate_attenuation(radar, z_offset, debug=False, doc=15, fzl=4000.0,
                          rhv_min=0.8, ncp_min=0.5, a_coef=0.06, beta=0.8):
    """
    Calculate the attenuation from a polarimetric radar using Z-PHI method.

    Parameters
    ----------
    radar : Radar
        Radar object to use for attenuation calculations.  Must have
        copol_coeff, norm_coherent_power, proc_dp_phase_shift,
        reflectivity_horizontal fields.
    z_offset : float
        Horizontal reflectivity offset in dBZ.
    debug : bool
        True to print debugging information, False supressed this printing.

    Returns
    -------
    spec_at : dict
        Field dictionary containing the specific attenuation.
    cor_z : dict
        Field dictionary containing the corrected reflectivity.

    Other Parameters
    ----------------

    doc : float
        Number of gates at the end of each ray to to remove from the
        calculation.
    fzl : float
        Freezing layer, gates above this point are not included in the
        correction.
    rhv_min : float
        Minimum copol_coeff value to consider valid.
    ncp_min : float
        Minimum norm_coherent_power to consider valid.
    a_coef : float
        A coefficient in attenuation calculation.
    beta : float
        Beta parameter in attenuation calculation.

    References
    ----------
    Gu et al. Polarimetric Attenuation Correction in Heavy Rain at C Band,
    JAMC, 2011, 50, 39-58.

    """
    # extract fields and parameters from radar
    norm_coherent_power = radar.fields['norm_coherent_power']['data']
    copol_coeff = radar.fields['copol_coeff']['data']
    reflectivity_horizontal = radar.fields['reflectivity_horizontal']['data']
    proc_dp_phase_shift = radar.fields['proc_dp_phase_shift']['data']
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

    for sweep in xrange(nsweeps):
        # loop over the sweeps
        if debug:
            print "Doing ", sweep
        end_gate, start_ray, end_ray = phase_proc.det_process_range(
            radar, sweep, fzl, doc=doc)

        for i in xrange(start_ray, end_ray):
            # perform attenuation calculation on a single ray

            # extract the ray's phase shift and init. refl. correction
            ray_phase_shift = proc_dp_phase_shift[i, 0:end_gate]
            ray_init_refl = init_refl_correct[i, 0:end_gate]

            # perform calculation
            last_six_good = np.where(is_good[i, 0:end_gate])[0][-6:]
            phidp_max = np.median(ray_phase_shift[last_six_good])
            sm_refl = phase_proc.smooth_and_trim(ray_init_refl, window_len=5)
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
    spec_at = dict()
    spec_at['data'] = specific_atten
    spec_at['valid_min'] = 0.0
    spec_at['valid_max'] = 1.0
    spec_at['standard_name'] = 'specific_attenuation'
    spec_at['long_name'] = 'specific_attenuation'
    spec_at['least_significant_digit'] = 4
    spec_at['units'] = 'dB/km'

    cor_z = copy.deepcopy(radar.fields['reflectivity_horizontal'])
    cor_z['data'] = atten + cor_z['data'] + z_offset
    cor_z['data'].mask = init_refl_correct.mask
    cor_z['least_significant_digit'] = 2
    cor_z['valid_max'] = 80.0

    return spec_at, cor_z
