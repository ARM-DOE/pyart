"""
pyart.retrieve.kdp_proc
=======================

Module for retrieving specific differential phase (KDP) from radar total
differential phase (PSIDP) measurements. Total differential phase is a function
of propagation differential phase (PHIDP), backscatter differential phase
(DELTAHV), and the system phase offset.

.. autosummary::
    :toctree: generated/

    kdp_schneebeli
    kdp_vulpiani
    kdp_maesaka
    filter_psidp
    boundary_conditions_maesaka

    _kdp_estimation_backward_fixed
    _kdp_kalman_profile
    _kdp_vulpiani_profile
    _cost_maesaka
    _jac_maesaka
    _forward_reverse_phidp
    _parse_range_resolution

"""

from functools import partial
import time
import warnings

import numpy as np
from scipy import optimize, stats, interpolate, linalg, signal

from . import _kdp_proc
from ..config import get_field_name, get_metadata, get_fillvalue
from ..util import rolling_window


# Constants in the Kalman filter retrieval method (generally no need to
# modify them)
SCALERS = [0.1, 10**(-0.8), 10**(-0.6), 10**(-0.4), 10**(-0.2), 1,
           10**(0.2), 10**(0.4), 10**(0.6), 10**(0.8), 1, 10]

PADDING = 50  # Noise padding of the psidp signal (before and after signal)
SHIFT = 13  # Shifting of the final signal


def kdp_schneebeli(radar, gatefilter=None, fill_value=None, psidp_field=None,
                   kdp_field=None, phidp_field=None, band='C', rcov=0, pcov=0,
                   prefilter_psidp=False, filter_opt=None, parallel=True):
    """
    Estimates Kdp with the Kalman filter method by Schneebeli and al. (2014)
    for a set of psidp measurements.

    Parameters
    ----------
    radar : Radar
        Radar containing differential phase field.
    gatefilter : GateFilter, optional
        A GateFilter indicating radar gates that should be excluded when
        analysing differential phase measurements.
    fill_value : float, optional
        Value indicating missing or bad data in differential phase field, if
        not specified, the default in the Py-ART configuration file will be
        used
    psidp_field : str, optional
        Total differential phase field. If None, the default field name must be
        specified in the Py-ART configuration file.
    kdp_field : str, optional
        Specific differential phase field. If None, the default field name must
        be specified in the Py-ART configuration file.
    phidp_field : str, optional
        Propagation differential phase field. If None, the default field name
        must be specified in the Py-ART configuration file.
    band : char, optional
       Radar frequency band string. Accepted "X", "C", "S" (capital
       or not). The band is used to compute intercepts -c and slope b of the
       delta = b*Kdp+c relation
    rcov : 3x3 float array, optional
        Measurement error covariance matrix
    pcov  : 4x4 float array, optional
        Scaled state transition error covariance matrix
    prefilter_psidp : bool, optional
        If set, the psidp measurements will first be filtered with the
        filter_psidp method, which can improve the quality of the final Kdp
    filter_opt : dict, optional
        The arguments for the prefilter_psidp method, if empty, the defaults
        arguments of this method will be used
    parallel : bool, optional
        Flag to enable parallel computation (one core for every psidp profile)

    Returns
    -------
    kdp_dict : dict
        Retrieved specific differential phase data and metadata.
    kdp_std_dict : dict
        Estimated specific differential phase standard dev. data and metadata.
    phidpr_dict,: dict
        Retrieved differential phase data and metadata.

    References
    ----------
    Schneebeli, M., Grazioli, J., and Berne, A.: Improved Estimation
    of the Specific Differential Phase SHIFT Using a Compilation of
    Kalman Filter Ensembles, IEEE T. Geosci. Remote Sens., 52,
    5137-5149, doi:10.1109/TGRS.2013.2287017, 2014.

    """

    # create parallel computing instance
    if parallel:
        import multiprocessing as mp

        pool = mp.Pool(processes=mp.cpu_count(), maxtasksperchild=1)

    # parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # parse field names
    if psidp_field is None:
        psidp_field = get_field_name('differential_phase')
    if kdp_field is None:
        kdp_field = get_field_name('specific_differential_phase')
    if phidp_field is None:
        phidp_field = get_field_name('differential_phase')

    # parse range resolution, length scale and low-pass filter constraint
    # weight
    dr = _parse_range_resolution(radar, check_uniform=True)

    # parse total differential phase measurements
    if prefilter_psidp:
        if filter_opt is None:
            filter_opt = {}
        # Assign psidp field to filter_psidp inputs
        filter_opt['psidp_field'] = psidp_field
        # Filter psidp
        psidp_o = filter_psidp(radar, **filter_opt)
    else:
        psidp_o = radar.fields[psidp_field]['data']

    # mask radar gates indicated by the gate filter
    if gatefilter is not None:
        psidp_o = np.ma.masked_where(gatefilter.gate_excluded, psidp_o)

    func = partial(_kdp_kalman_profile, dr=dr, band=band, rcov=rcov,
                   pcov=pcov)

    all_psidp_prof = list(psidp_o)

    if parallel:
        list_est = pool.map(func, all_psidp_prof)
    else:
        list_est = map(func, all_psidp_prof)

    kdp = np.zeros(psidp_o.shape) * np.nan
    kdp = np.ma.masked_array(kdp, fill_value=fill_value)

    kdp_stdev = np.zeros(psidp_o.shape) * np.nan
    kdp_stdev = np.ma.masked_array(kdp_stdev, fill_value=fill_value)

    phidp_rec = np.zeros(psidp_o.shape) * np.nan
    phidp_rec = np.ma.masked_array(phidp_rec, fill_value=fill_value)

    for i, l in enumerate(list_est):
        kdp[i, 0:len(l[0])] = l[0]
        kdp_stdev[i, 0:len(l[1])] = l[1]
        phidp_rec[i, 0:len(l[2])] = l[2]

    # Mask the estimated Kdp and reconstructed Phidp with the mask of original
    # psidp
    if isinstance(psidp_o, np.ma.masked_array):
        masked = psidp_o.mask
        kdp = np.ma.array(kdp, mask=masked, fill_value=fill_value)
        kdp_stdev = np.ma.array(kdp_stdev, mask=masked, fill_value=fill_value)
        phidp_rec = np.ma.array(phidp_rec, mask=masked, fill_value=fill_value)

    # create specific differential phase field dictionary and store data
    kdp_dict = get_metadata(kdp_field)
    kdp_dict['data'] = kdp
    # kdp_dict['valid_min'] = -1.0

    # create reconstructed differential phase field dictionary and store data
    phidpr_dict = get_metadata(phidp_field)
    phidpr_dict['data'] = phidp_rec
    # phidpr_dict['valid_min'] = 0.0

    # create specific phase stdev field dictionary and store data
    kdp_stdev_dict = {}
    kdp_stdev_dict['units'] = 'degrees/km'
    kdp_stdev_dict['standard_name'] = 'estimated KDP stdev'
    kdp_stdev_dict['long_name'] = 'Estimated stdev of spec. diff. phase (KDP)'
    kdp_stdev_dict['coordinates'] = 'elevation azimuth range'
    kdp_stdev_dict['data'] = kdp_stdev
    kdp_stdev_dict['valid_min'] = 0.0

    if parallel:
        pool.close()

    return kdp_dict, kdp_stdev_dict, phidpr_dict


def _kdp_estimation_backward_fixed(
        psidp_in, rcov, pcov_scale, f, f_transposed, h_plus,
        c1, c2, b1, b2, kdp_th, mpsidp):
    """
    Processing one profile of Psidp and estimating Kdp and Phidp
    with the KFE algorithm described in Schneebeli et al, 2014
    IEEE_TGRS. This routine estimates Kdp in the backward
    direction given a set of matrices that define the Kalman
    filter.

    Parameters
    ----------
    psidp_in : ndarray
        one-dimensional vector of length -nrg- containining the input psidp
        [degrees]
    rcov : 3x3 float array##
        Measurement error covariance matrix
    pcov_scale  : 4x4 float array
        Scaled state transition error covariance matrix
    f : 4x4 float array
        Forward state prediction matrix [4x4]
    f_transposed: 4x4 float array
        Transpose of F
    h_plus : 4x3 float array
        Measurement prediction matrix [4x3]
    c1, c2,b1,b2: floats
        the values of the intercept of the relation c  = b*Kdp - delta.
        This relation uses b1, c1 IF kdp is lower than a kdp_th and b2, c2
        otherwise kdp_th
    kdp_th: float
        the kdp threshold which separates the two Kdp - delta regime
        i.e. the power law relating delta to Kdp will be different if Kdp is
        larger or smaller than kdp_th
    mpsidp: float
        final observed value of psidp along the radial (usually also
        the max value), needed for inverting the psidp vector


    Returns
    -------
    kdp: ndarray
        filtered Kdp [degrees/km]. Same length as Psidp
    error_kdp: ndarray
        estimated error on Kdp values

    """

    # Define the input
    psidp = psidp_in

    # invert Psidp (backward estimation)
    psidp = mpsidp - psidp[::-1]
    nrg_new = len(psidp)

    # Initialize the state vector to 0
    s = np.zeros([4, 1])  # first state estimate

    # define measurement vector
    z = np.zeros([3, 1])

    # Define the identity matrix
    identity_i = np.eye(4)
    p = identity_i * 4.

    kdp = np.zeros([nrg_new])
    kdp_error = np.zeros([nrg_new])

    # Loop on all the gates and apply the filter

    for ii in range(0, nrg_new - 1):
        z[0] = psidp[ii]
        z[1] = psidp[ii + 1]

        s_pred = np.dot(f, s)  # state prediciton

        p_pred = np.dot(f, np.dot(p, f_transposed)) + \
            pcov_scale  # error prediction

        if s_pred[0] > kdp_th:
            h_plus[2, 0] = b2
            z[2] = c2
        else:
            h_plus[2, 0] = b1
            z[2] = c1

        # as far as i see aludc is symmetrical, so i do not transpose it
        aludc = (np.dot(h_plus, (np.dot(p_pred, h_plus.T))) + rcov)

        # below we get the transposed of B_mat directly
        b_mat = np.dot(h_plus, p_pred)

        # Cholesky decomposition
        cho = linalg.cho_factor(aludc)
        k = linalg.cho_solve(cho, b_mat, check_finite=False,
                             overwrite_b=True).T

        # Update state and error
        s = np.dot(k, (np.dot(-h_plus, s_pred) + z)) + s_pred
        p = np.dot((identity_i - np.dot(k, h_plus)), p_pred)

        # Fill the output
        kdp[ii] = s[0]
        kdp_error[ii] = p[0, 0]

    # Shift
    n = len(kdp)
    dummy = np.copy(kdp)
    kdp[np.arange(SHIFT) + len(kdp) - SHIFT] = 0
    kdp[np.arange(n - 1 - SHIFT)] = dummy[np.arange(n - 1 - SHIFT) + SHIFT]

    # Reverse the estimates (backward direction)
    kdp = kdp[::-1]

    kdp_error = kdp_error[::-1]
    return kdp, kdp_error


def _kdp_estimation_forward_fixed(
        psidp_in, rcov, pcov_scale, f, f_transposed, h_plus,
        c1, c2, b1, b2, kdp_th):
    """
    Processing one profile of Psidp and estimating Kdp and Phidp
    with the KFE algorithm described in Schneebeli et al, 2014
    IEEE_TGRS. This routine estimates Kdp in the forward
    direction given a set of matrices that define the Kalman
    filter.

    Parameters
    ----------
    psidp_in : ndarray
        one-dimensional vector of length -nrg- containining the input psidp
        [degrees]
    rcov : 3x3 float array
        Measurement error covariance matrix
    pcov_scale  : 4x4 float array
        Scaled state transition error covariance matrix
    f : 4x4 float array
        Forward state prediction matrix [4x4]
    f_transposed: 4x4 float array
        Transpose of F
    h_plus : 4x3 float array*np.nan
        Measurement prediction matrix [4x3]
    c1, c2,b1,b2: floats
        the values of the intercept of the relation c  = b*Kdp - delta.
        This relation uses b1, c1 IF kdp is lower than a kdp_th and b2, c2
        otherwise kdp_th.

    Returns
    -------
    kdp: ndarray
        filtered Kdp [degrees/km]. Same length as Psidp
    phidp: ndarray
        estimated phidp (smooth psidp)
    error_kdp: ndarray
        estimated error on Kdp values

    """

    # Define the input
    psidp = np.ma.filled(psidp_in)
    nrg_new = len(psidp)

    # Initialize the state vector to 0
    s = np.zeros([4, 1])  # first state estimate

    # define measurement vector
    z = np.zeros([3, 1])

    # Define the identity matrix
    identity_i = np.eye(4)
    p = identity_i * 4.

    phidp = np.zeros([nrg_new])
    kdp = np.zeros([nrg_new])
    kdp_error = np.zeros([nrg_new])

    # Loop on all the gates and apply the filter
    for ii in range(0, nrg_new - 1):
        z[0] = psidp[ii]
        z[1] = psidp[ii + 1]

        s_pred = np.dot(f, s)  # state prediciton

        p_pred = np.dot(f, np.dot(p, f_transposed)) + \
            pcov_scale  # error prediction

        if s_pred[0] > kdp_th:
            h_plus[2, 0] = b2
            z[2] = c2
        else:
            h_plus[2, 0] = b1
            z[2] = c1

        # as far as i see aludc is symmetrical, so i do not transpose it
        aludc = (np.dot(h_plus, (np.dot(p_pred, h_plus.T))) + rcov)

        # below we get the transposed of B_mat directly
        b_mat = np.dot(h_plus, p_pred)

        # Cholesky decomposition
        cho = linalg.cho_factor(aludc)
        k = linalg.cho_solve(
            cho,
            b_mat,
            check_finite=False,
            overwrite_b=True).T

        # Update state and error
        s = np.dot(k, (np.dot(-h_plus, s_pred) + z)) + s_pred
        p = np.dot((identity_i - np.dot(k, h_plus)), p_pred)

        # Fill the output
        kdp[ii] = s[0]
        kdp_error[ii] = p[0, 0]
        phidp[ii] = s[2]

    # Shift
    dummy = np.copy(kdp)
    kdp[np.arange(SHIFT) + len(kdp) - SHIFT] = 0
    kdp[np.arange(len(kdp) - SHIFT)
        ] = dummy[np.arange(len(kdp) - SHIFT) + SHIFT]

    return kdp, phidp, kdp_error


def _kdp_kalman_profile(psidp_in, dr, band='X', rcov=0, pcov=0):
    """
    Estimates Kdp with the Kalman filter method by Schneebeli and al. (2014)
    for a set of psidp measurements.

    Parameters
    ----------
    psidp_in : ndarray
        one-dimensional vector of length -nrg- containining the input psidp
        [degrees]
    dr : float
        Range resolution in meters.
    band : char, optional
       Radar frequency band string. Accepted "X", "C", "S" (capital
       or not). The band is used to compute intercepts -c and slope b of the
       delta = b*Kdp+c relation
    rcov : 3x3 float array, optional
        Measurement error covariance matrix
    pcov  : 4x4 float array, optional
        Scaled state transition error covariance matrix

    Returns
    -------
    kdp_dict : ndarray
        Retrieved specific differential phase data
    kdp_std_dict : ndarray
        Estimated specific differential phase standard dev. data
    phidpr_dict,: ndarray
        Retrieved differential phase data

    References
    ----------
    Schneebeli, M., Grazioli, J., and Berne, A.: Improved Estimation
    of the Specific Differential Phase Shift Using a Compilation of
    Kalman Filter Ensembles, IEEE T. Geosci. Remote Sens., 52,
    5137-5149, doi:10.1109/TGRS.2013.2287017, 2014.

    """

    dr = dr / 1000.  # Convert rad. res. to km

    # NOTE! Parameters are not checked to save as much time as possible

    # Replace missing values with nans
    psidp_in = np.ma.filled(psidp_in, np.nan)
    # Check if psidp has at least one finite value
    if not np.isfinite(psidp_in).any():
        return psidp_in, psidp_in, psidp_in  # Return the NaNs...

    # Set default of the error covariance matrices
    if not isinstance(pcov, np.ndarray):
        pcov = np.array([[(0.11 + 1.56 * dr)**2,
                          (0.11 + 1.85 * dr)**2,
                          0,
                          (0.01 + 1.1 * dr)**2],
                         [(0.11 + 1.85 * dr)**2,
                          (0.18 + 3.03 * dr)**2,
                          0,
                          (0.01 + 1.23 * dr)**2],
                         [0,
                          0,
                          0,
                          0],
                         [(0.01 + 1.1 * dr)**2,
                          (0.01 + 1.23 * dr)**2,
                          0,
                          (-0.04 + 1.27 * dr)**2]])

    if not isinstance(rcov, np.ndarray):
        rcov = np.array([[4.10625, -0.0498779, -0.0634192],
                         [-0.0498779, 4.02369, -0.0421455],
                         [-0.0634192, -0.0421455, 1.44300]])
    '''
    Define default parameters
    '''
    # Intercepts -c and slope b of delta=b*Kdp+c
    # According to the Kdp threshold selected

    if band == 'X':
        c1 = -0.054
        c2 = -6.155
        b1 = 2.3688
        b2 = 0.2734
        kdp_th = 2.5
    elif band == 'C':
        c1 = -0.036
        c2 = -1.03
        b1 = 0.53
        b2 = 0.15
        kdp_th = 2.5
    elif band == 'S':
        c1 = -0.024
        c2 = -0.15
        b1 = 0.19
        b2 = 0.019
        kdp_th = 1.1

    # Parameters for the final selection from the KF ensemble members
    fac1 = 1.2
    fac2 = 3.

    th1_comp = -0.15
    th2_comp = 0.15
    th1_final = -0.25

    # Kalman matrices
    # State matrix
    f = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1],
                  [2 * dr, 0, 0, 1]], dtype=float)
    f_transposed = f.T

    # Measurement prediction matrix--------------------------------
    # J. Grazioli modification 07.2015 --previous H_plus buggy--
    h_plus = np.array(
        [[-2 * dr, 1, 0, 1], [2 * dr, 1, 1, 0], [0, -1, 0, 0]], dtype=float)

    # Define the input
    psidp = psidp_in
    # Get indices of finite data
    real_data_ind = np.where(np.isfinite(psidp.ravel()))[0]
    offset = real_data_ind[0]

    if len(real_data_ind):
        mpsidp = psidp.ravel()[real_data_ind[-1]]
    else:
        mpsidp = np.nan

    psidp = psidp[offset:real_data_ind[-1] + 1]

    nrg = len(psidp)

    # Define the output
    kdp_filter_out = np.zeros([nrg, ])

    kdp_mat = np.zeros([nrg, 2 * len(SCALERS)])
    kdp_sim = np.zeros([nrg, len(SCALERS)])
    phidp_filter_out = np.zeros([nrg])

    '''
    Prepare a longer array with some extra gates on each side
    '''
    # add  values at the beginning and at the end of the profile
    psidp_long = np.zeros([nrg + PADDING * 2, ]) * np.nan

    nn = nrg + PADDING * 2

    noise = 2 * np.random.randn(PADDING)
    psidp_long[0:PADDING] = noise + psidp[0]
    psidp_long[nrg + PADDING: nrg + 2 * PADDING] = mpsidp + noise
    psidp_long[PADDING:nrg + PADDING] = psidp

    psidp = psidp_long

    # Get information of valid and non valid points in psidp the new psidp
    nonan = np.where(np.isfinite(psidp))[0]
    nan = np.where(np.isnan(psidp))[0]

    ranged = np.arange(0, nn)

    psidp_interp = psidp

    # interpolate
    if len(nan):
        interp = interpolate.interp1d(ranged[nonan], psidp[nonan], kind='zero')
        psidp_interp[nan] = interp(ranged[nan])
    else:
        psidp_interp = psidp

    # add noise
    if len(nan):
        psidp_interp[nan] = psidp_interp[nan] + 2 * np.random.randn(len(nan))

    # Define the final input and output
    psidp = psidp_interp

    '''
    Smallest scaler
    '''
    scaler = 10 ** (-2.)

    # Backward
    kdp_dummy_b2, _ = _kdp_estimation_backward_fixed(psidp, rcov,
                                                     pcov * scaler, f,
                                                     f_transposed, h_plus,
                                                     c1, c2, b1, b2,
                                                     kdp_th, mpsidp)
    kdp002 = kdp_dummy_b2[PADDING:nrg + PADDING]

    # Forward
    kdp_dummy_b2, _, _ = _kdp_estimation_forward_fixed(psidp, rcov,
                                                       pcov * scaler, f,
                                                       f_transposed, h_plus,
                                                       c1, c2, b1, b2, kdp_th)
    kdp002f = kdp_dummy_b2[PADDING:nrg + PADDING]

    '''
    Generate the ensemble of Kalman filters estimates in backward and
    Forward directions
    '''
    for i, sc in enumerate(SCALERS):  # Loop on scalers
        # Forward
        kdp_dummy_f2, _, _ = _kdp_estimation_forward_fixed(psidp,
                                                           rcov, pcov * sc,
                                                           f, f_transposed,
                                                           h_plus, c1, c2,
                                                           b1, b2, kdp_th)
        kdp_mat[:, 2 * i] = kdp_dummy_f2[PADDING:nrg + PADDING]

        # Backward
        kdp_dummy_b2, _ = _kdp_estimation_backward_fixed(psidp, rcov,
                                                         pcov * sc, f,
                                                         f_transposed, h_plus,
                                                         c1, c2, b1, b2,
                                                         kdp_th, mpsidp)
        kdp_mat[:, 2 * i + 1] = kdp_dummy_b2[PADDING:nrg + PADDING]

    '''
    Compile the final estimate
    '''
    # Get some reference mean values
    kdp_mean = np.nanmean(kdp_mat, axis=1)
    kdp_mean_shift = np.roll(kdp_mean, -1)
    diff_mean = kdp_mean - kdp_mean_shift

    kdp_std = np.nanstd(kdp_mat, axis=1)

    if len(diff_mean) < 4:
        size_filt = len(diff_mean)
    else:
        size_filt = 4

    diff_mean_smooth = np.convolve(diff_mean,
                                   np.ones((size_filt,)) / size_filt,
                                   mode='same')

    # Backward estimate if diff_mean greater than a defined threshold
    condi = np.where(diff_mean_smooth > th2_comp)[0]
    if len(condi):
        kdp_dummy = kdp_mat[:, np.arange(len(SCALERS)) * 2 + 1]
        kdp_sim[condi, :] = kdp_dummy[condi, :]

    # Forward estimate if diff_mean lower than a defined threshold
    condi = np.where(diff_mean_smooth < th1_comp)[0]
    if len(condi):
        kdp_dummy = kdp_mat[:, np.arange(len(SCALERS)) * 2]
        kdp_sim[condi, :] = kdp_dummy[condi, :]

    # Combination of the two in the middle
    condi = np.where(np.logical_and(diff_mean_smooth >= th1_comp,
                                    diff_mean_smooth <= th2_comp))[0]
    if len(condi):
        weight2 = (-0.5 / 0.15) * diff_mean_smooth + 0.5

        weight2 = np.tile(weight2, (len(SCALERS), 1)).T

        kdp_dummy = (1-weight2)*kdp_mat[:, np.arange(len(SCALERS))*2+1] \
            + weight2 * kdp_mat[:, np.arange(len(SCALERS)) * 2]
        kdp_sim[condi, :] = kdp_dummy[condi, :]

    # Now we reduced to 11 ensemble members: compile the final one
    kdp_mean_sim = np.nanmean(kdp_sim, axis=1)
    kdp_std_sim = np.nanstd(kdp_sim, axis=1)
    kdp_low_mean2 = np.nanmean(np.vstack((kdp002, kdp002f)).T, axis=1)

    # Get the range of ensemble members that compile
    # a final estimate

    # Lower bounds
    lower_bound = np.round(kdp_mean_sim * fac1) - np.round(kdp_std_sim * fac2)
    lower_bound = np.maximum(lower_bound, 0)
    lower_bound = np.minimum(lower_bound, len(SCALERS) - 1)

    # Upper bounds
    upper_bound = np.round(kdp_mean_sim * fac1) + np.round(kdp_std_sim * fac2)
    upper_bound = np.maximum(upper_bound, 0)
    upper_bound = np.minimum(upper_bound, len(SCALERS) - 1)

    # Final selection of the ensemble members
    for uu in range(0, nrg - 1):
        selection_vector = np.arange(
            upper_bound[uu] - lower_bound[uu] + 1) + lower_bound[uu]
        selection_vector = selection_vector.astype(int)
        kdp_filter_out[uu] = np.mean(kdp_sim[uu, selection_vector])

    # Final filtering of excessively negative values:
    # TO DO: It would be better to get rid of this filtering
    ind_lt_0 = np.where(kdp_filter_out < th1_final)[0]

    if len(ind_lt_0):
        kdp_filter_out[ind_lt_0] = kdp_low_mean2[ind_lt_0]

    # Compute phidp from Kdp
    phidp_filter_out = np.cumsum(kdp_filter_out) * 2. * dr

    phinan = np.where(np.isnan(psidp))[0]
    if len(phinan):
        phidp_filter_out[phinan] = np.nan
        kdp_filter_out[phinan] = np.nan
    phidp_filter_out[nrg - 1] = np.nan
    kdp_filter_out[nrg - 1] = np.nan

    # Pad with nan if offset > 0
    kdp_filter_out = np.pad(kdp_filter_out, (offset, 0), mode='constant',
                            constant_values=np.nan)
    kdp_std = np.pad(kdp_std, (offset, 0), mode='constant',
                     constant_values=np.nan)
    phidp_filter_out = np.pad(phidp_filter_out, (offset, 0), mode='constant',
                              constant_values=np.nan)

    return kdp_filter_out, kdp_std, phidp_filter_out


def kdp_vulpiani(radar, gatefilter=None, fill_value=None, psidp_field=None,
                 kdp_field=None, phidp_field=None, band='C', windsize=10,
                 n_iter=10, interp=False, prefilter_psidp=False,
                 filter_opt=None, parallel=False):
    """
    Estimates Kdp with the Vulpiani method for a 2D array of psidp measurements
    with the first dimension being the distance from radar and the second
    dimension being the angles (azimuths for PPI, elev for RHI).The input psidp
    is assumed to be pre-filtered (for ex. with the filter_psidp function)

    Parameters
    ----------
    radar : Radar
        Radar containing differential phase field.
    gatefilter : GateFilter, optional
        A GateFilter indicating radar gates that should be excluded when
        analysing differential phase measurements.
    fill_value : float, optional
        Value indicating missing or bad data in differential phase field, if
        not specified, the default in the Py-ART configuration file will be
        used
    psidp_field : str, optional
        Total differential phase field. If None, the default field name must be
        specified in the Py-ART configuration file.
    kdp_field : str, optional
        Specific differential phase field. If None, the default field name must
        be specified in the Py-ART configuration file.
    phidp_field : str, optional
        Propagation differential phase field. If None, the default field name
        must be specified in the Py-ART configuration file.
    band : char, optional
        Radar frequency band string. Accepted "X", "C", "S" (capital
        or not). It is used to set default boundaries for expected
        values of Kdp.
    windsize : int, optional
        Size in # of gates of the range derivative window. Should be even.
    n_iter : int, optional
        Number of iterations of the method. Default is 10.
    interp : bool, optional
        If True, all the nans are interpolated.The advantage is that less data
        are lost (the iterations in fact are "eating the edges") but some
        non-linear errors may be introduced.
    prefilter_psidp : bool, optional
        If set, the psidp measurements will first be filtered with the
        filter_psidp method, which can improve the quality of the final Kdp.
    filter_opt : dict, optional
        The arguments for the prefilter_psidp method, if empty, the defaults
        arguments of this method will be used.
    parallel : bool, optional
        Flag to enable parallel computation (one core for every psidp profile).

    Returns
    -------
    kdp_dict : dict
        Retrieved specific differential phase data and metadata.
    phidpr_dict,: dict
        Retrieved differential phase data and metadata.

    References
    ----------
    Gianfranco Vulpiani, Mario Montopoli, Luca Delli Passeri, Antonio G. Gioia,
    Pietro Giordano, and Frank S. Marzano, 2012: On the Use of Dual-Polarized
    C-Band Radar for Operational Rainfall Retrieval in Mountainous Areas.
    J. Appl. Meteor. Climatol., 51, 405-425, doi: 10.1175/JAMC-D-10-05024.1.

    """
    if np.mod(windsize, 2):
        warnings.warn('In the Vulpiani method, the windsize should be even. '
                      + 'Using default value, windsize = 10')
        windsize = 10

    if parallel:
        import multiprocessing as mp
        pool = mp.Pool(processes=mp.cpu_count(), maxtasksperchild=1)

    # parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # parse field names
    if psidp_field is None:
        psidp_field = get_field_name('differential_phase')
    if kdp_field is None:
        kdp_field = get_field_name('specific_differential_phase')
    if phidp_field is None:
        phidp_field = get_field_name('differential_phase')

    # parse range resolution, length scale and low-pass filter constraint
    # weight
    dr = _parse_range_resolution(radar, check_uniform=True)

    # parse total differential phase measurements
    if prefilter_psidp:
        if filter_opt is None:
            filter_opt = {}
        # Assign psidp field to filter_psidp inputs
        filter_opt['psidp_field'] = psidp_field
        # Filter psidp
        psidp_o = filter_psidp(radar, **filter_opt)
    else:
        psidp_o = radar.fields[psidp_field]['data']

    # mask radar gates indicated by the gate filter
    if gatefilter is not None:
        psidp_o = np.ma.masked_where(gatefilter.gate_excluded, psidp_o)

    func = partial(_kdp_vulpiani_profile, dr=dr, windsize=windsize,
                   band=band, n_iter=n_iter, interp=interp)

    all_psidp_prof = list(psidp_o)

    if parallel:
        list_est = pool.map(func, all_psidp_prof)
    else:
        list_est = map(func, all_psidp_prof)

    kdp = np.ma.zeros(psidp_o.shape)
    kdp[:] = np.ma.masked
    kdp.set_fill_value(fill_value)

    phidp_rec = np.ma.zeros(psidp_o.shape)
    phidp_rec[:] = np.ma.masked
    phidp_rec.set_fill_value(fill_value)

    for i, l in enumerate(list_est):
        kdp[i, 0:len(l[0])] = l[0]
        phidp_rec[i, 0:len(l[1])] = l[1]

    # Mask the estimated Kdp and reconstructed Phidp with the mask of original
    # psidp
    if isinstance(psidp_o, np.ma.masked_array):
        masked = psidp_o.mask
        kdp = np.ma.array(kdp, mask=masked, fill_value=fill_value)
        phidp_rec = np.ma.array(phidp_rec, mask=masked, fill_value=fill_value)

    # create specific differential phase field dictionary and store data
    kdp_dict = get_metadata(kdp_field)
    kdp_dict['data'] = kdp
    # kdp_dict['valid_min'] = 0.0

    # create reconstructed differential phase field dictionary and store data
    phidpr_dict = get_metadata(phidp_field)
    phidpr_dict['data'] = phidp_rec
    # phidpr_dict['valid_min'] = 0.0

    if parallel:
        pool.close()

    return kdp_dict, phidpr_dict


def _kdp_vulpiani_profile(psidp_in, dr, windsize=10,
                          band='X', n_iter=10, interp=False):
    """
    Estimates Kdp with the Vulpiani method for a single profile of psidp
    measurements

    Parameters
    ----------
    psidp_in : ndarray
        Total differential phase measurements.
    dr : float
        Range resolution in meters.
    windsize : int, optional
        Size in # of gates of the range derivative window.
    band : char, optional
        Radar frequency band string. Accepted "X", "C", "S" (capital
        or not). It is used to set default boundaries for expected
        values of Kdp
    n_iter : int, optional
        Number of iterations of the method. Default is 10.
    interp : bool, optional
        If set all the nans are interpolated.The advantage is that less data
        are lost (the iterations in fact are "eating the edges") but some
        non-linear errors may be introduced

    Returns
    -------
    kdp_calc : ndarray
        Retrieved specific differential profile
    phidp_rec,: ndarray
        Retrieved differential phase profile

    """
    mask = np.ma.getmaskarray(psidp_in)
    l = windsize
    l2 = int(l/2)
    drm = dr/1000.

    if mask.all() is True:
        # Check if all elements are masked
        return psidp_in, psidp_in, psidp_in  # Return the NaNs...

    # Thresholds in kdp calculation
    if band == 'X':
        th1 = -2.
        th2 = 40.
        std_th = 5.
    elif band == 'C':
        th1 = -2.
        th2 = 20.
        std_th = 5.
    elif band == 'S':
        th1 = -2.
        th2 = 14.
        std_th = 5.
    else:
        print('Unexpected value set for the band keyword ')
        print(band)
        return None

    psidp = psidp_in
    nn = len(psidp_in)

    # Get information of valid and non valid points in psidp the new psidp
    valid = np.logical_not(mask)
    if interp:
        ranged = np.arange(0, nn)
        psidp_interp = psidp
        # interpolate
        if np.ma.is_masked(psidp):
            interp = interpolate.interp1d(ranged[valid], psidp[valid],
                                          kind='zero', bounds_error=False,
                                          fill_value=np.nan)
            psidp_interp[mask] = interp(ranged[mask])

        psidp = psidp_interp

    psidp = np.ma.filled(psidp, np.nan)
    kdp_calc = np.zeros([nn])

    # first guess
    # In the core of the profile
    kdp_calc[l2:nn-l2] = (psidp[l:nn]-psidp[0:nn-l])/(2.*l*drm)

    # set ray extremes to 0
    kdp_calc[0:l2] = 0.
    kdp_calc[nn-l2:] = 0.

    # apply thresholds
    kdp_calc[kdp_calc <= th1] = 0.
    kdp_calc[kdp_calc >= th2] = 0.

    # set all non-valid data to 0
    kdp_calc[np.isnan(kdp_calc)] = 0.

    # Remove bins with texture higher than treshold
    tex = np.ma.zeros(kdp_calc.shape)
    # compute the local standard deviation
    # (make sure that it is and odd window)
    tex_aux = np.ma.std(rolling_window(kdp_calc, l2*2+1), -1)
    tex[l2:-l2] = tex_aux
    kdp_calc[tex > std_th] = 0.

    # Loop over iterations
    for i in range(0, n_iter):
        phidp_rec = np.ma.cumsum(kdp_calc)*2.*drm

        # In the core of the profile
        kdp_calc[l2:nn-l2] = (phidp_rec[l:nn]-phidp_rec[0:nn-l])/(2.*l*drm)

        # set ray extremes to 0
        kdp_calc[0:l2] = 0.
        kdp_calc[nn-l2:] = 0.

        # apply thresholds
        kdp_calc[kdp_calc <= th1] = 0.
        kdp_calc[kdp_calc >= th2] = 0.

    # Censor Kdp where Psidp was not defined
    kdp_calc = np.ma.masked_where(mask, kdp_calc)

    # final reconstructed PhiDP from KDP
    phidp_rec = np.ma.cumsum(kdp_calc)*2.*drm

    return kdp_calc, phidp_rec


def filter_psidp(radar, psidp_field=None, rhohv_field=None, minsize_seq=5,
                 median_filter_size=7, thresh_rhohv=0.65, max_discont=90):
    """
    Filter measured psidp to remove spurious data in four steps:
         1. Censor it where Rhohv is lower than threshold
         2. Unravel angles when strong discontinuities are detected
         3. Remove very short sequences of valid data
         4. Apply a median filter on every profile

    Parameters
    ----------
    radar : Radar
        Radar containing differential phase field.
    psidp_field : str, optional
        Total differential phase field. If None, the default field name must be
        specified in the Py-ART configuration file.
    rhohv_field : str, optional
        Cross correlation ratio field. If None, the default field name must
        be specified in the Py-ART configuration file.
    minsize_seq  : integer, optional
        Minimal len (in radar gates) of sequences of valid data to be accepted
    median_filter_size : integer, optional
        Size (in radar gates) of the median filter to be applied on psidp
    thresh_rhohv : float, optional
        Censoring threshold in rhohv (gates with rhohv < thresh_rhohv)
        will be rejected
    max_discont : int, optional
        Maximum discontinuity between psidp values, default is 90 deg

    Returns
    -------
    psidp_filt : ndarray
        Filtered psidp field

    """
    # parse field names
    if psidp_field is None:
        psidp_field = get_field_name('differential_phase')
    if rhohv_field is None:
        rhohv_field = get_field_name('cross_correlation_ratio')

    # parse total differential phase and cross corr. measurements
    psidp_o = radar.fields[psidp_field]['data']
    rhohv = radar.fields[rhohv_field]['data']

    if not isinstance(psidp_o, np.ma.masked_array):
        psidp_o = np.ma.array(psidp_o, mask=np.isnan(psidp_o))

    # Initialize mask
    mask = np.ones(psidp_o.shape) * False

    # Condition on rhohv
    mask += rhohv < thresh_rhohv
    # Get original mask
    mask += psidp_o.mask

    # Remove short sequences and unwrap
    psidp_filt = np.zeros(psidp_o.shape)
    for i, psi_row in enumerate(psidp_o):
        idx = np.where(~psi_row.mask)[0]  # idx of last valid gate
        if len(idx):
            psi_row = psi_row[0:idx[-1] + 1]
            psi_row[~psi_row.mask] = np.rad2deg(
                np.unwrap(np.deg2rad(psi_row[~psi_row.mask]),
                          np.deg2rad(max_discont)))
            psi_row_with_nan = np.ma.filled(psi_row, np.nan)
            # To be sure to always have a left and right neighbour,
            # we need to pad the signal with NaN
            psi_row_with_nan = np.pad(
                psi_row_with_nan, (1, 1), 'constant', constant_values=(
                    np.nan,))
            idx = np.where(np.isfinite(psi_row_with_nan))[0]
            nan_left = idx[np.where(np.isnan(psi_row_with_nan[idx - 1]))[0]]
            nan_right = idx[np.where(np.isnan(psi_row_with_nan[idx + 1]))[0]]

            len_sub = nan_right - nan_left

            for j, l in enumerate(len_sub):
                if l < minsize_seq:
                    mask[i, nan_left[j] - 1:nan_right[j] + 1] = True

            # median filter
            psi_row = signal.medfilt(psi_row_with_nan, median_filter_size)
            psidp_filt[i, 0:len(psi_row[1:-1])] = psi_row[1:-1]

    psidp_filt = np.ma.masked_array(
        psidp_filt, mask=mask, fill_value=psidp_o.fill_value)

    return psidp_filt


# Necessary and/or potential future improvements to the KDP module:
#
# * The near and far range gate boundary conditions necessary for the Maesaka
#   et al. (2012) variational method are sensitive to the number of range gates
#   used to define both boundaries as well as outliers in these samples. This
#   may have to be further investigated and potentially improved upon.
#
# * It would be beneficial for the GateFilter class to have a method which can
#   flag radar gates above the melting layer, maybe through the user supplying
#   an estimate of the 0 deg Celsius level above mean sea level.
#
# * The backscatter differential phase parameter will likely have to be updated
#   in the future to handle various parameterizations.
#
# * The addition of an azimuthal low-pass filter (smoothness) constraint would
#   be beneficial to the variational Maesaka et al. (2012) algorithm.
#
# * Add the ability to record and/or store minimization data, e.g., current
#   iteration number, functional value and functional gradient norm value as a
#   function of iteration.


def kdp_maesaka(radar, gatefilter=None, method='cg', backscatter=None,
                Clpf=1.0, length_scale=None, first_guess=0.01,
                finite_order='low', fill_value=None, proc=1, psidp_field=None,
                kdp_field=None, phidp_field=None, debug=False, verbose=False,
                **kwargs):
    """
    Compute the specific differential phase (KDP) from corrected (e.g.,
    unfolded) total differential phase data based on the variational method
    outlined in Maesaka et al. (2012). This method assumes a monotonically
    increasing propagation differential phase (PHIDP) with increasing range
    from the radar, and therefore is limited to rainfall below the melting
    layer and/or warm clouds at weather radar frequencies (e.g., S-, C-, and
    X-band). This method currently only supports radar data with constant range
    resolution.

    Following the notation of Maesaka et al. (2012), the primary control
    variable k is proportional to KDP,

                                k**2 = 2 * KDP * dr

    which, because of the square, assumes that KDP always takes a positive
    value.

    Parameters
    ----------
    radar : Radar
        Radar containing differential phase field.
    gatefilter : GateFilter
        A GateFilter indicating radar gates that should be excluded when
        analysing differential phase measurements.
    method : str, optional
        Type of scipy.optimize method to use when minimizing the cost
        functional. The default method uses a nonlinear conjugate gradient
        algorithm. In Maesaka et al. (2012) they use the Broyden-Fletcher-
        Goldfarb-Shanno (BFGS) algorithm, however for large functional size
        (e.g., 100K+ variables) this algorithm is considerably slower than a
        conjugate gradient algorithm.
    backscatter : optional
        Define the backscatter differential phase. If None, the backscatter
        differential phase is set to zero for all range gates. Note that
        backscatter differential phase can be parameterized using attentuation
        corrected differential reflectivity.
    Clpf : float, optional
        The low-pass filter (radial smoothness) constraint weight as in
        equation (15) of Maesaka et al. (2012).
    length_scale : float, optional
        Length scale in meters used to bring the dimension and magnitude of the
        low-pass filter cost functional in line with the observation cost
        functional. If None, the length scale is set to the range resolution.
    first_guess : float, optional
        First guess for control variable k. Since k is proportional to the
        square root of KDP, the first guess should be close to zero to signify
        a KDP field close to 0 deg/km everywhere. However, the first guess
        should not be exactly zero in order to avoid convergence criteria after
        the first iteration. In fact it is recommended to use a value closer to
        one than zero.
    finite_order : 'low' or 'high', optional
        The finite difference accuracy to use when computing derivatives.
    maxiter : int, optional
        Maximum number of iterations to perform during cost functional
        minimization. The maximum number of iterations are only performed if
        convergence criteria are not met. For variational schemes such as this
        one, it is generally not recommended to try and achieve convergence
        criteria since the values of the cost functional and/or its gradient
        norm are somewhat arbitrary.
    fill_value : float, optional
        Value indicating missing or bad data in differential phase field.
    proc : int, optional
        The number of parallel threads (CPUs) to use. Currently no
        multiprocessing capability exists.
    psidp_field : str, optional
        Total differential phase field. If None, the default field name must be
        specified in the Py-ART configuration file.
    kdp_field : str, optional
        Specific differential phase field. If None, the default field name must
        be specified in the Py-ART configuration file.
    phidp_field : str, optional
        Propagation differential phase field. If None, the default field name
        must be specified in the Py-ART configuration file.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    Returns
    -------
    kdp_dict : dict
        Retrieved specific differential phase data and metadata.
    phidpf_dict, phidpr_dict : dict
        Retrieved forward and reverse direction propagation differential phase
        data and metadata.

    References
    ----------
    Maesaka, T., Iwanami, K. and Maki, M., 2012: "Non-negative KDP Estimation
    by Monotone Increasing PHIDP Assumption below Melting Layer". The Seventh
    European Conference on Radar in Meteorology and Hydrology.

    """
    # parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # parse field names
    if psidp_field is None:
        psidp_field = get_field_name('differential_phase')
    if kdp_field is None:
        kdp_field = get_field_name('specific_differential_phase')
    if phidp_field is None:
        phidp_field = get_field_name('differential_phase')

    # parse range resolution, length scale and low-pass filter constraint
    # weight
    dr = _parse_range_resolution(radar, check_uniform=True, verbose=verbose)

    # parse length scale
    if length_scale is None:
        length_scale = dr

    # parse low-pass filter constraint weight
    # this brings the dimensions and magnitude of the low-pass filter
    # constraint in line with the differential phase measurement constraint
    Clpf *= length_scale**4

    # parse near and far range gate propagation differential phase boundary
    # conditions
    bcs = boundary_conditions_maesaka(
        radar, gatefilter=gatefilter, psidp_field=psidp_field, debug=debug,
        verbose=verbose, **kwargs)
    phi_near, phi_far = bcs[:2]
    idx_near, idx_far = bcs[4:]

    # parse total differential phase measurements
    psidp_o = radar.fields[psidp_field]['data']

    if debug:
        N = np.ma.count(psidp_o)
        print('Sample size before filtering: {}'.format(N))

    # mask radar gates indicated by the gate filter
    if gatefilter is not None:
        psidp_o = np.ma.masked_where(gatefilter.gate_excluded, psidp_o)

    # mask any radar gates which are closer (further) than the near (far)
    # boundary condition ranges
    for ray in range(radar.nrays):
        psidp_o[ray, :idx_near[ray]] = np.ma.masked
        psidp_o[ray, idx_far[ray] + 1:] = np.ma.masked

    if debug:
        N = np.ma.count(psidp_o)
        print('Sample size after filtering: {}'.format(N))

    # parse differential phase measurement weight
    Cobs = np.logical_not(np.ma.getmaskarray(psidp_o)).astype(psidp_o.dtype)
    psidp_o = np.ma.filled(psidp_o, fill_value)

    # parse backscattered differential phase data
    if backscatter is None:
        dhv = np.zeros_like(psidp_o, subok=False)

    # parse solver options
    options = {
        'maxiter': kwargs.get('maxiter', 50),
        'gtol': kwargs.get('gtol', 1.0e-5),
        'disp': verbose,
    }
    if debug:
        optimize.show_options(solver='minimize', method=method)

    # parse initial conditions (first guess)
    x0 = np.zeros_like(psidp_o, subok=False).flatten()
    x0.fill(first_guess)

    if verbose:
        print('Cost functional size: {}'.format(x0.size))

    # define arguments for cost functional and its Jacobian (gradient)
    args = (psidp_o, [phi_near, phi_far],
            dhv, dr, Cobs, Clpf,
            finite_order, fill_value,
            proc, debug, verbose)

    if debug:
        start = time.time()

    # minimize the cost functional
    xopt = optimize.minimize(
        _cost_maesaka, x0, args=args, method=method, jac=_jac_maesaka,
        hess=None, hessp=None, bounds=None, constraints=None, callback=None,
        options=options)

    if debug:
        elapsed = time.time() - start
        print('Elapsed time for minimization: {:.0f} sec'.format(elapsed))

    # parse control variables from optimized result
    k = xopt.x.reshape(psidp_o.shape)

    # compute specific differential phase from control variable k in deg/km
    kdp = k**2 / (2.0 * dr) * 1000.0

    if debug:
        print('Min retrieved KDP: {:.2f} deg/km'.format(kdp.min()))
        print('Max retrieved KDP: {:.2f} deg/km'.format(kdp.max()))
        print('Mean retrieved KDP: {:.2f} deg/km'.format(kdp.mean()))

    # create specific differential phase field dictionary and store data
    kdp_dict = get_metadata(kdp_field)
    kdp_dict['data'] = kdp
    kdp_dict['valid_min'] = 0.0
    kdp_dict['Clpf'] = Clpf

    # compute forward and reverse direction propagation differential phase
    phidp_f, phidp_r = _forward_reverse_phidp(
        k, [phi_near, phi_far], verbose=verbose)

    # create forward direction propagation differential phase field dictionary
    # and store data
    phidpf_dict = get_metadata(phidp_field)
    phidpf_dict['data'] = phidp_f
    phidpf_dict['comment'] = 'Retrieved in forward direction'

    # create reverse direction propagation differential phase field dictionary
    # and store data
    phidpr_dict = get_metadata(phidp_field)
    phidpr_dict['data'] = phidp_r
    phidpr_dict['comment'] = 'Retrieved in reverse direction'

    return kdp_dict, phidpf_dict, phidpr_dict


def boundary_conditions_maesaka(
        radar, gatefilter=None, n=20, psidp_field=None, debug=False,
        verbose=False, **kwargs):
    """
    Determine near range gate and far range gate propagation differential phase
    boundary conditions. This follows the method outlined in Maesaka et al.
    (2012), except instead of using the mean we use the median which is less
    susceptible to outliers. This function can also be used to estimate the
    system phase offset.

    Parameters
    ----------
    radar : Radar
        Radar containing total differential phase measurements.
    gatefilter : GateFilter
        A GateFilter indicating radar gates that should be excluded when
        analysing differential phase measurements.
    n : int, optional
        The number of range gates necessary to define the near and far range
        gate boundary conditions. Maesaka et al. (2012) uses a value of 30. If
        this value is too small then a spurious spike in specific differential
        phase close to the radar may be retrieved.
    check_outliers : bool, optional
        True to check for near range gate boundary condition outliers. Outliers
        near the radar are primarily the result of ground clutter returns.
    psidp_field : str, optional
        Field name of total differential phase. If None, the default field name
        must be specified in the Py-ART configuration file.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    Returns
    -------
    phi_near : ndarray
        The near range differential phase boundary condition for each ray.
    phi_far : ndarray
        The far range differential phase boundary condition for each ray.
    range_near : ndarray
        The near range gate in meters for each ray.
    range_far : ndarray
        The far range gate in meters for each ray.
    idx_near : ndarray
        Index of nearest range gate for each ray.
    idx_far : ndarray
        Index of furthest range gate for each ray.

    """
    # parse field names
    if psidp_field is None:
        psidp_field = get_field_name('differential_phase')

    # parse differential phase measurements
    psidp = radar.fields[psidp_field]['data']
    if gatefilter is not None:
        psidp = np.ma.masked_where(gatefilter.gate_excluded, psidp)

    # find contiguous unmasked data along each ray
    slices = np.ma.notmasked_contiguous(psidp, axis=1)
    if debug:
        N = sum([len(slc) for slc in slices if hasattr(slc, '__len__')])
        print('Total number of unique non-masked regions: {}'.format(N))

    phi_near = np.zeros(radar.nrays, dtype=psidp.dtype)
    phi_far = np.zeros(radar.nrays, dtype=psidp.dtype)
    range_near = np.zeros(radar.nrays, dtype=psidp.dtype)
    range_far = np.zeros(radar.nrays, dtype=psidp.dtype)
    idx_near = np.zeros(radar.nrays, dtype=np.int32)
    idx_far = np.zeros(radar.nrays, dtype=np.int32)

    for ray, regions in enumerate(slices):

        # check if all range gates are missing
        if regions is None:
            continue

        # check if all range gates available, i.e., no masked data
        if isinstance(regions, slice):
            regions = [regions]

        # near range gate boundary condition
        for slc in regions:

            # check if enough samples exist in slice
            if slc.stop - slc.start >= n:

                # parse index and range of nearest range gate
                idx = slice(slc.start, slc.start + n)
                idx_near[ray] = idx.start
                range_near[ray] = radar.range['data'][idx][0]

                # parse data for linear regression and compute slope
                x = radar.range['data'][idx]
                y = psidp[ray, idx]
                slope = stats.linregress(x, y)[0]

                # if linear regression slope is positive, set near range gate
                # boundary condition to first differential phase measurement,
                # otherwise use the median value
                if slope > 0.0:
                    phi_near[ray] = y[0]
                else:
                    phi_near[ray] = np.median(y)
                break

            else:
                continue

        # far range gate boundary condition
        for slc in reversed(regions):

            if slc.stop - slc.start >= n:

                # parse index and range of furthest range gate
                idx = slice(slc.stop - n, slc.stop)
                idx_far[ray] = idx.stop
                range_far[ray] = radar.range['data'][idx][-1]

                # parse data for linear regression and compute slope
                x = radar.range['data'][idx]
                y = psidp[ray, idx]
                slope = stats.linregress(x, y)[0]

                # if linear regression slope is positive, set far range gate
                # boundary condition to last differential phase measurement,
                # otherwise use the median value
                if slope > 0.0:
                    phi_far[ray] = y[-1]
                else:
                    phi_far[ray] = np.median(y)
                break

            else:
                continue

    # check for outliers in the near range boundary conditions, e.g., ground
    # clutter can introduce spurious values in certain rays

    # do not include missing values in the analysis
    phi_near_valid = phi_near[phi_near != 0.0]

    # skip the check if there are no valid values
    if kwargs.get('check_outliers', True) and (phi_near_valid.size != 0):

        # bin and count near range boundary condition values, i.e., create
        # a distribution of values
        # the default bin width is 5 deg
        counts, edges = np.histogram(
            phi_near_valid, bins=144, range=(-360, 360), density=False)

        # assume that the maximum counts corresponds to the maximum in the
        # system phase distribution
        # this assumption breaks down if there are a significant amount of
        # near range boundary conditions characterized by ground clutter
        system_phase_peak_left = edges[counts.argmax()]
        system_phase_peak_right = edges[counts.argmax() + 1]

        if debug:
            print('Peak of system phase distribution: {:.0f} deg'.format(
                system_phase_peak_left))

        # determine left edge location of system phase distribution
        # we consider five counts or less to be insignificant
        is_left_side = np.logical_and(
            edges[:-1] < system_phase_peak_left, counts <= 5)
        left_edge = edges[:-1][is_left_side][-1]

        # determine right edge location of system phase distribution
        # we consider five counts or less to be insignificant
        is_right_side = np.logical_and(
            edges[1:] > system_phase_peak_right, counts <= 5)
        right_edge = edges[1:][is_right_side][0]

        if debug:
            print('Left edge of system phase distribution: {:.0f} deg'.format(
                left_edge))
            print('Right edge of system phase distribution: {:.0f} deg'.format(
                right_edge))

        # define the system phase offset as the median value of the system
        # phase distriubion
        is_system_phase = np.logical_or(
            phi_near_valid >= left_edge, phi_near_valid <= right_edge)
        system_phase_offset = np.median(phi_near_valid[is_system_phase])

        if debug:
            print('Estimated system phase offset: {:.0f} deg'.format(
                system_phase_offset))

        for ray, bc in enumerate(phi_near):

            # if near range boundary condition does not draw from system phase
            # distribution then set it to the system phase offset
            if bc < left_edge or bc > right_edge:
                phi_near[ray] = system_phase_offset

    # check for unphysical boundary conditions, i.e., propagation differential
    # phase should monotonically increase from the near boundary to the far
    # boundary
    is_unphysical = phi_far - phi_near < 0.0
    phi_far[is_unphysical] = phi_near[is_unphysical]
    if verbose:
        N = is_unphysical.sum()
        print('Rays with unphysical boundary conditions: {}'.format(N))

    return phi_near, phi_far, range_near, range_far, idx_near, idx_far


def _cost_maesaka(x, psidp_o, bcs, dhv, dr, Cobs, Clpf, finite_order,
                  fill_value, proc, debug=False, verbose=False):
    """
    Compute the value of the cost functional similar to equations (12)-(15) in
    Maesaka et al. (2012).

    Parameters
    ----------
    x : ndarray
        Analysis vector containing control variable k.
    psidp_o : ndarray
        Total differential phase measurements.
    bcs : array_like
        The near and far range gate propagation differential phase boundary
        conditions.
    dhv : ndarray
        Backscatter differential phase.
    dr : float
        Range resolution in meters.
    Cobs : ndarray
        The differential phase measurement constraint weights. The weight
        should vanish where no differential phase measurements are available.
    Clpf : float
        The low-pass filter (radial smoothness) constraint weight as in
        equation (15) of Maesaka et al. (2012).
    finite_order : 'low' or 'high'
        The finite difference accuracy to use when computing derivatives.
    fill_value : float
        Value indicating missing or bad data in radar field data.
    proc : int
        The number of parallel threads (CPUs) to use.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print progress information, False to suppress.

    Returns
    -------
    J : float
        Value of total cost functional.

    """
    # parse control variable k from analysis vector
    nr, ng = psidp_o.shape
    k = x.reshape(nr, ng)

    # parse near and far range gate boundary conditions
    phi_near, phi_far = bcs

    # compute forward direction propagation differential phase from control
    # variable k
    phi_fa = np.zeros_like(k, subok=False)
    phi_fa[:, 1:] = np.cumsum(k[:, :-1]**2, axis=1)

    # compute reverse direction propagation differential phase from control
    # variable k
    phi_ra = np.zeros_like(k, subok=False)
    phi_ra[:, :-1] = np.cumsum(k[:, :0:-1]**2, axis=1)[:, ::-1]

    # compute forward and reverse propagation differential phase
    # from total differential phase observations
    phi_fo = psidp_o - dhv - phi_near[:, np.newaxis].repeat(ng, axis=1)
    phi_ro = phi_far[:, np.newaxis].repeat(ng, axis=1) - psidp_o + dhv

    # cost: forward direction differential phase observations
    Jof = 0.5 * np.sum(Cobs * (phi_fa - phi_fo)**2)

    # cost: reverse direction differential phase observations
    Jor = 0.5 * np.sum(Cobs * (phi_ra - phi_ro)**2)

    # prepare control variable k for Cython function
    k = np.ascontiguousarray(k, dtype=np.float64)

    # compute low-pass filter term, i.e., second order derivative of k with
    # respect to range
    d2kdr2 = np.empty_like(k)
    _kdp_proc.lowpass_maesaka_term(k, dr, finite_order, d2kdr2)

    # cost: low-pass filter, i.e., radial smoothness
    Jlpf = 0.5 * np.sum(Clpf * (d2kdr2)**2)

    # compute value of total cost functional
    J = Jof + Jor + Jlpf

    if verbose:
        print('Forward direction observation cost : {:1.3e}'.format(Jof))
        print('Reverse direction observation cost : {:1.3e}'.format(Jor))
        print('Low-pass filter cost ............. : {:1.3e}'.format(Jlpf))
        print('Total cost ....................... : {:1.3e}'.format(J))

    return J


def _jac_maesaka(x, psidp_o, bcs, dhv, dr, Cobs, Clpf, finite_order,
                 fill_value, proc, debug=False, verbose=False):
    """
    Compute the Jacobian (gradient) of the cost functional similar to equations
    (16)-(18) in Maesaka et al. (2012).

    Parameters
    ----------
    x : ndarray
        Analysis vector containing control variable k.
    psidp_o : ndarray
        Total differential phase measurements.
    bcs : array_like
        The near and far range gate propagation differential phase boundary
        conditions.
    dhv : ndarray
        Backscatter differential phase.
    dr : float
        Range resolution in meters.
    Cobs : ndarray
        The differential phase measurement constraint weights. The weight
        should vanish where no differential phase measurements are available.
    Clpf : float
        The low-pass filter (radial smoothness) constraint weight as in
        equation (15) of Maesaka et al. (2012).
    finite_order : 'low' or 'high'
        The finite difference accuracy to use when computing derivatives.
    fill_value : float
        Value indicating missing or bad data in radar field data.
    proc : int
        The number of parallel threads (CPUs) to use.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print progress information, False to suppress.

    Returns
    -------
    jac : ndarray
        Jacobian of the cost functional.

    """
    # parse control variable k from analysis vector
    nr, ng = psidp_o.shape
    k = x.reshape(nr, ng)

    # parse near and far range gate boundary conditions
    phi_near, phi_far = bcs

    # compute forward direction propagation differential phase from control
    # variable k
    phi_fa = np.zeros_like(k, subok=False)
    phi_fa[:, 1:] = np.cumsum(k[:, :-1]**2, axis=1)

    # compute reverse direction propagation differential phase from control
    # variable k
    phi_ra = np.zeros_like(k, subok=False)
    phi_ra[:, :-1] = np.cumsum(k[:, :0:-1]**2, axis=1)[:, ::-1]

    # compute forward and reverse propagation differential phase
    # from total differential phase observations
    phi_fo = psidp_o - dhv - phi_near[:, np.newaxis].repeat(ng, axis=1)
    phi_ro = phi_far[:, np.newaxis].repeat(ng, axis=1) - psidp_o + dhv

    # prepare control variable k for Cython functions
    k = np.ascontiguousarray(k, dtype=np.float64)

    # cost: forward direction differential phase observations
    dJofdk = np.zeros_like(k, subok=False)
    dJofdk[:, :-1] = 2.0 * k[:, :-1] * np.cumsum(
        (Cobs[:, 1:] * (phi_fa[:, 1:] - phi_fo[:, 1:]))[:, ::-1],
        axis=1)[:, ::-1]

    # cost: reverse direction differential phase observations
    dJordk = np.zeros_like(k, subok=False)
    dJordk[:, 1:] = 2.0 * k[:, 1:] * np.cumsum(
        (Cobs[:, :-1] * (phi_ra[:, :-1] - phi_ro[:, :-1])), axis=1)

    # compute low-pass filter term, i.e., second order derivative of k with
    # respect to range
    d2kdr2 = np.empty_like(k)
    _kdp_proc.lowpass_maesaka_term(k, dr, finite_order, d2kdr2)

    # compute gradients of Jlpf with respect to the control variable k
    dJlpfdk = np.empty_like(d2kdr2)
    _kdp_proc.lowpass_maesaka_jac(d2kdr2, dr, Clpf, finite_order, dJlpfdk)

    # sum control variable derivative components
    dJdk = dJofdk + dJordk + dJlpfdk
    jac = dJdk.flatten()

    # compute the vector norm of the Jacobian
    if verbose:
        mag = np.linalg.norm(jac, ord=None, axis=None)
        print('Vector norm of Jacobian: {:1.3e}'.format(mag))

    return jac


def _forward_reverse_phidp(k, bcs, verbose=False):
    """
    Compute the forward and reverse direction propagation differential phases
    from the control variable k and boundary conditions following equations (1)
    and (7) in Maesaka et al. (2012).

    Parameters
    ----------
    k : ndarray
        Control variable k of the Maesaka et al. (2012) method. The control
        variable k is proportional to the square root of specific differential
        phase.
    bcs : array_like
        The near and far range gate boundary conditions.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    Returns
    -------
    phidp_f : ndarray
        Forward direction propagation differential phase.
    phidp_r : ndarray
        Reverse direction propagation differential phase.

    """
    # parse near and far range gate boundary conditions
    nr, ng = k.shape
    phi_near, phi_far = bcs

    # compute forward direction propagation differential phase
    phi_f = np.zeros_like(k, subok=False)
    phi_f[:, 1:] = np.cumsum(k[:, :-1]**2, axis=1)
    phidp_f = phi_f + phi_near[:, np.newaxis].repeat(ng, axis=1)

    # compute reverse direction propagation differential phase
    phi_r = np.zeros_like(k, subok=False)
    phi_r[:, :-1] = np.cumsum(k[:, :0:-1]**2, axis=1)[:, ::-1]
    phidp_r = phi_far[:, np.newaxis].repeat(ng, axis=1) - phi_r

    # check quality of retrieval by comparing forward and reverse directions
    if verbose:
        phidp_mbe = np.ma.mean(phidp_f - phidp_r)
        phidp_mae = np.ma.mean(np.abs(phidp_f - phidp_r))
        print('Forward-reverse PHIDP MBE: {:.2f} deg'.format(phidp_mbe))
        print('Forward-reverse PHIDP MAE: {:.2f} deg'.format(phidp_mae))

    return phidp_f, phidp_r


def _parse_range_resolution(
        radar, check_uniform=True, atol=1.0, verbose=False):
    """
    Parse the radar range gate resolution.

    Parameters
    ----------
    radar : Radar
        Radar containing range data.
    check_uniform : bool, optional
        True to check if all range gates are equally spaced, and if so return
        a scalar value for range resolution. If False, the resolution between
        each range gate is returned.
    atol : float, optional
        The absolute tolerance in meters allowed for discrepancies in range
        gate spacings. Only applicable when check_uniform is True. This
        parameter may be necessary to catch instances where range gate spacings
        differ by a few meters or so.
    verbose : bool, optional
        True to print the range gate resolution. Only valid if check_uniform is
        True.

    Returns
    -------
    dr : float or ndarray
        The radar range gate spacing in meters.

    """
    # parse radar range gate spacings
    dr = np.diff(radar.range['data'], n=1)

    # check for uniform resolution
    if check_uniform and np.allclose(np.diff(dr), 0.0, atol=atol):
        dr = dr[0]
        if verbose:
            print('Range resolution: {:.2f} m'.format(dr))
    else:
        raise ValueError('Radar gate spacing is not uniform')

    return dr
