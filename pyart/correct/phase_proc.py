"""
Utilities for working with phase data.

Code based upon algorithm descriped in:
S. E. Giangrande et al, J. of Atmos. and Ocean. Tech., 2013, 30, 1716.

Adapted by Scott Collis and Scott Giangrande, refactored by Jonathan Helmus.

"""

import copy
from time import time

import numpy as np
from numpy import ma
import scipy.ndimage

from ..config import get_fillvalue, get_field_name, get_metadata
from ..filters import GateFilter
from ..util import rolling_window


def det_sys_phase(radar, ncp_lev=0.4, rhohv_lev=0.6,
                  ncp_field=None, rhv_field=None, phidp_field=None):
    """
    Determine the system phase.

    Parameters
    ----------
    radar : Radar
        Radar object for which to determine the system phase.
    ncp_lev : float, optional
        Miminum normal coherent power level. Regions below this value will
        not be included in the phase calculation.
    rhohv_lev : float, optional
        Miminum copolar coefficient level. Regions below this value will not
        be included in the phase calculation.
    ncp_field, rhv_field, phidp_field : str, optional
        Field names within the radar object which represent the normal
        coherent power, the copolar coefficient, and the differential phase
        shift. A value of None for any of these parameters will use the
        default field name as defined in the Py-ART configuration file.

    Returns
    -------
    sys_phase : float or None
        Estimate of the system phase. None is not estimate can be made.

    """
    # parse the field parameters
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')
    if rhv_field is None:
        rhv_field = get_field_name('cross_correlation_ratio')
    if phidp_field is None:
        phidp_field = get_field_name('differential_phase')

    ncp = radar.fields[ncp_field]['data'][:, 30:]
    rhv = radar.fields[rhv_field]['data'][:, 30:]
    phidp = radar.fields[phidp_field]['data'][:, 30:]
    last_ray_idx = radar.sweep_end_ray_index['data'][0]
    return _det_sys_phase(ncp, rhv, phidp, last_ray_idx, ncp_lev,
                          rhohv_lev)


def _det_sys_phase(ncp, rhv, phidp, last_ray_idx, ncp_lev=0.4,
                   rhv_lev=0.6):
    """ Determine the system phase, see :py:func:`det_sys_phase`. """
    good = False
    phases = []
    for radial in range(last_ray_idx + 1):
        meteo = np.logical_and(ncp[radial, :] > ncp_lev,
                               rhv[radial, :] > rhv_lev)
        mpts = np.where(meteo)
        if len(mpts[0]) > 25:
            good = True
            msmth_phidp = smooth_and_trim(phidp[radial, mpts[0]], 9)
            phases.append(msmth_phidp[0:25].min())
    if not good:
        return None
    return np.median(phases)


def fzl_index(fzl, ranges, elevation, radar_height):
    """
    Return the index of the last gate below a given altitude.

    Parameters
    ----------
    fzl : float
        Maximum altitude.
    ranges : array
        Range to measurement volume/gate in meters.
    elevation : float
        Elevation of antenna in degrees.
    radar_height :
        Altitude of radar in meters.

    Returns
    -------
    idx : int
        Index of last gate which has an altitude below `fzl`.

    Notes
    -----
    Standard atmosphere is assumed, R = 4 / 3 * Re

    """
    Re = 6371.0 * 1000.0
    p_r = 4.0 * Re / 3.0
    z = radar_height + (ranges ** 2 + p_r ** 2 + 2.0 * ranges * p_r *
                        np.sin(elevation * np.pi / 180.0)) ** 0.5 - p_r
    # Make sure the freezing level isn't under the radar!
    # Return the minimum window size for the 5-pt filter
    if np.all(z > fzl):
        return 6
    else:
        return np.where(z < fzl)[0].max()


def det_process_range(radar, sweep, fzl, doc=10):
    """
    Determine the processing range for a given sweep.

    Queues the radar and returns the indices which can be used to slice
    the radar fields and select the desired sweep with gates which are
    below a given altitude.

    Parameters
    ----------
    radar : Radar
        Radar object from which ranges will be determined.
    sweep : int
        Sweep (0 indexed) for which to determine processing ranges.
    fzl : float
        Maximum altitude in meters. The determined range will not include
        gates which are above this limit.
    doc : int, optional
        Minimum number of gates which will be excluded from the determined
        range.

    Returns
    -------
    gate_end : int
        Index of last gate below `fzl` and satisfying the `doc` parameter.
    ray_start : int
        Ray index which defines the start of the region.
    ray_end : int
        Ray index which defined the end of the region.

    """
    # determine the index of the last valid gate
    ranges = radar.range['data']
    elevation = radar.fixed_angle['data'][sweep]
    radar_height = radar.altitude['data']
    gate_end = fzl_index(fzl, ranges, elevation, radar_height)
    if doc is not None:
        gate_end = min(gate_end, len(ranges) - doc)
    else:
        gate_end = min(gate_end, len(ranges))

    ray_start = radar.sweep_start_ray_index['data'][sweep]
    ray_end = radar.sweep_end_ray_index['data'][sweep] + 1
    return gate_end, ray_start, ray_end


def snr(line, wl=11):
    """ Return the signal to noise ratio after smoothing. """
    signal = smooth_and_trim(line, window_len=wl)
    _noise = smooth_and_trim(np.sqrt((line - signal) ** 2), window_len=wl)
    return abs(signal) / _noise


def smooth_masked(raw_data, wind_len=11, min_valid=6, wind_type='median'):
    """
    Smoothes the data using a rolling window.
    data with less than n valid points is masked.

    Parameters
    ----------
    raw_data : float masked array
        The data to smooth.
    win_len : float
        Length of the moving window.
    min_valid : float
        Minimum number of valid points for the smoothing to be valid.
    wind_type : str
        Type of window. Can be median or mean.

    Returns
    -------
    data_smooth : float masked array
        Smoothed data.

    """
    valid_wind = ['median', 'mean']
    if wind_type not in valid_wind:
        raise ValueError(
            "Window " + win_type + " is none of " + ' '.join(valid_wind))

    # we want an odd window
    if wind_len % 2 == 0:
        wind_len += 1
    half_wind = int((wind_len-1)/2)

    # initialize smoothed data
    nrays, nbins = np.shape(raw_data)
    data_smooth = np.ma.zeros((nrays, nbins))
    data_smooth[:] = np.ma.masked
    data_smooth.set_fill_value(get_fillvalue())

    mask = np.ma.getmaskarray(raw_data)
    valid = np.logical_not(mask)

    mask_wind = rolling_window(mask, wind_len)
    valid_wind = np.logical_not(mask_wind).astype(int)
    nvalid = np.sum(valid_wind, -1)

    data_wind = rolling_window(raw_data, wind_len)

    # check which gates are valid
    ind_valid = np.logical_and(
        nvalid >= min_valid, valid[:, half_wind:-half_wind]).nonzero()

    data_smooth[ind_valid[0], ind_valid[1]+half_wind] = (
        eval('np.ma.' + wind_type + '(data_wind, axis=-1)')[ind_valid])

    return data_smooth


def unwrap_masked(lon, centered=False, copy=True):
    """
    Unwrap a sequence of longitudes or headings in degrees.

    Parameters
    ----------
    lon : array
        Longtiudes or heading in degress. If masked output will also be
        masked.
    centered : bool, optional
        Center the unwrapping as close to zero as possible.
    copy : bool, optional.
        True to return a copy, False will avoid a copy when possible.

    Returns
    -------
    unwrap : array
        Array of unwrapped longtitudes or headings, in degrees.

    """
    masked_input = ma.isMaskedArray(lon)
    if masked_input:
        fill_value = lon.fill_value
        # masked_invalid loses the original fill_value (ma bug, 2011/01/20)
    lon = np.ma.masked_invalid(lon).astype(float)
    if lon.ndim != 1:
        raise ValueError("Only 1-D sequences are supported")
    if lon.shape[0] < 2:
        return lon
    x = lon.compressed()
    if len(x) < 2:
        return lon
    w = np.zeros(x.shape[0] - 1, int)
    ld = np.diff(x)
    np.putmask(w, ld > 180, -1)
    np.putmask(w, ld < -180, 1)
    x[1:] += (w.cumsum() * 360.0)
    if centered:
        x -= 360 * np.round(x.mean() / 360.0)
    if lon.mask is ma.nomask:
        lon[:] = x
    else:
        lon[~lon.mask] = x
    if masked_input:
        lon.fill_value = fill_value
        return lon
    else:
        return lon.filled(np.nan)


# this function adapted from the Scipy Cookbook:
# http://www.scipy.org/Cookbook/SignalSmooth
def smooth_and_trim(x, window_len=11, window='hanning'):
    """
    Smooth data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    Parameters
    ----------
    x : array
        The input signal.
    window_len : int, optional
        The dimension of the smoothing window; should be an odd integer.
    window : str
        The type of window from 'flat', 'hanning', 'hamming', 'bartlett',
        'blackman' or 'sg_smooth'. A flat window will produce a moving
        average smoothing.

    Returns
    -------
    y : array
        The smoothed signal with length equal to the input signal.

    """
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len < 3:
        return x
    valid_windows = ['flat', 'hanning', 'hamming', 'bartlett', 'blackman',
                     'sg_smooth']
    if window not in valid_windows:
        raise ValueError("Window is on of " + ' '.join(valid_windows))

    s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]

    if window == 'flat':  # moving average
        w = np.ones(int(window_len), 'd')
    elif window == 'sg_smooth':
        w = np.array([0.1, .25, .3, .25, .1])
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')

    return y[int(window_len / 2):len(x) + int(window_len / 2)]


# adapted smooth and trim function to work with 2dimensional arrays
def smooth_and_trim_scan(x, window_len=11, window='hanning'):
    """
    Smooth data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    Parameters
    ----------
    x : ndarray
        The input signal.
    window_len : int, optional
        The dimension of the smoothing window; should be an odd integer.
    window : str, optional
        The type of window from 'flat', 'hanning', 'hamming', 'bartlett',
        'blackman' or 'sg_smooth'. A flat window will produce a moving
        average smoothing.

    Returns
    -------
    y : ndarray
        The smoothed signal with length equal to the input signal.

    """
    from scipy.ndimage.filters import convolve1d

    if x.ndim != 2:
        raise ValueError("smooth only accepts 2 dimension arrays.")
    if x.shape[1] < window_len:
        mess = "Input dimension 1 needs to be bigger than window size."
        raise ValueError(mess)
    if window_len < 3:
        return x
    valid_windows = ['flat', 'hanning', 'hamming', 'bartlett', 'blackman',
                     'sg_smooth']
    if window not in valid_windows:
        raise ValueError("Window is on of " + ' '.join(valid_windows))

    if window == 'flat':  # moving average
        w = np.ones(int(window_len), 'd')
    elif window == 'sg_smooth':
        w = np.array([0.1, .25, .3, .25, .1])
    else:
        w = eval('np.' + window + '(window_len)')

    y = convolve1d(x, w / w.sum(), axis=1)

    return y


def noise(line, wl=11):
    """ Return the noise after smoothing. """
    signal = smooth_and_trim(line, window_len=wl)
    _noise = np.sqrt((line - signal) ** 2)
    return _noise


def get_phidp_unf(radar, ncp_lev=0.4, rhohv_lev=0.6, debug=False, ncpts=20,
                  doc=-10, overide_sys_phase=False, sys_phase=-135,
                  nowrap=None, refl_field=None, ncp_field=None,
                  rhv_field=None, phidp_field=None):
    """
    Get Unfolded Phi differential phase

    Parameters
    ----------
    radar : Radar
        The input radar.
    ncp_lev : float, optional
        Miminum normal coherent power level. Regions below this value will
        not be included in the calculation.
    rhohv_lev : float, optional
        Miminum copolar coefficient level. Regions below this value will not
        be included in the calculation.
    debug : bool, optional
        True to print debugging information, False to supress printing.
    ncpts : int, optional
        Minimum number of points in a ray. Regions within a ray smaller than
        this or beginning before this gate number are excluded from
        calculations.
    doc : int or None, optional
        Index of first gate not to include in field data, None include all.
    overide_sys_phase : bool, optional
        True to use `sys_phase` as the system phase. False will determine a
        value automatically.
    sys_phase : float, optional
        System phase, not used if overide_sys_phase is False.
    nowrap : int or None, optional
        Gate number where unwrapping should begin. `None` will unwrap all
        gates.
    refl_field ncp_field, rhv_field, phidp_field : str, optional
        Field names within the radar object which represent the horizonal
        reflectivity, normal coherent power, the copolar coefficient, and the
        differential phase shift. A value of None for any of these parameters
        will use the default field name as defined in the Py-ART
        configuration file.

    Returns
    -------
    cordata : array
        Unwrapped phi differential phase.

    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')
    if rhv_field is None:
        rhv_field = get_field_name('cross_correlation_ratio')
    if phidp_field is None:
        phidp_field = get_field_name('differential_phase')

    if doc is not None:
        my_phidp = radar.fields[phidp_field]['data'][:, 0:doc]
        my_rhv = radar.fields[rhv_field]['data'][:, 0:doc]
        my_ncp = radar.fields[ncp_field]['data'][:, 0:doc]
        my_z = radar.fields[refl_field]['data'][:, 0:doc]
    else:
        my_phidp = radar.fields[phidp_field]['data']
        my_rhv = radar.fields[rhv_field]['data']
        my_ncp = radar.fields[ncp_field]['data']
        my_z = radar.fields[refl_field]['data']
    t = time()
    if overide_sys_phase:
        system_zero = sys_phase
    else:
        system_zero = det_sys_phase(
            radar, ncp_field=ncp_field, rhv_field=rhv_field,
            phidp_field=phidp_field)
        if system_zero is None:
            system_zero = sys_phase
    cordata = np.zeros(my_rhv.shape, dtype=float)
    for radial in range(my_rhv.shape[0]):
        my_snr = snr(my_z[radial, :])
        notmeteo = np.logical_or(np.logical_or(
            my_ncp[radial, :] < ncp_lev,
            my_rhv[radial, :] < rhohv_lev), my_snr < 10.0)
        x_ma = ma.masked_where(notmeteo, my_phidp[radial, :])
        try:
            ma.notmasked_contiguous(x_ma)
            for slc in ma.notmasked_contiguous(x_ma):
                # so trying to get rid of clutter and small things that
                # should not add to phidp anyway
                if slc.stop - slc.start < ncpts or slc.start < ncpts:
                    x_ma.mask[slc.start - 1:slc.stop + 1] = True
            c = 0
        except TypeError:  # non sequence, no valid regions
            c = 1  # ie do nothing
            x_ma.mask = True
        except AttributeError:
            # sys.stderr.write('No Valid Regions, ATTERR \n ')
            # sys.stderr.write(myfile.times['time_end'].strftime('%Y-%m-%dT%H:%M:%SZ') + '\n')
            # print x_ma
            # print x_ma.mask
            c = 1  # also do nothing
            x_ma.mask = True
        if nowrap is not None:
            # Start the unfolding a bit later in order to avoid false
            # jumps based on clutter
            unwrapped = copy.deepcopy(x_ma)
            end_unwrap = unwrap_masked(x_ma[nowrap::], centered=False)
            unwrapped[nowrap::] = end_unwrap
        else:
            unwrapped = unwrap_masked(x_ma, centered=False)
        # end so no clutter expected
        system_max = unwrapped[np.where(np.logical_not(
            notmeteo))][-10:-1].mean() - system_zero
        unwrapped_fixed = np.zeros(len(x_ma), dtype=float)
        based = unwrapped-system_zero
        based[0] = 0.0
        notmeteo[0] = False
        based[-1] = system_max
        notmeteo[-1] = False
        unwrapped_fixed[np.where(np.logical_not(based.mask))[0]] = \
            based[np.where(np.logical_not(based.mask))[0]]
        if len(based[np.where(np.logical_not(based.mask))[0]]) > 11:
            unwrapped_fixed[np.where(based.mask)[0]] = \
                np.interp(np.where(based.mask)[0],
                          np.where(np.logical_not(based.mask))[0],
                          smooth_and_trim(based[np.where(
                              np.logical_not(based.mask))[0]]))
        else:
            unwrapped_fixed[np.where(based.mask)[0]] = \
                np.interp(np.where(based.mask)[0],
                          np.where(np.logical_not(based.mask))[0],
                          based[np.where(np.logical_not(based.mask))[0]])
        if c != 1:
            cordata[radial, :] = unwrapped_fixed
        else:
            cordata[radial, :] = np.zeros(my_rhv.shape[1])
    if debug:
        print("Exec time: ", time() - t)
    return cordata


def construct_A_matrix(n_gates, filt):
    """
    Construct a row-augmented A matrix. Equation 5 in Giangrande et al, 2012.

    A is a block matrix given by:

    .. math::

        \\bf{A} = \\begin{bmatrix} \\bf{I} & \\bf{-I} \\\\\\\\
                  \\bf{-I} & \\bf{I} \\\\\\\\ \\bf{Z}
                  & \\bf{M} \\end{bmatrix}

    where
        :math:`\\bf{I}` is the identity matrix
        :math:`\\bf{Z}` is a matrix of zeros
        :math:`\\bf{M}` contains our differential constraints.

    Each block is of shape n_gates by n_gates making
    shape(:math:`\\bf{A}`) = (3 * n, 2 * n).

    Note that :math:`\\bf{M}` contains some side padding to deal with edge
    issues.

    Parameters
    ----------
    n_gates : int
        Number of gates, determines size of identity matrix.
    filt : array
        Input filter.

    Returns
    -------
    a : matrix
        Row-augmented A matrix.

    """
    Identity = np.eye(n_gates)
    filter_length = len(filt)
    M_matrix_middle = np.diag(np.ones(n_gates - filter_length + 1), k=0) * 0.0
    posn = np.linspace(-1.0 * (filter_length - 1) / 2, (filter_length - 1)/2,
                       filter_length)
    for diag in range(filter_length):
        M_matrix_middle = M_matrix_middle + np.diag(
            np.ones(int(n_gates - filter_length + 1 - np.abs(posn[diag]))),
            k=int(posn[diag])) * filt[diag]
    side_pad = (filter_length - 1) // 2
    M_matrix = np.bmat(
        [np.zeros([n_gates-filter_length + 1, side_pad], dtype=float),
         M_matrix_middle, np.zeros(
             [n_gates-filter_length+1, side_pad], dtype=float)])
    Z_matrix = np.zeros([n_gates - filter_length + 1, n_gates])
    return np.bmat([[Identity, -1.0 * Identity], [Identity, Identity],
                    [Z_matrix, M_matrix]])


def construct_B_vectors(phidp_mod, z_mod, filt, coef=0.914, dweight=60000.0):
    """
    Construct B vectors. See Giangrande et al, 2012.

    Parameters
    ----------
    phidp_mod : 2D array
        Phi differential phases.
    z_mod : 2D array.
       Reflectivity, modified as needed.
    filt : array
        Input filter.
    coef : float, optional.
        Cost coefficients.
    dweight : float, optional.
        Weights.

    Returns
    -------
    b : matrix
        Matrix containing B vectors.

    """
    n_gates = phidp_mod.shape[1]
    n_rays = phidp_mod.shape[0]
    filter_length = len(filt)
    side_pad = (filter_length - 1) // 2
    top_of_B_vectors = np.bmat([[-phidp_mod, phidp_mod]])
    data_edges = np.bmat([phidp_mod[:, 0:side_pad],
                          np.zeros([n_rays, n_gates-filter_length+1]),
                          phidp_mod[:, -side_pad:]])
    ii = filter_length - 1
    jj = data_edges.shape[1] - 1
    list_corrl = np.zeros([n_rays, jj - ii + 1])
    for count in range(list_corrl.shape[1]):
        list_corrl[:, count] = -1.0 * (
            np.array(filt) * (np.asarray(
                data_edges))[:, count:count+ii+1]).sum(axis=1)

    sct = (((10.0 ** (0.1 * z_mod)) ** coef / dweight))[:, side_pad: -side_pad]
    sct[np.where(sct < 0.0)] = 0.0
    sct[:, 0:side_pad] = list_corrl[:, 0:side_pad]
    sct[:, -side_pad:] = list_corrl[:, -side_pad:]
    B_vectors = np.bmat([[top_of_B_vectors, sct]])
    return B_vectors


def LP_solver_cvxopt(A_Matrix, B_vectors, weights, solver='glpk'):
    """
    Solve the Linear Programming problem given in Giangrande et al, 2012 using
    the CVXOPT module.

    Parameters
    ----------
    A_Matrix : matrix
        Row augmented A matrix, see :py:func:`construct_A_matrix`
    B_vectors : matrix
        Matrix containing B vectors, see :py:func:`construct_B_vectors`
    weights : array
        Weights.
    solver : str or None
        LP solver backend to use, choices are 'glpk', 'mosek' or None to use
        the conelp function in CVXOPT.  'glpk' and 'mosek' are only available
        if they are installed and CVXOPT was build with the correct bindings.

    Returns
    -------
    soln : array
        Solution to LP problem.

    See Also
    --------
    LP_solver_pyglpk : Solve LP problem using the PyGLPK module.
    LP_solver_cylp : Solve LP problem using the cylp module.
    LP_solver_cylp_mp : Solve LP problem using the cylp module
                        using multi processes.

    """
    from cvxopt import matrix, solvers
    n_gates = weights.shape[1] // 2
    n_rays = B_vectors.shape[0]
    mysoln = np.zeros([n_rays, n_gates])

    G = matrix(np.bmat([[-A_Matrix], [-np.eye(2 * n_gates)]]))
    h_array = np.zeros(5 * n_gates - 4)
    for raynum in range(n_rays):
        c = matrix(weights[raynum]).T
        h_array[:3 * n_gates - 4] = -B_vectors[raynum]
        h = matrix(h_array)
        sol = solvers.lp(c, G, h, solver=solver)
        # XXX when a solution is not found sol is None, need to check and
        # deal with this...

        # extract the solution
        this_soln = np.zeros(n_gates)
        for i in range(n_gates):
            this_soln[i] = sol['x'][i + n_gates]

        # apply smoothing filter and record in output array
        mysoln[raynum, :] = smooth_and_trim(this_soln, window_len=5,
                                            window='sg_smooth')
    return mysoln


def LP_solver_pyglpk(A_Matrix, B_vectors, weights, it_lim=7000, presolve=True,
                     really_verbose=False):
    """
    Solve the Linear Programming problem given in Giangrande et al, 2012 using
    the PyGLPK module.

    Parameters
    ----------
    A_Matrix : matrix
        Row augmented A matrix, see :py:func:`construct_A_matrix`
    B_vectors : matrix
        Matrix containing B vectors, see :py:func:`construct_B_vectors`
    weights : array
        Weights.
    it_lim : int, optional
        Simplex iteration limit.
    presolve : bool, optional
        True to use the LP presolver.
    really_verbose : bool, optional
        True to print LPX messaging. False to suppress.

    Returns
    -------
    soln : array
        Solution to LP problem.

    See Also
    --------
    LP_solver_cvxopt : Solve LP problem using the CVXOPT module.
    LP_solver_cylp : Solve LP problem using the cylp module.
    LP_solver_cylp_mp : Solve LP problem using the cylp module
                        using multi processes.

    """
    import glpk

    if really_verbose:
        message_state = glpk.LPX.MSG_ON
    else:
        message_state = glpk.LPX.MSG_OFF
    n_gates = weights.shape[1] // 2
    n_rays = B_vectors.shape[0]
    mysoln = np.zeros([n_rays, n_gates])
    lp = glpk.LPX()  # Create empty problem instance
    lp.name = 'LP_MIN'  # Assign symbolic name to problem
    lp.obj.maximize = False  # Set this as a maximization problem
    lp.rows.add(2 * n_gates + n_gates - 4)  # Append rows
    lp.cols.add(2 * n_gates)
    glpk.env.term_on = True
    for cur_row in range(2 * n_gates + n_gates - 4):
        lp.rows[cur_row].matrix = list(np.squeeze(np.asarray(
            A_Matrix[cur_row, :])))
    for i in range(2 * n_gates):
        lp.cols[i].bounds = 0.0, None
    for raynum in range(n_rays):
        this_soln = np.zeros(n_gates)
        for i in range(2 * n_gates + n_gates - 4):
            lp.rows[i].bounds = B_vectors[raynum, i], None
        for i in range(2 * n_gates):
            lp.obj[i] = weights[raynum, i]
        lp.simplex(msg_lev=message_state, meth=glpk.LPX.PRIMAL,
                   it_lim=it_lim, presolve=presolve)
        for i in range(n_gates):
            this_soln[i] = lp.cols[i+n_gates].primal
        mysoln[raynum, :] = smooth_and_trim(this_soln, window_len=5,
                                            window='sg_smooth')
    return mysoln


def solve_cylp(model, B_vectors, weights, ray, chunksize):
    """
    Worker process for LP_solver_cylp_mp.

    Parameters
    ----------
    model : CyClpModel
        Model of the LP Problem, see :py:func:`LP_solver_cylp_mp`
    B_vectors : matrix
        Matrix containing B vectors, see :py:func:`construct_B_vectors`
    weights : array
        Weights.
    ray : int
        Starting ray.
    chunksize : int
        Number of rays to process.

    Returns
    -------
    soln : array
        Solution to LP problem.

    See Also
    --------
    LP_solver_cylp_mp : Parent function.
    LP_solver_cylp : Single Process Solver.

    """
    from cylp.cy.CyClpSimplex import CyClpSimplex
    from cylp.py.modeling.CyLPModel import CyLPModel, CyLPArray

    n_gates = weights.shape[1] // 2
    n_rays = B_vectors.shape[0]
    soln = np.zeros([chunksize, n_gates])

    # import LP model in solver
    s = CyClpSimplex(model)

    # disable logging in multiprocessing anyway
    s.logLevel = 0

    i = 0
    for raynum in range(ray, ray + chunksize):
        # set new B_vector values for actual ray
        s.setRowLowerArray(np.squeeze(np.asarray(B_vectors[raynum])))
        # set new weights (objectives) for actual ray
        s.setObjectiveArray(np.squeeze(np.asarray(weights[raynum])))
        # solve with dual method, it is faster
        s.dual()
        # extract primal solution
        soln[i, :] = s.primalVariableSolution['x'][n_gates: 2 * n_gates]
        i = i + 1

    return soln


def LP_solver_cylp_mp(A_Matrix, B_vectors, weights, really_verbose=False,
                      proc=1):
    """
    Solve the Linear Programming problem given in Giangrande et al, 2012 using
    the CyLP module using multiple processes.

    Parameters
    ----------
    A_Matrix : matrix
        Row augmented A matrix, see :py:func:`construct_A_matrix`
    B_vectors : matrix
        Matrix containing B vectors, see :py:func:`construct_B_vectors`
    weights : array
        Weights.
    really_verbose : bool, optional
        True to print CLP messaging. False to suppress.
    proc : int, optional
        Number of worker processes.

    Returns
    -------
    soln : array
        Solution to LP problem.

    See Also
    --------
    LP_solver_cvxopt : Solve LP problem using the CVXOPT module.
    LP_solver_pyglpk : Solve LP problem using the PyGLPK module.
    LP_solver_cylp : Solve LP problem using the CyLP module using single
                     process.

    """
    from cylp.cy.CyClpSimplex import CyClpSimplex
    from cylp.py.modeling.CyLPModel import CyLPModel, CyLPArray
    import multiprocessing as mp

    n_gates = weights.shape[1] // 2
    n_rays = B_vectors.shape[0]
    soln = np.zeros([n_rays, n_gates])

    # Create CyLPModel and initialize it
    model = CyLPModel()
    G = np.matrix(A_Matrix)
    h = CyLPArray(np.empty(B_vectors.shape[1]))
    x = model.addVariable('x', G.shape[1])
    model.addConstraint(G * x >= h)
    c = CyLPArray(np.empty(weights.shape[1]))
    model.objective = c * x

    chunksize = int(n_rays/proc)
    # check if equal sized chunks can be distributed to worker processes
    if n_rays % chunksize != 0:
        print("Problem of %d rays cannot be split to %d worker processes!\n\r"
              "Fallback to 1 process!" % (n_rays, proc))
        chunksize = n_rays  # fall back to one process
        proc = 1

    print("Calculating with %d processes, %d rays per chunk" %
          (proc, chunksize))

    def worker(model, B_vectors, weights, ray, chunksize, out_q):
        """
        The worker function, invoked in a process.
        The results are placed in a dictionary that's pushed to a queue.
        """
        outdict = {}
        iray = int(ray/chunksize)
        outdict[iray] = solve_cylp(model, B_vectors, weights, ray, chunksize)
        out_q.put(outdict)

    # Queue for LP solutions
    out_q = mp.Queue()
    procs = []

    # fire off worker processes
    for raynum in range(0, n_rays, chunksize):
        p = mp.Process(target=worker, args=(
            model, B_vectors, weights, raynum, chunksize, out_q))
        procs.append(p)
        p.start()

    # collecting results
    resultdict = {}
    for raynum in range(0, n_rays, chunksize):
        resultdict.update(out_q.get())

    # Wait for all worker processes to finish
    for p in procs:
        p.join()

    # copy results in output array
    for raynum in range(0, int(n_rays / chunksize)):
        soln[raynum * chunksize:raynum * chunksize + chunksize, :] = (
            resultdict[raynum])

    # apply smoothing filter to output array
    soln = smooth_and_trim_scan(soln, window_len=5, window='sg_smooth')

    return soln


def LP_solver_cylp(A_Matrix, B_vectors, weights, really_verbose=False):
    """
    Solve the Linear Programming problem given in Giangrande et al, 2012 using
    the CyLP module.

    Parameters
    ----------
    A_Matrix : matrix
        Row augmented A matrix, see :py:func:`construct_A_matrix`
    B_vectors : matrix
        Matrix containing B vectors, see :py:func:`construct_B_vectors`
    weights : array
        Weights.
    really_verbose : bool, optional
        True to print CLP messaging. False to suppress.

    Returns
    -------
    soln : array
        Solution to LP problem.

    See Also
    --------
    LP_solver_cvxopt : Solve LP problem using the CVXOPT module.
    LP_solver_pyglpk : Solve LP problem using the PyGLPK module.

    """
    from cylp.cy.CyClpSimplex import CyClpSimplex
    from cylp.py.modeling.CyLPModel import CyLPModel, CyLPArray

    n_gates = weights.shape[1] // 2
    n_rays = B_vectors.shape[0]
    soln = np.zeros([n_rays, n_gates])

    # Create CyLPModel and initialize it
    model = CyLPModel()
    G = np.matrix(A_Matrix)
    h = CyLPArray(np.empty(B_vectors.shape[1]))
    x = model.addVariable('x', G.shape[1])
    model.addConstraint(G * x >= h)
    c = CyLPArray(np.squeeze(weights[0]))
    model.objective = c * x

    # import model in solver
    s = CyClpSimplex(model)
    # disable logging
    if not really_verbose:
        s.logLevel = 0

    for raynum in range(n_rays):

        # set new B_vector values for actual ray
        s.setRowLowerArray(np.squeeze(np.asarray(B_vectors[raynum])))
        # set new weights (objectives) for actual ray
        # solve with dual method, it is faster
        s.dual()
        # extract primal solution
        soln[raynum, :] = s.primalVariableSolution['x'][n_gates: 2 * n_gates]

    # apply smoothing filter on a per scan basis
    soln = smooth_and_trim_scan(soln, window_len=5, window='sg_smooth')
    return soln


def phase_proc_lp(radar, offset, debug=False, self_const=60000.0,
                  low_z=10.0, high_z=53.0, min_phidp=0.01, min_ncp=0.5,
                  min_rhv=0.8, fzl=4000.0, sys_phase=0.0,
                  overide_sys_phase=False, nowrap=None, really_verbose=False,
                  LP_solver='cylp', refl_field=None, ncp_field=None,
                  rhv_field=None, phidp_field=None, kdp_field=None,
                  unf_field=None, window_len=35, proc=1, coef=0.914):
    """
    Phase process using a LP method [1].

    Parameters
    ----------
    radar : Radar
        Input radar.
    offset : float
        Reflectivity offset in dBz.
    debug : bool, optional
        True to print debugging information.
    self_const : float, optional
        Self consistency factor.
    low_z : float, optional
        Low limit for reflectivity. Reflectivity below this value is set to
        this limit.
    high_z : float, optional
        High limit for reflectivity. Reflectivity above this value is set to
        this limit.
    min_phidp : float, optional
        Minimum Phi differential phase.
    min_ncp : float, optional
        Minimum normal coherent power.
    min_rhv : float, optional
        Minimum copolar coefficient.
    fzl : float, optional
        Maximum altitude.
    sys_phase : float, optional
        System phase in degrees.
    overide_sys_phase : bool, optional
        True to use `sys_phase` as the system phase. False will calculate a
        value automatically.
    nowrap : int or None, optional
        Gate number to begin phase unwrapping. None will unwrap all phases.
    really_verbose : bool, optional
        True to print LPX messaging. False to suppress.
    LP_solver : 'pyglpk' or 'cvxopt', 'cylp', or 'cylp_mp', optional
        Module to use to solve LP problem. Default is 'cylp'.
    refl_field, ncp_field, rhv_field, phidp_field, kdp_field : str, optional
        Name of field in radar which contains the horizonal reflectivity,
        normal coherent power, copolar coefficient, differential phase shift,
        and differential phase. A value of None for any of these parameters
        will use the default field name as defined in the Py-ART configuration
        file.
    unf_field : str, optional
        Name of field which will be added to the radar object which will
        contain the unfolded differential phase. Metadata for this field
        will be taken from the phidp_field. A value of None will use
        the default field name as defined in the Py-ART configuration file.
    window_len : int, optional
        Length of Sobel window applied to PhiDP field when prior to
        calculating KDP.
    proc : int, optional
        Number of worker processes, only used when `LP_solver` is 'cylp_mp'.
    coef : float, optional
        Exponent linking Z to KDP in self consistency. kdp=(10**(0.1z))*coef

    Returns
    -------
    reproc_phase : dict
        Field dictionary containing processed differential phase shifts.
    sob_kdp : dict
        Field dictionary containing recalculated differential phases.

    References
    ----------
    [1] Giangrande, S.E., R. McGraw, and L. Lei. An Application of
    Linear Programming to Polarimetric Radar Differential Phase Processing.
    J. Atmos. and Oceanic Tech, 2013, 30, 1716.

    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')
    if rhv_field is None:
        rhv_field = get_field_name('cross_correlation_ratio')
    if phidp_field is None:
        phidp_field = get_field_name('differential_phase')
    if kdp_field is None:
        kdp_field = get_field_name('specific_differential_phase')
    if unf_field is None:
        unf_field = get_field_name('unfolded_differential_phase')

    # prepare reflectivity field
    refl = copy.deepcopy(radar.fields[refl_field]['data']) + offset
    is_low_z = (refl) < low_z
    is_high_z = (refl) > high_z
    refl[np.where(is_high_z)] = high_z
    refl[np.where(is_low_z)] = low_z
    z_mod = refl

    # unfold Phi_DP
    if debug:
        print('Unfolding')
    my_unf = get_phidp_unf(radar, ncp_lev=min_ncp, rhohv_lev=min_rhv,
                           debug=debug, ncpts=2, doc=None,
                           sys_phase=sys_phase, nowrap=nowrap,
                           overide_sys_phase=overide_sys_phase,
                           refl_field=refl_field, ncp_field=ncp_field,
                           rhv_field=rhv_field, phidp_field=phidp_field)
    my_new_ph = copy.deepcopy(radar.fields[phidp_field])
    my_unf[:, -1] = my_unf[:, -2]
    my_new_ph['data'] = my_unf
    radar.fields.update({unf_field: my_new_ph})

    phidp_mod = copy.deepcopy(radar.fields[unf_field]['data'])
    phidp_neg = phidp_mod < min_phidp
    phidp_mod[np.where(phidp_neg)] = min_phidp

    # process
    proc_ph = copy.deepcopy(radar.fields[phidp_field])
    proc_ph['data'] = phidp_mod
    St_Gorlv_differential_5pts = [-.2, -.1, 0, .1, .2]
    for sweep in range(len(radar.sweep_start_ray_index['data'])):
        if debug:
            print("Doing ", sweep)
        end_gate, start_ray, end_ray = det_process_range(
            radar, sweep, fzl, doc=15)
        start_gate = 0

        A_Matrix = construct_A_matrix(
            len(radar.range['data'][start_gate:end_gate]),
            St_Gorlv_differential_5pts)

        B_vectors = construct_B_vectors(
            phidp_mod[start_ray:end_ray, start_gate:end_gate],
            z_mod[start_ray:end_ray, start_gate:end_gate],
            St_Gorlv_differential_5pts, dweight=self_const,
            coef=coef)

        weights = np.ones(
            phidp_mod[start_ray:end_ray, start_gate:end_gate].shape)

        nw = np.bmat([weights, np.zeros(weights.shape)])

        if LP_solver == 'pyglpk':
            mysoln = LP_solver_pyglpk(A_Matrix, B_vectors, nw,
                                      really_verbose=really_verbose)
        elif LP_solver == 'cvxopt':
            mysoln = LP_solver_cvxopt(A_Matrix, B_vectors, nw)
        elif LP_solver == 'cylp':
            mysoln = LP_solver_cylp(A_Matrix, B_vectors, nw,
                                    really_verbose=really_verbose)
        elif LP_solver == 'cylp_mp':
            mysoln = LP_solver_cylp_mp(A_Matrix, B_vectors, nw,
                                       really_verbose=really_verbose,
                                       proc=proc)
        else:
            raise ValueError('unknown LP_solver:' + LP_solver)

        proc_ph['data'][start_ray:end_ray, start_gate:end_gate] = mysoln

    last_gates = proc_ph['data'][start_ray:end_ray, -16]
    proc_ph['data'][start_ray:end_ray, -16:] = \
        np.meshgrid(np.ones([16]), last_gates)[1]
    proc_ph['valid_min'] = 0.0          # XXX is this correct?
    proc_ph['valid_max'] = 400.0        # XXX is this correct?

    # prepare output
    sobel = 2. * np.arange(window_len)/(window_len - 1.0) - 1.0
    sobel = sobel/(abs(sobel).sum())
    sobel = sobel[::-1]
    gate_spacing = (radar.range['data'][1] - radar.range['data'][0]) / 1000.
    kdp = (scipy.ndimage.filters.convolve1d(proc_ph['data'], sobel, axis=1) /
           ((window_len / 3.0) * 2.0 * gate_spacing))

    # copy the KDP metadata from existing field or create anew
    if kdp_field in radar.fields:
        sob_kdp = copy.deepcopy(radar.fields[kdp_field])
    else:
        sob_kdp = get_metadata(kdp_field)

    sob_kdp['data'] = kdp
    sob_kdp['_FillValue'] = get_fillvalue()

    return proc_ph, sob_kdp


def phase_proc_lp_gf(radar, gatefilter=None, debug=False, self_const=60000.0,
                     low_z=10.0, high_z=53.0, min_phidp=0.01, fzl=4000.0,
                     system_phase=None, nowrap=None, really_verbose=False,
                     LP_solver='cylp', refl_field=None, phidp_field=None,
                     kdp_field=None, unf_field=None, window_len=35, proc=1,
                     coef=0.914, ncpts=None, first_gate_sysp=None, offset=0.0,
                     doc=0):
    """
    Phase process using a LP method [1] using Py-ART's Gatefilter.

    Parameters
    ----------
    radar : Radar
        Input radar.
    gatefilter : Gatefilter, optional
        Py-ART gatefilter object indicating where processing should be
        carried out
    debug : bool, optional
        True to print debugging information.
    self_const : float, optional
        Self consistency factor.
    low_z : float, optional
        Low limit for reflectivity. Reflectivity below this value is set to
        this limit.
    high_z : float, optional
        High limit for reflectivity. Reflectivity above this value is set to
        this limit.
    fzl : float, optional
        Maximum altitude.
    system_phase : float, optional
        System phase in degrees.
    nowrap : int or None, optional
        Gate number to begin phase unwrapping. None will unwrap all phases.
    really_verbose : bool, optional
        True to print LPX messaging. False to suppress.
    LP_solver : 'pyglpk' or 'cvxopt', 'cylp', or 'cylp_mp', optional
        Module to use to solve LP problem. Default is 'cylp'.
    refl_field, ncp_field, rhv_field, phidp_field, kdp_field : str, optional
        Name of field in radar which contains the horizonal reflectivity,
        normal coherent power, copolar coefficient, differential phase shift,
        and differential phase. A value of None for any of these parameters
        will use the default field name as defined in the Py-ART configuration
        file.
    unf_field : str, optional
        Name of field which will be added to the radar object which will
        contain the unfolded differential phase. Metadata for this field
        will be taken from the phidp_field. A value of None will use
        the default field name as defined in the Py-ART configuration file.
    window_len : int, optional
        Length of Sobel window applied to PhiDP field when prior to
        calculating KDP.
    proc : int, optional
        Number of worker processes, only used when `LP_solver` is 'cylp_mp'.
    coef : float, optional
        Exponent linking Z to KDP in self consistency. kdp=(10**(0.1z))*coef
    ncpts : int, optional
        Minimum number of points in a ray. Regions within a ray smaller than
        this or beginning before this gate number are excluded from unfolding.
    offset : float, optional
        Reflectivity offset to add in dBz.
    doc : int, optional
        Number of gates to "doc" off the end of a ray.

    Returns
    -------
    reproc_phase : dict
        Field dictionary containing processed differential phase shifts.
    sob_kdp : dict
        Field dictionary containing recalculated differential phases.

    References
    ----------
    [1] Giangrande, S.E., R. McGraw, and L. Lei. An Application of
    Linear Programming to Polarimetric Radar Differential Phase Processing.
    J. Atmos. and Oceanic Tech, 2013, 30, 1716.

    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if phidp_field is None:
        phidp_field = get_field_name('differential_phase')
    if kdp_field is None:
        kdp_field = get_field_name('specific_differential_phase')
    if unf_field is None:
        unf_field = get_field_name('unfolded_differential_phase')

    # if there is no gatefilter included, include all

    if gatefilter is None:
        gatefilter = GateFilter(radar)
        gatefilter.include_all()

    # prepare reflectivity field
    refl = copy.deepcopy(radar.fields[refl_field]['data']) + offset
    is_low_z = (refl) < low_z
    is_high_z = (refl) > high_z
    refl[np.where(is_high_z)] = high_z
    refl[np.where(is_low_z)] = low_z
    z_mod = refl

    # unfold Phi_DP
    if debug:
        print('Unfolding')

    my_unf = get_phidp_unf_gf(radar, gatefilter,
                              debug=debug, ncpts=ncpts,
                              sys_phase=system_phase, nowrap=nowrap,
                              phidp_field=phidp_field,
                              first_gate_sysp=first_gate_sysp)

    my_new_ph = copy.deepcopy(radar.fields[phidp_field])
    my_unf[:, -1] = my_unf[:, -2]
    my_new_ph['data'] = my_unf
    radar.fields.update({unf_field: my_new_ph})

    phidp_mod = copy.deepcopy(radar.fields[unf_field]['data'])
    phidp_neg = phidp_mod < min_phidp
    phidp_mod[np.where(phidp_neg)] = min_phidp

    # process
    proc_ph = copy.deepcopy(radar.fields[phidp_field])
    proc_ph['data'] = phidp_mod
    St_Gorlv_differential_5pts = [-.2, -.1, 0, .1, .2]

    for sweep in range(len(radar.sweep_start_ray_index['data'])):
        if debug:
            print("Doing ", sweep)

        end_gate, start_ray, end_ray = det_process_range(
            radar, sweep, fzl, doc=doc)

        start_gate = 0

        A_Matrix = construct_A_matrix(
            len(radar.range['data'][start_gate:end_gate]),
            St_Gorlv_differential_5pts)

        B_vectors = construct_B_vectors(
            phidp_mod[start_ray:end_ray, start_gate:end_gate],
            z_mod[start_ray:end_ray, start_gate:end_gate],
            St_Gorlv_differential_5pts, dweight=self_const,
            coef=coef)

        weights = np.ones(
            phidp_mod[start_ray:end_ray, start_gate:end_gate].shape)

        nw = np.bmat([weights, np.zeros(weights.shape)])

        if LP_solver == 'pyglpk':
            mysoln = LP_solver_pyglpk(A_Matrix, B_vectors, nw,
                                      really_verbose=really_verbose)
        elif LP_solver == 'cvxopt':
            mysoln = LP_solver_cvxopt(A_Matrix, B_vectors, nw)
        elif LP_solver == 'cylp':
            mysoln = LP_solver_cylp(A_Matrix, B_vectors, nw,
                                    really_verbose=really_verbose)
        elif LP_solver == 'cylp_mp':
            mysoln = LP_solver_cylp_mp(A_Matrix, B_vectors, nw,
                                       really_verbose=really_verbose,
                                       proc=proc)
        else:
            raise ValueError('unknown LP_solver:' + LP_solver)

        proc_ph['data'][start_ray:end_ray, start_gate:end_gate] = mysoln

    last_gates = proc_ph['data'][start_ray:end_ray, -16]
    proc_ph['data'][start_ray:end_ray, -16:] = \
        np.meshgrid(np.ones([16]), last_gates)[1]
    proc_ph['valid_min'] = 0.0  # XXX is this correct?
    proc_ph['valid_max'] = 400.0  # XXX is this correct?

    # prepare output
    sobel = 2. * np.arange(window_len) / (window_len - 1.0) - 1.0
    sobel = sobel / (abs(sobel).sum())
    sobel = sobel[::-1]
    gate_spacing = (radar.range['data'][1] - radar.range['data'][0]) / 1000.
    kdp = (scipy.ndimage.filters.convolve1d(proc_ph['data'], sobel, axis=1) /
           ((window_len / 3.0) * 2.0 * gate_spacing))

    # copy the KDP metadata from existing field or create anew
    if kdp_field in radar.fields:
        sob_kdp = copy.deepcopy(radar.fields[kdp_field])
    else:
        sob_kdp = get_metadata(kdp_field)

    sob_kdp['data'] = kdp
    sob_kdp['_FillValue'] = get_fillvalue()

    return proc_ph, sob_kdp


def get_phidp_unf_gf(radar, gatefilter, debug=False, ncpts=2, sys_phase=None,
                     nowrap=None, phidp_field=None, first_gate_sysp=None):
    """
    Get Unfolded Phi differential phase in areas not gatefiltered.

    Parameters
    ----------
    radar : Radar
        The input radar.
    gatefilter : GateFilter
        Only apply on areas incuded in the gatefilter
    debug : bool, optioanl
        True to print debugging information, False to supress printing.
    ncpts : int, optional
        Minimum number of points in a ray. Regions within a ray smaller than
        this or beginning before this gate number are excluded from
        calculations.
    sys_phase : float, optional
        System phase overide.
    nowrap : int or None, optional
        Gate number where unwrapping should begin. `None` will unwrap all
        gates.
    refl_field ncp_field, rhv_field, phidp_field : str, optional
        Field names within the radar object which represent the horizontal
        reflectivity, normal coherent power, the copolar coefficient, and the
        differential phase shift. A value of None for any of these parameters
        will use the default field name as defined in the Py-ART
        configuration file.

    Returns
    -------
    cordata : array
        Unwrapped phi differential phase.

    """
    # parse the field parameters
    if phidp_field is None:
        phidp_field = get_field_name('differential_phase')

    t = time()
    my_phidp = radar.fields[phidp_field]['data']
    if sys_phase is not None:
        system_zero = sys_phase
    else:
        system_zero = det_sys_phase_gf(radar, gatefilter,
                                       phidp_field=phidp_field,
                                       first_gate=first_gate_sysp)
        if system_zero is None:
            system_zero = -135
    cordata = np.zeros(my_phidp.shape, dtype=float)
    all_non_meteo = gatefilter.gate_excluded
    for radial in range(my_phidp.shape[0]):

        # deterimine non meteo from gatefilter CHANGED
        notmeteo = all_non_meteo[radial, :]
        x_ma = ma.masked_where(notmeteo, my_phidp[radial, :])
        try:
            ma.notmasked_contiguous(x_ma)
            for slc in ma.notmasked_contiguous(x_ma):
                # so trying to get rid of clutter and small things that
                # should not add to phidp anyway
                if slc.stop - slc.start < ncpts or slc.start < ncpts:
                    x_ma.mask[slc.start - 1:slc.stop + 1] = True
            c = 0
        except TypeError:  # non sequence, no valid regions
            c = 1  # ie do nothing
            x_ma.mask = True
        except AttributeError:
            c = 1  # also do nothing
            x_ma.mask = True
        if nowrap is not None:
            # Start the unfolding a bit later in order to avoid false
            # jumps based on clutter
            unwrapped = copy.deepcopy(x_ma)
            end_unwrap = unwrap_masked(x_ma[nowrap::], centered=False)
            unwrapped[nowrap::] = end_unwrap
        else:
            unwrapped = unwrap_masked(x_ma, centered=False)

        # end so no clutter expected
        system_max = unwrapped[np.where(np.logical_not(
            notmeteo))][-10:-1].mean() - system_zero
        unwrapped_fixed = np.zeros(len(x_ma), dtype=float)
        based = unwrapped-system_zero
        based[0] = 0.0
        notmeteo[0] = False
        based[-1] = system_max
        notmeteo[-1] = False
        unwrapped_fixed[np.where(np.logical_not(based.mask))[0]] = \
            based[np.where(np.logical_not(based.mask))[0]]
        if len(based[np.where(np.logical_not(based.mask))[0]]) > 11:
            unwrapped_fixed[np.where(based.mask)[0]] = \
                np.interp(np.where(based.mask)[0],
                          np.where(np.logical_not(based.mask))[0],
                          smooth_and_trim(based[np.where(
                              np.logical_not(based.mask))[0]]))
        else:
            unwrapped_fixed[np.where(based.mask)[0]] = \
                np.interp(np.where(based.mask)[0],
                          np.where(np.logical_not(based.mask))[0],
                          based[np.where(np.logical_not(based.mask))[0]])
        if c != 1:
            cordata[radial, :] = unwrapped_fixed
        else:
            cordata[radial, :] = np.zeros(my_phidp.shape[1])
    if debug:
        print("Exec time: ", time() - t)
    return cordata


def det_sys_phase_gf(radar, gatefilter, phidp_field=None, first_gate=30.):
    """
    Determine the system phase.

    Parameters
    ----------
    radar : Radar
        Radar object for which to determine the system phase.
    gatefilter : Gatefilter
        Gatefilter object highlighting valid gates.
    phidp_field : str, optional
        Field name within the radar object which represents
        differential phase shift. A value of None will use the default
        field name as defined in the Py-ART configuration file.
    first_gate : int, optional
        Gate index for where to being applying the gatefilter.

    Returns
    -------
    sys_phase : float or None
        Estimate of the system phase. None is not estimate can be made.

    """
    # parse the field parameters
    if phidp_field is None:
        phidp_field = get_field_name('differential_phase')
    phidp = radar.fields[phidp_field]['data'][:, first_gate:]
    last_ray_idx = radar.sweep_end_ray_index['data'][0]
    is_meteo = gatefilter.gate_included[:, first_gate:]
    return _det_sys_phase_gf(phidp, last_ray_idx, is_meteo)


def _det_sys_phase_gf(phidp, last_ray_idx, radar_meteo):
    """ Determine the system phase, see :py:func:`det_sys_phase`. """
    good = False
    phases = []
    for radial in range(last_ray_idx + 1):
        meteo = radar_meteo[radial, :]
        mpts = np.where(meteo)
        if len(mpts[0]) > 25:
            good = True
            msmth_phidp = smooth_and_trim(phidp[radial, mpts[0]], 9)
            phases.append(msmth_phidp[0:25].min())
    if not good:
        return None
    return np.median(phases)
