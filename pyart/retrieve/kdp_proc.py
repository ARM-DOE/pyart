"""
pyart.retrieve.kdp_proc
=======================

Module for retrieving specific differential phase (KDP) from radar total
differential phase (PSIDP) measurements. Total differential phase is a function
of propagation differential phase (PHIDP), backscatter differential phase
(DELTAHV), and the system phase offset.

.. autosummary::
    :toctree: generated/

    kdp_maesaka
    boundary_conditions_maesaka
    _cost_maesaka
    _jac_maesaka
    _forward_reverse_phidp
    _parse_range_resolution

"""

import time
import numpy as np

from scipy import optimize, stats

from . import _kdp_proc
from ..config import get_field_name, get_metadata, get_fillvalue

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
        psidp_o[ray, idx_far[ray]+1:] = np.ma.masked

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
    if kwargs.get('check_outliers', True):

        # do not include missing values in the analysis
        phi_near_valid = phi_near[phi_near != 0.0]

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
