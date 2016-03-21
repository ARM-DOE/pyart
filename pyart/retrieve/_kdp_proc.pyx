"""
pyart.retrieve._kdp_proc
========================

Cython routines for specific differential phase retrievals.

.. autosummary::
    :toctree: generated/

    lowpass_maesaka_term
    lowpass_maesaka_jac

"""

# Necessary and/or potential future improvements to this module:
#
# * High order finite difference schemes still need to be added. Note that the
#   added value of using a higher order scheme has yet to be determined and may
#   be fruitless in the end.


cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def lowpass_maesaka_term(
        double[:, ::1] k, dr, finite_order, double[:, ::1] d2kdr2):
    """
    Compute the filter term.

    Compute the low-pass filter term found in Maesaka et al. (2012). This term
    represents the second-order derivative of the control variable k with
    respect to range. This subroutine does not currently support radars with
    variable range resolution.

    Parameters
    ----------
    k : 2D array of float64
        Control variable k defined in Maesaka et al. (2012). This variable is
        proportional to the square root of specific differential phase.
    dr : float
        The range resolution in meters.
    finite_order : str, 'low' or 'high'
        The finite difference accuracy to use when computing the second-order
        range derivative of the control variable k.
    d2kdr2 : 2D array of float64
        Second-order derivative of k with respect to range. Updated in place.

    """

    # define local variables
    cdef int nr, ng, g, r
    cdef double dr2

    nr = k.shape[0]
    ng = k.shape[1]
    dr2 = dr**2.

    # Use a low order finite difference scheme to compute the second-order
    # range derivative
    if finite_order == 'low':
        for r in range(nr):
            for g in range(ng):
                # For interior range gates, i.e. g = [1, ng-1]
                # use a centered difference scheme where p = 2.

                # When at ray boundaires, i.e., g = 0 or ng-1 use either a
                # forward or backward difference scheme where p = 1 depending
                # on which boundary

                # Computing --> d2k/dr2
                if g > 0 and g < ng-1:
                    d2kdr2[r, g] = (k[r, g+1] - 2.*k[r, g] + k[r, g-1]) / dr2
                elif g == 0:
                    d2kdr2[r, g] = (k[r, g] - 2.*k[r, g+1] + k[r, g+2]) / dr2
                else:
                    d2kdr2[r, g] = (k[r, g] - 2.*k[r, g-1] + k[r, g-2]) / dr2
    else:
        raise ValueError("Invalid finite_order")
    return


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def lowpass_maesaka_jac(
        double[:, ::1] d2kdr2, dr, double Clpf, finite_order,
        double[:, ::1] dJlpfdk):
    """
    Compute the Jacobian of the filter cost functional.

    Compute the Jacobian of the low-pass filter cost functional similar to
    equation (18) in Maesaka et al. (2012). This function does not currently
    support radars with variable range resolution.

    Parameters
    ----------
    d2kdr2 : 2D array of float64
       Second-order derivative of the control variable k with respect to range.
       The control variable k is proportional to the square root of specific
       differential phase.
    dr : float
       The range resolution in meters.
    Clpf : float
       The low-pass filter (radial smoothness) constraint weight.
    finite_order :  str, 'low' or 'high'
       The finite difference accuracy used to compute the second-order range
       derivative of the control variable k.
    dJlpfdk : 2D array of float64
       The Jacobian of the low-pass filter cost functional with respect to the
       control variable k.  Updated in place.

    """

    # Define local variables
    cdef int nr, ng, g, r
    cdef double dr2

    nr = d2kdr2.shape[0]
    ng = d2kdr2.shape[1]

    # The low-pass filter cost is defined as,
    # Jlpf = 0.5 * Clpf * sum[ (d2k/dr2)**2 ] ,
    # where the sum is over all range gates for all rays.
    dr2 = dr**2

    # The Jacobian of Jlpf when a low finite order has been used to compute the
    # second-order range derivative of the control variable k
    if finite_order == 'low':
        for r in range(nr):
            for g in range(ng):
                if g > 2 and g < ng - 3:
                    dJlpfdk[r, g] = Clpf * (
                        d2kdr2[r, g-1] - 2.*d2kdr2[r, g] +
                        d2kdr2[r, g+1]) / dr2
                elif g == 2:
                    dJlpfdk[r, g] = Clpf * (
                        d2kdr2[r, g-2] + d2kdr2[r, g-1] - 2.*d2kdr2[r, g] +
                        d2kdr2[r, g+1]) / dr2
                elif g == 1:
                    dJlpfdk[r, g] = Clpf * (
                        d2kdr2[r, g+1] - 2.*d2kdr2[r, g] -
                        2.*d2kdr2[r, g-1]) / dr2
                elif g == 0:
                    dJlpfdk[r, g] = Clpf * (
                        d2kdr2[r, g] + d2kdr2[r, g+1]) / dr2
                elif g == ng - 3:
                    dJlpfdk[r, g] = Clpf * (
                        d2kdr2[r, g+2] + d2kdr2[r, g+1] - 2.*d2kdr2[r, g] +
                        d2kdr2[r, g-1]) / dr2
                elif g == ng - 2:
                    dJlpfdk[r, g] = Clpf * (
                        d2kdr2[r, g-1] - 2.*d2kdr2[r, g] -
                        2.*d2kdr2[r, g+1]) / dr2
                else:
                    dJlpfdk[r, g] = Clpf * (
                        d2kdr2[r, g] + d2kdr2[r, g-1]) / dr2
    else:
        raise ValueError("Invalid finite_order")
    return
