"""
pyart.correct._fourdd_interface
===============================

Cython wrapper around the University of Washington FourDD algorithm.

.. autosummary::
    :toctree: generated/

    create_soundvolume
    fourdd_dealias

"""

cimport _fourdd_h
cimport numpy as np
from pyart.io._rsl_interface cimport _RslVolume

from ..io import _rsl_interface


cpdef create_soundvolume(radialVelVolume,
                         np.ndarray[np.float32_t, ndim=1] hc,
                         np.ndarray[np.float32_t, ndim=1] sc,
                         np.ndarray[np.float32_t, ndim=1] dc,
                         maxshear=0.05, sign=1):
    """
    Create a RSL Volume containing sounding data.

    Parameters
    ----------
    radialVelVolume : _RslVolume
        Radial velocities which will be dealiased, shape used to create
        soundvolume.
    hc : ndarray
        Sounding heights in meters.  Must be a contiguous one-dimensional
        float32 array.
    sc : ndarray
        Sounding wind speed in m/s.  Must be a contiguous one-dimensional
        float32 array.
    dc : ndarray
        Sounding wind direction in degrees.  Must be a contiguous
        one-dimensional float32 array.
    maxshear : float
        Maximum vertical shear which will be incorperated into the created
        volume.
    sign : int
        Sign convention which the radial velocities in the created volume
        will follow.  A value of 1 represents when positive values
        velocities are towards the radar, -1 represents when negative
        velocities are towards the radar.

    Returns
    -------
    usuccess : int
        Flag indicating if loading of data was successful, 1 = yes, 0 = no.
    soundvolume : _RslVolume
        RslVolume containing sounding data.

    """

    cdef _RslVolume soundVolume
    cdef unsigned short success = 0
    cdef float MISSINGVEL = 131072.0
    soundVolume = _rsl_interface.copy_volume(radialVelVolume)
    success = _fourdd_h.sounding_to_volume(
        soundVolume._Volume, MISSINGVEL,
        <float *> hc.data, <float *> sc.data, <float *> dc.data,
        <int> len(hc), maxshear, sign)
    return success, soundVolume


cpdef fourdd_dealias(
    _RslVolume radialVelVolume, _RslVolume lastVelVolume,
    _RslVolume soundVolume, filt,
    compthresh=0.25, compthresh2=0.49,
    thresh=0.4, ckval=1.0, stdthresh=0.8, epsilon=0.00001, maxcount=10,
    pass2=1, rm=0, proximity=5, mingood=5, ba_mincount=5, ba_edgecount=3,
    debug=False):
    """
    fourdd_dealias(
        radialVelVolume, lastVelVolume, soundVolume, filt,
        compthresh=0.25, compthresh2=0.49, thresh=0.4,
        epsilon=0.00001, ckval=1.0, stdthresh=0.8, maxcount=10, pass2=1,
        rm=0, proximity=5, mingood=5, ba_mincount=5, ba_edgecount=3,
        debug=False)

    Dealias using the FourDD algorithm.

    Parameters
    ----------
    radialVelVolume : _RslVolume
        Radial velocities which will be dealiased.
    lastVelVolume : _RslVolume or None
        Radial velocities from a previously dealiased radar volume. For best
        results, this radar should represent the previous volume scan in time.
        If the last velocity volume is unavailable, set this to None.
    soundVolume : _RslVolume or None
        Volume created from sounding data.  If unavailable, set this to None.
        soundVolume and lastVelVolume cannot both be None.
    filt : int
        Flag controlling Bergen and Albers filter, 1 = yes, 0 = no.

    Other Parameters
    ----------------
    compthresh : float
        Fraction of the Nyquist velocity to use as a threshold when performing
        continity (initial) dealiasing.  Velocities differences above this
        threshold will not be marked as gate from which to begin unfolding
        during spatial dealiasing.
    compthresh2 : float
        The same as compthresh but the value used during the second pass of
        dealasing.  This second pass is only performed in both a sounding
        and last volume are provided.
    thresh : float
        Fraction of the Nyquist velocity to use as a threshold when performing
        spatial dealiasing.  Horizontally adjacent gates with velocities above
        this theshold will count against assigning the gate in question the
        velocity value being tested.
    ckval : float
        When the absolute value of the velocities are below this value they
        will not be marked as gates from which to begin unfolding during
        spatial dealiasing.
    stdthresh : float
       Fraction of the Nyquist velocity to use as a standard deviation
       threshold in the window dealiasing portion of the algorithm.
    epsilon : float
        Difference used when comparing a value to missing value, changing this
        from the default is not recommended.
    maxcount : int
        Maximum allowed number of fold allowed when unfolding velocities.
    pass2 : int
        Controls weather unfolded gates should be removed (a value of 0)
        or retained for unfolding during the second pass (a value of 1) when
        both a sounding volume and last volume are provided.
    rm : int
        Determines what should be done with gates that are left unfolded
        after the first pass of dealiasing.  A value of 1 will remove these
        gates, a value of 0 sets these gates to their initial velocity.  If
        both a sounding volume and last volume are provided this parameter is
        ignored.
    proximity : int
        Number of gates and rays to include of either side of the current gate
        during window dealiasing.  This value may be doubled in cases where
        a standard sized window does not capture a sufficient number of
        good valued gates.
    mingood : int
        Number of good valued gates required within the window before the
        current gate will be unfolded.
    ba_mincount : int
        Number of neighbors required during Bergen and Albers filter for
        a given gate to be included, must be between 1 and 8, 5 recommended.
    ba_edgecount : int
        Same as ba_mincount but used at ray edges, must be between 1 and 5,
        3 recommended.
    debug : bool
        True to return RSL Volume objects for debugging:
        usuccess, radialVelVolume, lastVelVolume, soundVolume, unfoldedVolume

    Returns
    -------
    usuccess : int
        Flag indicating if the unfolding was successful, 1 = yes, 0 = no.
    data : np.ndarray
        Array of unfolded velocities.

    References
    ----------
    C. N. James and R. A Houze Jr, A Real-Time Four-Dimensional Doppler
    Dealising Scheme, Journal of Atmospheric and Oceanic Technology, 2001, 18,
    1674.

    """
    cdef _RslVolume unfoldedVolume
    cdef float MISSINGVEL = 131072.0

    if lastVelVolume is None and soundVolume is None:
        raise ValueError('lastVelVolume or soundVolume must be defined')

    # The following closely follows that in FourDD.c starting at line 142

    # copy the radial velocity data to unfoldedVolume
    unfoldedVolume = _rsl_interface.copy_volume(radialVelVolume)

    # unfold the velocity fields in unfoldedVolume
    if lastVelVolume is None:   # only soundVolume
        usuccess  = _fourdd_h.dealias_fourdd(
            unfoldedVolume._Volume, soundVolume._Volume, NULL,
            MISSINGVEL, compthresh, compthresh2, thresh,
            ckval, stdthresh, epsilon,
            maxcount, pass2, rm, proximity, mingood,
            filt, ba_mincount, ba_edgecount)
    elif soundVolume is None:   # only lastVelVolume
        usuccess = _fourdd_h.dealias_fourdd(
            unfoldedVolume._Volume, NULL, lastVelVolume._Volume,
            MISSINGVEL, compthresh, compthresh2, thresh,
            ckval, stdthresh, epsilon,
            maxcount, pass2, rm, proximity, mingood,
            filt, ba_mincount, ba_edgecount)
    else:   # both soundVolume and lastVelVolume
        usuccess = _fourdd_h.dealias_fourdd(
            unfoldedVolume._Volume, soundVolume._Volume, lastVelVolume._Volume,
            MISSINGVEL, compthresh, compthresh2, thresh,
            ckval, stdthresh, epsilon,
            maxcount, pass2, rm, proximity, mingood,
            filt, ba_mincount, ba_edgecount)
    if debug:
        return (usuccess, radialVelVolume, lastVelVolume, soundVolume,
                unfoldedVolume)

    data = unfoldedVolume.get_data()
    return usuccess, data
