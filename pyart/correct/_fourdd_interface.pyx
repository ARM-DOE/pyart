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
                         vad_time):
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
    vad_time : int
        Time of sounding in YYDDDHHMM format.  Note that this is not
        acutally used in the underlying C code so the value can be anything.

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
    _fourdd_h.firstGuessNoRead(soundVolume._Volume, MISSINGVEL,
                               <float *> hc.data, <float *> sc.data,
                               <float *> dc.data, <int> len(hc),
                               vad_time, &success)
    return success, soundVolume


cpdef fourdd_dealias(_RslVolume radialVelVolume,
                     _RslVolume lastVelVolume,
                     _RslVolume soundVolume,
                     _RslVolume DZvolume,
                     prep, filt, debug):
    """
    fourdd_dealias(DZvolume, radialVelVolume, lastVelVolume hc, sc, dc,
                   vad_time, prep, filt)

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
    DZvolume : _RslVolume or None
        Reflectivity to use when thresholding is selected.  This parameter
        is not used when prep is 0 and can then be set to None.
    prep : int
        Flag controlling thresholding of DZvolume, 1 = yes, 0 = no.
    filt : int
        Flag controlling Bergen and Albers filter, 1 = yes, 0 = no.
    debug : bool
        True to return RSL Volume objects for debugging:
        usuccess, radialVelVolume, lastVelVolume, soundVolume, DZvolume,
        unfoldedVolume

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
    cdef unsigned short usuccess = 0

    if lastVelVolume is None and soundVolume is None:
        raise ValueError('lastVelVolume or soundVolume must be defined')

    # The following closely follows that in FourDD.c starting at line 142

    # copy the radial velocity data to unfoldedVolume
    unfoldedVolume = _rsl_interface.copy_volume(radialVelVolume)

    if prep:
        if DZvolume is None:
            raise ValueError('DZvolume must be defined if prep is True')
        # remove any bins where reflectivity is missing or outside
        # the accepted interval.
        _fourdd_h.prepVolume(DZvolume._Volume, unfoldedVolume._Volume,
                             MISSINGVEL)

    # unfold the velocity fields in unfoldedVolume
    if lastVelVolume is None:   # only soundVolume
        _fourdd_h.unfoldVolume(unfoldedVolume._Volume, soundVolume._Volume,
                               NULL,
                               MISSINGVEL, filt, &usuccess)
    elif soundVolume is None:   # only lastVelVolume
        _fourdd_h.unfoldVolume(unfoldedVolume._Volume, NULL,
                               lastVelVolume._Volume,
                               MISSINGVEL, filt, &usuccess)
    else:   # both soundVolume and lastVelVolume
        _fourdd_h.unfoldVolume(unfoldedVolume._Volume, soundVolume._Volume,
                               lastVelVolume._Volume,
                               MISSINGVEL, filt, &usuccess)

    if debug:
        return (usuccess, radialVelVolume, lastVelVolume, soundVolume,
                DZvolume, unfoldedVolume)

    data = unfoldedVolume.get_data()
    return usuccess, data
