"""
pyart.correct._fourdd_interface
===============================

Cython wrapper around the University of Washington FourDD algorithm.

.. autosummary::
    :toctree: generated/

    fourdd_dealias

"""

cimport _fourdd_h
cimport numpy as np
from pyart.io._rsl_interface cimport _RslVolume

from ..io import _rsl_interface


cpdef fourdd_dealias(_RslVolume DZvolume,
                     _RslVolume radialVelVolume,
                     _RslVolume lastVelVolume,
                     np.ndarray[np.float32_t, ndim=1] hc,
                     np.ndarray[np.float32_t, ndim=1] sc,
                     np.ndarray[np.float32_t, ndim=1] dc,
                     vad_time, prep, filt, debug):
    """
    fourdd_dealias(DZvolume, radialVelVolume, lastVelVolume hc, sc, dc,
                   vad_time, prep, filt)

    Dealias using the FourDD algorithm.

    Parameters
    ----------
    DZvolume : _RslVolume
        Reflectivity to use when thresholding is selected.
    radialVelVolume : _RslVolume
        Radial velocities which will be dealiased.
    lastVelVolume : _RslVolume
        Radial velocities from a previously dealiased radar volume. For best
        results, this radar should represent the previous volume scan in time.
        If the last velocity volume is unavailable, set this to None.
    hc : np.ndarray
        Sounding heights in meters.  Must be a contiguous one-dimensional
        float32 array.
    sc : np.ndarray
        Sounding wind speed in m/s.  Must be a contiguous one-dimensional
        float32 array.
    dc : np.ndarray
        Sounding wind direction in degrees.  Must be a contiguous
        one-dimensional float32 array.
    vad_time : int
        Time of sounding in YYDDDHHMM format.
    prep : int
        Flag controlling thresholding of DZvolume, 1 = yes, 0 = no.
    filt : int
        Flag controlling Bergen and Albers filter, 1 = yes, 0 = no.
    debug : bool
        True to return debugging RSL Volume objects for debugging:
        usuccess, DZvolume, radialVelVolume, unfoldedVolume, sondVolume


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
    
    # TODO: version which does not need DZvolume, i.e. prep = 0
    # See FourDD.c for addition details
    
    cdef _RslVolume unfoldedVolume
    cdef _RslVolume sondVolume
    cdef float MISSINGVEL = 131072.0
    cdef unsigned short success = 0
    cdef unsigned short usuccess = 0

    unfoldedVolume = _rsl_interface.copy_volume(radialVelVolume)
    sondVolume = _rsl_interface.copy_volume(radialVelVolume)

    # may not always be needed...
    _fourdd_h.firstGuessNoRead(sondVolume._Volume, MISSINGVEL,
                        <float *> hc.data, <float *> sc.data,
                        <float *> dc.data, <int> len(hc),
                        vad_time, &success)

    if success != 1:
        raise ValueError('Error with sounding data')

    # prepare the velocity volume
    # any gates where reflectivity is bad "should" be removed
    if prep == 1:
        _fourdd_h.prepVolume(DZvolume._Volume, unfoldedVolume._Volume,
                             MISSINGVEL)
        
    # dealias
    if lastVelVolume is not None:
        _fourdd_h.unfoldVolume(unfoldedVolume._Volume, sondVolume._Volume,
                               lastVelVolume._Volume, MISSINGVEL, filt,
                               &usuccess)
    else:
        _fourdd_h.unfoldVolume(unfoldedVolume._Volume, sondVolume._Volume,
                               NULL, MISSINGVEL, filt, &usuccess)
        
    if debug:
        return usuccess, DZvolume, radialVelVolume, unfoldedVolume, sondVolume
    
    data = unfoldedVolume.get_data()
    
    return usuccess, data
