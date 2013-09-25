"""
pyart.correct._fourdd_interface
===============================

Cython wrapper around the Univ. of Washington FourDD algorithm.

.. autosummary::
    :toctree: generated/

    fourdd_dealias

"""

# TODO
# better FourDD dealias

cimport _fourdd_h
cimport numpy as np
from pyart.io._rsl_interface cimport _RslVolume

from ..io import _rsl_interface


cpdef fourdd_dealias(_RslVolume DZvolume,
                     _RslVolume radialVelVolume,
                     np.ndarray[np.float32_t, ndim=1] hc,
                     np.ndarray[np.float32_t, ndim=1] sc,
                     np.ndarray[np.float32_t, ndim=1] dc,
                     vad_time, prep, filt, debug):
    """
    fourdd_dealias(DZvolume, radialVelVolume, hc, sc, dc, vad_time, prep,
                   filt)

    Dealias using the FourDD algorithm.

    Parameters
    ----------
    DZvolume : _RslVolume
        Reflectivity to use when thresholding is selected.
    radialVelVolume : _RslVolume
        Radial velocities which will be dealiased.
    hc : array
        Sounding heights in meters.  Must be a contiguous one-dimensional
        float32 array.
    sc : array
        Sounding wind speed in m/s.  Must be a contiguous one-dimensional
        float32 array.
    dc : array
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
    data : array
        Array of unfolded velocities.

    References
    ----------
    C. N. James and R. A Houze Jr, A Real-Time Four-Dimensional Doppler
    Dealising Scheme, Journal of Atmospheric and Oceanic Technology, 2001, 18,
    1674.

    """
    # TODO version which does not need DZvolume (prep is 0)
    # TODO version which uses
    # See FourDD.c for addition details
    cdef _RslVolume unfoldedVolume
    cdef _RslVolume sondVolume
    cdef float MISSINGVEL = 131072.0
    cdef unsigned short success = 0
    cdef unsigned short usuccess = 0

    unfoldedVolume = _rsl_interface.copy_volume(radialVelVolume)
    sondVolume = _rsl_interface.copy_volume(radialVelVolume)

    # May not always be needed...
    _fourdd_h.firstGuessNoRead(
        sondVolume._Volume, MISSINGVEL, <float *> hc.data, <float *> sc.data,
        <float *> dc.data, <int> len(hc), vad_time, &success)

    if success != 1:
        raise ValueError

    # dealias
    if prep == 1:
        _fourdd_h.prepVolume(DZvolume._Volume, unfoldedVolume._Volume,
                             MISSINGVEL)
    _fourdd_h.unfoldVolume(unfoldedVolume._Volume, sondVolume._Volume, NULL,
                           MISSINGVEL, filt, &usuccess)
    data = unfoldedVolume.get_data()
    if debug:
        return usuccess, DZvolume, radialVelVolume, unfoldedVolume, sondVolume
    return usuccess, data
