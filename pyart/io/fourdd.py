""" Python wrapper around the Univ. Washington FourDD dealiasing code. """

import os
from ctypes import c_int, c_ushort, c_float, POINTER

import numpy as np

import _rsl


c_float_p = POINTER(c_float)


def dealias_radar(radar, vad_time, last_radar=None, sondfile=None,
                  prep=1, DZindex=_rsl.fieldTypes().DZ, filt=1):
    """
    Dealias a RSL Radar object.

    Parameters
    ----------
    radar : pointer to a Radar object (LP_Radar)
        RSL Radar object to dealias.
    vad_time : int
        Sounding time, format must be YYDDDHHMM.
    last_radar : pointer to Radar object (LP_Radar) or None
        RSL Radar object containg last measurements. Either last_radar or
        sondfile must be supplied.
    sondfile : str or None.
        Filename of sounding file.  Either sondfile or last_radar must be
        supplied.
    prep : int, optional
        Thresholding flag, 1 = yes, 0 = no.
    DZindex : int, optional
        Index of DZ measurements in RSL Radar object.
    filt : int. optional
        Flag controlling Bergen and Albers filter, 1 = yes, 0 = no.

    Returns
    -------
    unfoldedVolume : pointer to Radar object (LP_Radar)
        Dealised RSL Radar object.

    """
    # fixed parameters
    MISSINGVEL = 131072.0
    MAX_RADAR_VOLUMES = 41

    # create a new radar object to hold the dealiased results
    dealiased = _rsl.RSL_new_radar(MAX_RADAR_VOLUMES)
    dealiased.contents.h = radar.contents.h
    dealiased.contents.h.nvolumes = MAX_RADAR_VOLUMES

    if radar.contents.volumes[DZindex] is not None:
        DZvolume = radar.contents.v[DZindex]
    else:
        print "no DZ field available"
        prep = 0

    if radar.contents.volumes[_rsl.fieldTypes().VR] is None:
        raise ValueError('No VR field avaialable to unfold')
    radialVelVolume = radar.contents.v[_rsl.fieldTypes().VR]

    # store last Velocity Volumes from last_radar if present
    if last_radar is not None:
        lastVelVolume = last_radar.contents.volumes[_rsl.fieldTypes().VE]
    else:
        lastVelVolume = None

    # read in sounding if present
    success = c_ushort(0)
    if os.path.isfile(sondfile):
        sondVolume = _rsl.RSL_copy_volume(radialVelVolume)
        _rsl.firstGuess(sondVolume, MISSINGVEL, sondfile, vad_time, success)
    else:
        sondVolume = None
        print "file does not exist"

    if success.value == 1 or lastVelVolume is not None:
        print "dealiasing"
        unfoldedVolume = _rsl.RSL_copy_volume(radialVelVolume)
        if prep == 1:
            _rsl.prepVolume(DZvolume, unfoldedVolume, MISSINGVEL)
        usuccess = c_ushort(0)
        _rsl.unfoldVolume(unfoldedVolume, sondVolume, lastVelVolume,
                          MISSINGVEL, filt, usuccess)
    return unfoldedVolume


def dealias_radar_array(radar, vad_time, height, speed, direction, last_radar,
                        prep=1, DZindex=_rsl.fieldTypes().DZ, filt=1):
    """
    Dealias a RSL Radar object using sond arrays.

    Parameters
    ----------
    radar : pointer to a Radar object (LP_Radar)
        RSL Radar object to dealias.
    vad_time : int
        Sounding time, format must be YYDDDHHMM.
    height : array_like
        Height of wind profile.
    speed : array_like
        Speed of wind profile.
    direction : array_like
        Direction of wind profile.
    last_radar : pointer to Radar object (LP_Radar) or None, optional.
        RSL Radar object containg last measurements.
    prep : int, optional
        Flag controlling thresholding, 1 = yes, 0 = no.
    DZindex : int, optional
        Index of DZ measurements in RSL Radar object.
    filt : int. optional
        Flag controlling Bergen and Albers filter, 1 = yes, 0 = no.


    Returns
    -------
    unfoldedVolume : pointer to Radar object (LP_Radar)
        Dealised RSL Radar object.
    sondVolume : XXX
        XXX

    """

    # prepare arrays
    hc = np.ascontiguousarray(height, dtype='float32')
    sc = np.ascontiguousarray(speed, dtype='float32')
    dc = np.ascontiguousarray(direction, dtype='float32')

    # fixed parameters
    MISSINGVEL = 131072.0
    MAX_RADAR_VOLUMES = 41

    # create a new radar object to hold dealised results
    dealiased = _rsl.RSL_new_radar(MAX_RADAR_VOLUMES)
    dealiased.contents.h = radar.contents.h
    dealiased.contents.h.nvolumes = MAX_RADAR_VOLUMES

    if radar.contents.volumes[DZindex] is not None:
        DZvolume = radar.contents.v[DZindex]
    else:
        print "no DZ field available"
        prep = 0

    if radar.contents.volumes[_rsl.fieldTypes().VR] is None:
        raise ValueError('No VR field avaialable to unfold')
    radialVelVolume = radar.contents.v[_rsl.fieldTypes().VR]

    # store last velocity volumes from last_data if present
    if last_radar is not None:
        lastVelVolume = last_radar.contents.volumes[_rsl.fieldTypes().VE]
    else:
        lastVelVolume = None

    # read in sounding if present
    success = c_ushort(0)
    sondVolume = _rsl.RSL_copy_volume(radialVelVolume)
    _rsl.firstGuessNoRead(sondVolume, MISSINGVEL, hc.ctypes.data_as(c_float_p),
                          sc.ctypes.data_as(c_float_p),
                          dc.ctypes.data_as(c_float_p), c_int(len(height)),
                          vad_time, success)

    if success.value == 1 or lastVelVolume is not None:
        print "dealiasing"
        unfoldedVolume = _rsl.RSL_copy_volume(radialVelVolume)
        if prep == 1:
            _rsl.prepVolume(DZvolume, unfoldedVolume, MISSINGVEL)
        usuccess = c_ushort(0)
        _rsl.unfoldVolume(unfoldedVolume, sondVolume, lastVelVolume,
                          MISSINGVEL, filt, usuccess)
    return unfoldedVolume, sondVolume
