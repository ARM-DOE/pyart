"""
pyart.io.fourdd
===============
Python wrapper around the Univ. Washington FourDD dealiasing code.

.. autosummary::
    :toctree: generated/

    dealias_radar
    dealias_radar_array

"""

import os
from ctypes import c_int, c_ushort, c_float, POINTER

import numpy as np

import _rsl, rsl_utils, rsl


c_float_p = POINTER(c_float)


def dealias_radar(pyradar, vad_time, py_last_radar=None, sondfile=None,
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
    radar = rsl_utils.radar_to_rsl(pyradar, {refl: 'DZ', vel: 'VR'})
    
    if py_last_radar is not None:
        last_radar = rsl_utils.radar_to_rsl(
            py_last_radar, {refl: 'DZ', vel:'VR'})

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


def dealias_radar_array(pyradar, vad_time, height, speed, direction, 
                        py_last_radar,
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
    
    
    refl='reflectivity_horizontal'
    vel='mean_doppler_velocity'
    fill_value = -9999.0
    rsl_badval = 131072

    radar = rsl_utils.radar_to_rsl(pyradar, {refl: 'DZ', vel: 'VR'})
    
    if py_last_radar is not None:
        last_radar = rsl_utils.radar_to_rsl(
            py_last_radar, {refl: 'DZ', vel:'VR'})
    else:
        last_radar = None


    # prepare arrays
    hc = np.ascontiguousarray(height, dtype='float32')
    sc = np.ascontiguousarray(speed, dtype='float32')
    dc = np.ascontiguousarray(direction, dtype='float32')

    # fixed parameters
    MISSINGVEL = 131072.0
    MAX_RADAR_VOLUMES = 41

    # create a new radar object to hold dealised results
    #dealiased = _rsl.RSL_new_radar(MAX_RADAR_VOLUMES)
    #dealiased.contents.h = radar.contents.h
    #dealiased.contents.h.nvolumes = MAX_RADAR_VOLUMES

    #if radar.contents.volumes[DZindex] is not None:
    #    DZvolume = radar.contents.v[DZindex]
    #else:
    #    print "no DZ field available"
    #    prep = 0

    #if radar.contents.volumes[_rsl.fieldTypes().VR] is None:
    #    raise ValueError('No VR field avaialable to unfold')
    #radialVelVolume = radar.contents.v[_rsl.fieldTypes().VR]

    # store last velocity volumes from last_data if present
    #if last_radar is not None:
    #    lastVelVolume = last_radar.contents.volumes[_rsl.fieldTypes().VE]
    #else:
    #    lastVelVolume = None

    # read in sounding if present
    success = c_ushort(0)
    
    DZvolume = radar.contents.v[DZindex]
    radialVelVolume = radar.contents.v[_rsl.fieldTypes().VR]
    sondVolume = _rsl.RSL_copy_volume(radialVelVolume)
    unfoldedVolume = _rsl.RSL_copy_volume(radialVelVolume)
    
    _rsl.firstGuessNoRead(sondVolume, MISSINGVEL, hc.ctypes.data_as(c_float_p),
                          sc.ctypes.data_as(c_float_p),
                          dc.ctypes.data_as(c_float_p), c_int(len(height)),
                          vad_time, success)

    if success.value == 1:
        print "dealiasing"
        if prep == 1:
            print "Prepping"
            _rsl.prepVolume(DZvolume, unfoldedVolume, MISSINGVEL)
        usuccess = c_ushort(0)
        _rsl.unfoldVolume(unfoldedVolume, sondVolume, None,
                          MISSINGVEL, filt, usuccess)
    
    dealiased_data = rsl.create_cube_array_lim(
        unfoldedVolume[0], unfoldedVolume.contents.h.nsweeps,
        unfoldedVolume.contents.sweeps[0].h.nrays)
    
    dealiased_data[np.where(np.isnan(dealiased_data))] = fill_value
    dealiased_data[np.where(dealiased_data == rsl_badval)] = fill_value
    dealiased_data = np.ma.masked_equal(dealiased_data, -9999.0)
    dealiased_data.shape = (-1, dealiased_data.shape[2])

    return dealiased_data
