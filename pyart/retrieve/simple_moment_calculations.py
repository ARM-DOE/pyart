"""
Simple moment calculations.

"""

import numpy as np

from scipy import ndimage
from ..config import get_metadata, get_field_name
from ..core.transforms import antenna_to_cartesian
from ..util import angular_texture_2d


def calculate_snr_from_reflectivity(
        radar, refl_field=None, snr_field=None, toa=25000.):
    """
    Calculate the signal to noise ratio, in dB, from the reflectivity field.

    Parameters
    ----------
    radar : Radar
        Radar object from which to retrieve reflectivity field.
    refl_field : str, optional
        Name of field in radar which contains the reflectivity.
        None will use the default field name in the Py-ART configuration file.
    snr_field : str, optional
        Name to use for snr metadata. None will use the default field name
        in the Py-ART configuration file.
    toa : float, optional
        Height above which to take noise floor measurements, in meters.

    Returns
    -------
    snr : field dictionary
        Field dictionary containing the signal to noise ratio.

    """
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if snr_field is None:
        snr_field = get_field_name('signal_to_noise_ratio')

    range_grid = np.meshgrid(radar.range['data'],
                             np.ma.ones(radar.time['data'].shape))[0] + 1.0

    # remove range scale.. This is basically the radar constant scaled dBm
    pseudo_power = (radar.fields[refl_field]['data'] -
                    20.0*np.log10(range_grid / 1000.0))

    # Noise floor estimate
    # 25km.. should be no scatterers, not even planes, this high
    # we could get undone by AP though.. also sun
    rg, azg = np.meshgrid(radar.range['data'], radar.azimuth['data'])
    rg, eleg = np.meshgrid(radar.range['data'], radar.elevation['data'])
    _, _, z = antenna_to_cartesian(rg / 1000.0, azg, eleg)  # XXX: need to fix

    points_above = np.where(z > toa)
    noise_floor_estimate = pseudo_power[points_above].mean()

    snr_dict = get_metadata(snr_field)
    snr_dict['data'] = pseudo_power - noise_floor_estimate
    return snr_dict


def compute_noisedBZ(nrays, noisedBZ_val, _range, ref_dist,
                     noise_field=None):
    """
    Computes noise in dBZ from reference noise value.

    Parameters
    ----------
    nrays : int
        Number of rays in the reflectivity field.
    noisedBZ_val : float
        Estimated noise value in dBZ at reference distance.
    _range : np array of floats
        Range vector in m.
    ref_dist : float
        Reference distance in Km.
    noise_field : str, optional
        Name of the noise field.

    Returns
    -------
    noisedBZ : dict
        The noise field.

    """
    # parse the field parameters
    if noise_field is None:
        noise_field = get_field_name('noisedBZ_hh')

    noisedBZ_vec = noisedBZ_val+20.*np.ma.log10(1e-3*_range/ref_dist)

    noisedBZ = get_metadata(noise_field)
    noisedBZ['data'] = np.tile(noisedBZ_vec, (nrays, 1))

    return noisedBZ


def compute_snr(radar, refl_field=None, noise_field=None, snr_field=None):
    """
    Computes SNR from a reflectivity field and the noise in dBZ.

    Parameters
    ----------
    radar : Radar
        Radar object
    refl_field : str, optional
        Name of the reflectivity field to use.
    noise_field : str, optional
        Name of the noise field to use.
    snr_field : str, optional
        Name of the SNR field.

    Returns
    -------
    snr : dict
        The SNR field.

    """
    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if noise_field is None:
        noise_field = get_field_name('noisedBZ_hh')
    if snr_field is None:
        snr_field = get_field_name('signal_to_noise_ratio')

    # extract fields from radar
    if refl_field in radar.fields:
        refl = radar.fields[refl_field]['data']
    else:
        raise KeyError('Field not available: ' + refl_field)
    if noise_field in radar.fields:
        noisedBZ = radar.fields[noise_field]['data']
    else:
        raise KeyError('Field not available: ' + noise_field)

    snr_data = refl-noisedBZ

    snr = get_metadata(snr_field)
    snr['data'] = snr_data

    return snr


def compute_l(radar, rhohv_field=None, l_field=None):
    """
    Computes Rhohv in logarithmic scale according to L=-log10(1-RhoHV).

    Parameters
    ----------
    radar : Radar
        Radar object.
    rhohv_field : str, optional
        Name of the RhoHV field to use.
    l_field : str, optional
        Name of the L field.

    Returns
    -------
    l : dict
        L field.

    """
    # parse the field parameters
    if rhohv_field is None:
        rhohv_field = get_field_name('cross_correlation_ratio')
    if l_field is None:
        l_field = get_field_name('logarithmic_cross_correlation_ratio')

    # extract rhohv field from radar
    if rhohv_field in radar.fields:
        rhohv = radar.fields[rhohv_field]['data']
    else:
        raise KeyError('Field not available: ' + rhohv_field)

    rhohv[rhohv >= 1.] = 0.9999
    l_data = -np.ma.log10(1.-rhohv)

    l = get_metadata(l_field)
    l['data'] = l_data

    return l


def compute_cdr(radar, rhohv_field=None, zdr_field=None, cdr_field=None):
    """
    Computes the Circular Depolarization Ratio.

    Parameters
    ----------
    radar : Radar
        Radar object.
    rhohv_field : str, optional
        Name of the RhoHV field.
    zdr_field : str, optional
        Name of the ZDR field.
    cdr_field : str, optional
        Name of the CDR field.

    Returns
    -------
    cdr : dict
        CDR field.

    """
    # parse the field parameters
    if rhohv_field is None:
        rhohv_field = get_field_name('cross_correlation_ratio')
    if zdr_field is None:
        zdr_field = get_field_name('differential_reflectivity')
    if cdr_field is None:
        cdr_field = get_field_name('circular_depolarization_ratio')

    # extract fields from radar
    if rhohv_field in radar.fields:
        rhohv = radar.fields[rhohv_field]['data']
    else:
        raise KeyError('Field not available: ' + rhohv_field)
    if zdr_field in radar.fields:
        zdrdB = radar.fields[zdr_field]['data']
    else:
        raise KeyError('Field not available: ' + zdr_field)

    zdr = np.ma.power(10., 0.1*zdrdB)

    cdr_data = (
        10.*np.ma.log10(
            (1.+1./zdr-2.*rhohv*np.ma.sqrt(1./zdr)) /
            (1.+1./zdr+2.*rhohv*np.ma.sqrt(1./zdr))))

    cdr = get_metadata(cdr_field)
    cdr['data'] = cdr_data

    return cdr


def calculate_velocity_texture(radar, vel_field=None, wind_size=4, nyq=None,
                               check_nyq_uniform=True):
    """
    Derive the texture of the velocity field.

    Parameters
    ----------
    radar: Radar
        Radar object from which velocity texture field will be made.
    vel_field : str, optional
        Name of the velocity field. A value of None will force Py-ART to
        automatically determine the name of the velocity field.
    wind_size : int, optional
        The size of the window to calculate texture from. The window is
        defined to be a square of size wind_size by wind_size.
    nyq : float, optional
        The nyquist velocity of the radar. A value of None will force Py-ART
        to try and determine this automatically.
    check_nyquist_uniform : bool, optional
        True to check if the Nyquist velocities are uniform for all rays
        within a sweep, False will skip this check. This parameter is ignored
        when the nyq parameter is not None.

    Returns
    -------
    vel_dict: dict
        A dictionary containing the field entries for the radial velocity
        texture.

    """
    # Parse names of velocity field
    if vel_field is None:
        vel_field = get_field_name('velocity')

    # Allocate memory for texture field
    vel_texture = np.zeros(radar.fields[vel_field]['data'].shape)

    # If an array of nyquist velocities is derived, use different
    # nyquist velocites for each sweep in texture calculation according to
    # the nyquist velocity in each sweep.

    if nyq is None:
        # Find nyquist velocity if not specified
        nyq = [radar.get_nyquist_vel(i, check_nyq_uniform) for i in
               range(radar.nsweeps)]
        for i in range(0, radar.nsweeps):
            start_ray, end_ray = radar.get_start_end(i)
            inds = range(start_ray, end_ray)
            vel_sweep = radar.fields[vel_field]['data'][inds]
            vel_texture[inds] = angular_texture_2d(
                vel_sweep, wind_size, nyq[i])
    else:
        vel_texture = angular_texture_2d(
            radar.fields[vel_field]['data'], wind_size, nyq)
    vel_texture_field = get_metadata('velocity')
    vel_texture_field['long_name'] = 'Doppler velocity texture'
    vel_texture_field['standard_name'] = (
        'texture_of_radial_velocity' + '_of_scatters_away_from_instrument')
    vel_texture_field['data'] = ndimage.filters.median_filter(vel_texture,
                                                              size=(wind_size,
                                                                    wind_size))
    return vel_texture_field
