"""
Function for mathematical, signal processing and numerical routines.

"""

import numpy as np
from scipy import signal


def angular_texture_2d(image, N, interval):
    """
    Compute the angular texture of an image. Uses convolutions
    in order to speed up texture calculation by a factor of ~50
    compared to using ndimage.generic_filter.

    Parameters
    ----------
    image : 2D array of floats
        The array containing the velocities in which to calculate
        texture from.
    N : int
        This is the window size for calculating texture. The texture will be
        calculated from an N by N window centered around the gate.
    interval : float
        The absolute value of the maximum velocity. In conversion to
        radial coordinates, pi will be defined to be interval
        and -pi will be -interval. It is recommended that interval be
        set to the Nyquist velocity.

    Returns
    -------
    std_dev : float array
        Texture of the radial velocity field.

    """
    # transform distribution from original interval to [-pi, pi]
    interval_max = interval
    interval_min = -interval
    half_width = (interval_max - interval_min) / 2.
    center = interval_min + half_width

    # Calculate parameters needed for angular std. dev
    im = (np.asarray(image) - center) / (half_width) * np.pi
    x = np.cos(im)
    y = np.sin(im)

    # Calculate convolution
    kernel = np.ones((N, N))
    xs = signal.convolve2d(x, kernel, mode="same", boundary="symm")
    ys = signal.convolve2d(y, kernel, mode="same", boundary="symm")
    ns = N**2

    # Calculate norm over specified window
    xmean = xs/ns
    ymean = ys/ns
    norm = np.sqrt(xmean**2 + ymean**2)
    std_dev = np.sqrt(-2 * np.log(norm)) * (half_width) / np.pi
    return std_dev


def rolling_window(a, window):
    """ Create a rolling window object for application of functions
    eg: result=np.ma.std(array, 11), 1). """
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1], )
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def texture(radar, var):
    """ Determine a texture field using an 11pt stdev
    texarray=texture(pyradarobj, field). """
    fld = radar.fields[var]['data']
    print(fld.shape)
    tex = np.ma.zeros(fld.shape)
    for timestep in range(tex.shape[0]):
        ray = np.ma.std(rolling_window(fld[timestep, :], 11), 1)
        tex[timestep, 5:-5] = ray
        tex[timestep, 0:4] = np.ones(4) * ray[0]
        tex[timestep, -5:] = np.ones(5) * ray[-1]
    return tex


def texture_along_ray(radar, var, wind_size=7):
    """
    Compute field texture along ray using a user specified
    window size.

    Parameters
    ----------
    radar : radar object
        The radar object where the field is.
    var : str
        Name of the field which texture has to be computed.
    wind_size : int, optional
        Optional. Size of the rolling window used.

    Returns
    -------
    tex : radar field
        The texture of the specified field.

    """
    half_wind = int((wind_size-1)/2)
    fld = radar.fields[var]['data']
    tex = np.ma.zeros(fld.shape)
    for timestep in range(tex.shape[0]):
        ray = np.ma.std(rolling_window(fld[timestep, :], wind_size), 1)
        tex[timestep, half_wind:-half_wind] = ray
        tex[timestep, 0:half_wind] = np.ones(half_wind) * ray[0]
        tex[timestep, -half_wind:] = np.ones(half_wind) * ray[-1]
    return tex
