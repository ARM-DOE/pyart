""" Mathematical, signal processing and numerical routines

TODO
----
Put more stuff in here

"""

from __future__ import print_function

import numpy as np


def rolling_window(a, window):
    """ create a rolling window object for application of functions
    eg: result=np.ma.std(array, 11), 1)"""
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1], )
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def texture(myradar, var):
    """Determine a texture field using an 11pt stdev
    texarray=texture(pyradarobj, field)
    """
    fld = myradar.fields[var]['data']
    print(fld.shape)
    tex = np.ma.zeros(fld.shape)
    for timestep in range(tex.shape[0]):
        ray = np.ma.std(rolling_window(fld[timestep, :], 11), 1)
        tex[timestep, 5:-5] = ray
        tex[timestep, 0:4] = np.ones(4) * ray[0]
        tex[timestep, -5:] = np.ones(5) * ray[-1]
    return tex


def texture_along_ray(myradar, var, wind_size=7):
    """
    Compute field texture along ray using a user specified
    window size.

    Parameters
    ----------
    myradar : radar object
        The radar object where the field is
    var : str
        Name of the field which texture has to be computed
    wind_size : int
        Optional. Size of the rolling window used

    Returns
    -------
    tex : radar field
        the texture of the specified field

    """
    half_wind = int((wind_size-1)/2)
    fld = myradar.fields[var]['data']
    tex = np.ma.zeros(fld.shape)
    for timestep in range(tex.shape[0]):
        ray = np.ma.std(rolling_window(fld[timestep, :], wind_size), 1)
        tex[timestep, half_wind:-half_wind] = ray
        tex[timestep, 0:half_wind] = np.ones(half_wind) * ray[0]
        tex[timestep, -half_wind:] = np.ones(half_wind) * ray[-1]
    return tex
