"""
pyart.util.hildebrand_sekhon
============================

Estimation of noise in Doppler spectra using the Hildebrand Sekhon method.

.. autosummary::
    :toctree: generated/

    estimate_noise_hs74

"""

import numpy as np


def estimate_noise_hs74(spectrum, navg=1):
    """
    Estimate noise parameters of a Doppler spectrum.

    Use the method of estimating the noise level in Doppler spectra outlined
    by Hildebrand and Sehkon, 1974.

    Parameters
    ----------
    spectrum : array like
        Doppler spectrum in linear units.
    navg : int, optional
        The number of spectral bins over which a moving average has been
        taken. Corresponds to the **p** variable from equation 9 of the
        article.  The default value of 1 is appropiate when no moving
        average has been applied to the spectrum.

    Returns
    -------
    mean : float-like
        Mean of points in the spectrum identified as noise.
    threshold : float-like
        Threshold separating noise from signal.  The point in the spectrum with
        this value or below should be considered as noise, above this value
        signal. It is possible that all points in the spectrum are identified
        as noise.  If a peak is required for moment calculation then the point
        with this value should be considered as signal.
    var : float-like
        Variance of the points in the spectrum identified as noise.
    nnoise : int
        Number of noise points in the spectrum.

    References
    ----------
    P. H. Hildebrand and R. S. Sekhon, Objective Determination of the Noise
    Level in Doppler Spectra. Journal of Applied Meteorology, 1974, 13,
    808-811.

    """
    sorted_spectrum = np.sort(spectrum)
    nnoise = len(spectrum)  # default to all points in the spectrum as noise
    for npts in range(1, len(sorted_spectrum)+1):
        partial = sorted_spectrum[:npts]
        mean = np.mean(partial)
        var = np.var(partial)
        if var * navg < mean**2.:
            nnoise = npts
        else:
            # partial spectrum no longer has characteristics of white noise
            break

    noise_spectrum = sorted_spectrum[:nnoise]
    mean = np.mean(noise_spectrum)
    threshold = sorted_spectrum[nnoise-1]
    var = np.var(noise_spectrum)
    return mean, threshold, var, nnoise
