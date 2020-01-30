"""
Estimation of noise in Doppler spectra using the Hildebrand Sekhon method.

"""

import numpy as np

def estimate_noise_hs74(spectrum, navg=1, nnoise_min=1):
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
        article. The default value of 1 is appropriate when no moving
        average has been applied to the spectrum.
    nnoise_min : int, optional
        Minimum number of noise samples to consider the estimation valid.

    Returns
    -------
    mean : float-like
        Mean of points in the spectrum identified as noise.
    threshold : float-like
        Threshold separating noise from signal. The point in the spectrum with
        this value or below should be considered as noise, above this value
        signal. It is possible that all points in the spectrum are identified
        as noise. If a peak is required for moment calculation then the point
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

    rtest = 1+1/navg
    sum1 = 0.
    sum2 = 0.
    for i, pwr in enumerate(sorted_spectrum):
        npts = i+1
        sum1 += pwr
        sum2 += pwr*pwr

        if npts < nnoise_min:
            continue

        if npts*sum2 < sum1*sum1*rtest:
            nnoise = npts
        else:
            # partial spectrum no longer has characteristics of white noise.
            sum1 -= pwr
            sum2 -= pwr*pwr
            break

    mean = sum1/nnoise
    var = sum2/nnoise-mean*mean
    threshold = sorted_spectrum[nnoise-1]
    return mean, threshold, var, nnoise
