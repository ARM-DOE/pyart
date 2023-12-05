"""
Create CFAD from a radar or grid field

"""

import numpy as np


def create_cfad(
    field_data,
    altitude_data,
    field_bins,
    altitude_bins,
    min_frac_thres=0.1,
):
    """
    This function returns a Contoured Frequency by Altitude Diagram (CFAD; Yuter et al. 1995), a 2-dimensional
    histogram that is normalized by the number of points at each altitude. Altitude bins are masked where the counts
    are less than a minimum fraction of the largest number of counts for any altitude row.

    field_data : array
        Array of radar data to use for CFAD calculation.
    altitude_data : array
        Array of corresponding altitude data to use for CFAD calculation.
        Note: must be the same shape as `field_data`
    field_bins : list
        List of bin edges for field values to use for CFAD creation.
    altitude_bins : list
        List of bin edges for height values to use for CFAD creation.
    min_frac_thres : float, optional
        Fraction of values to remove in CFAD normalization (default 0.1). If an altitude row has a total count that
        is less than min_frac_thres of the largest number of total counts for any altitude row, the bins in that
        altitude row are masked.
    Returns
    -------
    freq_norm : array
        Array of normalized frequency.
    height_edges : array
        Array of bin edges for height data.
    field_edges : array of x coordinates
        Array of bin edges for field data.
    References
    ----------
    Yuter, S. E., and R. A. Houze, 1995: Three-Dimensional Kinematic and
    Microphysical Evolution of Florida Cumulonimbus. Part II: Frequency Distributions
    of Vertical Velocity, Reflectivity, and Differential Reflectivity. Mon. Wea. Rev.
    123, 1941-1963. https://doi.org/10.1175/1520-0493(1995)123%3C1941:TDKAME%3E2.0.CO;2


    """
    # get raw bin counts
    freq, height_edges, field_edges = np.histogram2d(altitude_data.compressed(), field_data.compressed(),
                                                     bins=[altitude_bins, field_bins])

    # sum counts over y axis (height)
    freq_sum = np.sum(freq, axis=1)
    # get threshold for normalizing
    point_thres = min_frac_thres * np.max(freq_sum)
    # repeat to create array same size as freq
    freq_sum_rep = np.repeat(freq_sum[..., np.newaxis], freq.shape[1], axis=1)
    # normalize
    freq_norm = freq / freq_sum_rep
    # mask data where there is not enough points
    freq_norm = np.ma.masked_where(freq_sum_rep < point_thres, freq_norm)

    return freq_norm, height_edges, field_edges
