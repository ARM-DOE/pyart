"""
Corrects polarimetric variables for noise

"""

import warnings

import dask
import numpy as np
import pint
from scipy import signal

from ..config import get_field_name, get_metadata
from ..core.radar import Radar


def correct_noise_rhohv(
    radar,
    urhohv_field=None,
    snr_field=None,
    zdr_field=None,
    nh_field=None,
    nv_field=None,
    rhohv_field=None,
):
    """
    Corrects RhoHV for noise according to eq. 6 in Gourley et al. 2006.
    This correction should only be performed if noise has not been subtracted
    from the signal during the moments computation.

    Parameters
    ----------
    radar : Radar
        Radar object.
    urhohv_field : str, optional
        Name of the RhoHV uncorrected for noise field.
    snr_field, zdr_field, nh_field, nv_field : str, optional
        Names of the SNR, ZDR, horizontal channel noise in dBZ and vertical
        channel noise in dBZ used to correct RhoHV.
    rhohv_field : str, optional
        Name of the rhohv field to output.

    Returns
    -------
    rhohv : dict
        Noise corrected RhoHV field.

    References
    ----------
    Gourley et al. Data Quality of the Meteo-France C-Band Polarimetric
    Radar, JAOT, 23, 1340-1356

    """
    # parse the field parameters
    if urhohv_field is None:
        urhohv_field = get_field_name("uncorrected_cross_correlation_ratio")
    if snr_field is None:
        snr_field = get_field_name("signal_to_noise_ratio")
    if zdr_field is None:
        zdr_field = get_field_name("differential_reflectivity")
    if nh_field is None:
        nh_field = get_field_name("noisedBZ_hh")
    if nv_field is None:
        nv_field = get_field_name("noisedBZ_vv")
    if rhohv_field is None:
        rhohv_field = get_field_name("cross_correlation_ratio")

    # extract fields from radar
    if urhohv_field in radar.fields:
        urhohv = radar.fields[urhohv_field]["data"]
    else:
        raise KeyError("Field not available: " + urhohv_field)
    if snr_field in radar.fields:
        snrdB_h = radar.fields[snr_field]["data"]
    else:
        raise KeyError("Field not available: " + snr_field)
    if zdr_field in radar.fields:
        zdrdB = radar.fields[zdr_field]["data"]
    else:
        raise KeyError("Field not available: " + zdr_field)
    if nh_field in radar.fields:
        nh = radar.fields[nh_field]["data"]
    else:
        raise KeyError("Field not available: " + nh_field)
    if nv_field in radar.fields:
        nv = radar.fields[nv_field]["data"]
    else:
        raise KeyError("Field not available: " + nv_field)

    snr_h = np.ma.power(10.0, 0.1 * snrdB_h)
    zdr = np.ma.power(10.0, 0.1 * zdrdB)
    alpha = np.ma.power(10.0, 0.1 * (nh - nv))

    rhohv_data = urhohv * np.ma.sqrt(
        (1.0 + 1.0 / snr_h) * (1.0 + zdr / (alpha * snr_h))
    )
    rhohv_data[rhohv_data > 1.0] = 1.0

    rhohv = get_metadata(rhohv_field)
    rhohv["data"] = rhohv_data

    return rhohv


def correct_bias(radar, bias=0.0, field_name=None):
    """
    Corrects a radar data bias. If field name is none the correction is
    applied to horizontal reflectivity by default.

    Parameters
    ----------
    radar : Radar
        Radar object.
    bias : float, optional
        The bias magnitude.
    field_name: str, optional
        Names of the field to be corrected.

    Returns
    -------
    corrected_field : dict
        The corrected field

    """
    # parse the field parameters
    if field_name is None:
        field_name = get_field_name("reflectivity")

    # extract fields from radar
    if field_name in radar.fields:
        field_data = radar.fields[field_name]["data"]
    else:
        raise KeyError("Field not available: " + field_name)

    corr_field_data = field_data - bias

    if field_name.startswith("corrected_"):
        corr_field_name = field_name
    else:
        corr_field_name = "corrected_" + field_name

    corr_field = get_metadata(corr_field_name)
    corr_field["data"] = corr_field_data

    return corr_field


def calc_zdr_offset(radar, gatefilter=None, height_range=None, zdr_var=None):
    """
    Function for calculating the ZDR bias from a VPT scan.

    Parameters
    ----------
    radar : PyART radar object
        Radar object with radar data
    gatefilter: PyART GateFilter
        Gatefilter for filtering out data for calculating ZDR bias
    height_range: tuple
        The minimum and maximum heights in meters for the scan.
    zdr_var: string or None
        The name of the ZDR variable. Set to None to have PyART try to determine this automatically.

    Returns
    -------
    profiles : dict
        The mean vertical profiles of each radar moment are extracted along with the ZDR bias.

    """
    if height_range is None:
        height_range = (0, 100000.0)

    if zdr_var is None:
        zdr_var = get_field_name("differential_reflectivity")

    height_mask = np.logical_and(
        radar.range["data"] >= height_range[0], radar.range["data"] <= height_range[1]
    )

    mask = gatefilter.gate_included
    Zdr = radar.fields[zdr_var]["data"]
    if isinstance(Zdr, np.ma.MaskedArray):
        Zdr = Zdr.filled(np.nan)
    Zdr = np.where(mask, Zdr, np.nan)
    bias = np.nanmean(Zdr[:, height_mask])
    height_mask_tiled = np.tile(height_mask, (radar.nrays, 1))
    Zdr = np.where(height_mask_tiled, Zdr, np.nan)
    results = {
        "bias": bias,
        "profile_zdr": np.nanmean(Zdr[:, :], axis=0),
        "range": radar.range["data"],
    }
    for k in radar.fields.keys():
        if k != "range":
            field_data = radar.fields[k]["data"].astype(float)
            if isinstance(field_data, np.ma.MaskedArray):
                field_data = field_data.filled(np.nan)

            field_data = np.where(mask, field_data, np.nan)
            field_data = np.where(height_mask_tiled, field_data, np.nan)
            results["profile_" + k] = np.nanmean(field_data[:, :], axis=0)
    return results


def calc_cloud_mask(
    radar,
    field,
    height=None,
    noise_threshold=-45.0,
    threshold_offset=5.0,
    counts_threshold=12,
):
    """
    Primary function for calculating the cloud mask.

    Parameters
    ----------
    radar : Radar
        Py-ART Radar object.
    field : string
        Reflectivity field name to calculate.
    height : string
        Height name to use for calculations.
    noise_threshold : float
        Threshold value used for noise detection. Greater than this value.
    threshold_offset : float
        Threshold offset value used for noise detection
    counts_threshold : int
        Threshold of counts used to determine mask. Greater than or equal to this value.

    Returns
    -------
    radar : Radar
        Returns an updated Radar object with cloud mask fields.

    References
    ----------
    Kollias, P., I. Jo, P. Borque, A. Tatarevic, K. Lamer, N. Bharadwaj, K. Widener,
    K. Johnson, and E.E. Clothiaux, 2014: Scanning ARM Cloud Radars. Part II: Data
    Quality Control and Processing. J. Atmos. Oceanic Technol., 31, 583–598,
    https://doi.org/10.1175/JTECH-D-13-00045.1

    """

    if not isinstance(radar, Radar):
        raise ValueError("Please use a valid Py-ART Radar object")

    if not isinstance(field, str):
        raise ValueError("Please specify a valid field name.")

    noise = calc_noise_floor(radar, field, height=height)

    noise_thresh = (
        np.nanmin(
            np.vstack(
                [
                    noise,
                    np.full(np.shape(radar.fields[field]["data"])[0], noise_threshold),
                ]
            ),
            axis=0,
        )
        + threshold_offset
    )

    data = range_correction(radar, field, height=height)

    task = []
    for i in range(np.shape(data)[0]):
        task.append(dask.delayed(_first_mask)(data[i, :], noise_thresh[i]))

    result = dask.compute(task)
    mask1 = np.array(result[0])

    counts = signal.convolve2d(mask1, np.ones((4, 4), dtype=int), mode="same")
    mask2 = np.zeros_like(data, dtype=np.int16)
    mask2[counts >= counts_threshold] = 1

    cloud_mask_1 = {}
    cloud_mask_1["long_name"] = "Cloud mask 1 (linear profile)"
    cloud_mask_1["units"] = "1"
    cloud_mask_1["comment"] = (
        "The mask is calculated with a " "linear mask along each time profile."
    )
    cloud_mask_1["flag_values"] = [0, 1]
    cloud_mask_1["flag_meanings"] = ["no_cloud", "cloud"]
    cloud_mask_1["data"] = mask1

    cloud_mask_2 = {}
    cloud_mask_2["long_name"] = "Cloud mask 2 (2D box)"
    cloud_mask_2["units"] = "1"
    cloud_mask_2["comment"] = "The mask uses a 2D box to " "filter out noise."
    cloud_mask_2["flag_values"] = [0, 1]
    cloud_mask_2["flag_meanings"] = ["no_cloud", "cloud"]
    cloud_mask_2["data"] = mask2

    radar.add_field("cloud_mask_1", cloud_mask_1, replace_existing=True)
    radar.add_field("cloud_mask_2", cloud_mask_2, replace_existing=True)

    return radar


def calc_noise_floor(radar, field, height):
    """
    Calculation for getting the noise floor

    Parameters
    ----------
    radar : Radar
        Py-ART Radar object containing data.
    field : string
        Reflectivity field name to correct.
    height : string
        Height name to use in correction.

    Returns
    -------
    noise : array
        Returns the noise floor value for each time sample.

    References
    ----------
    Kollias, P., I. Jo, P. Borque, A. Tatarevic, K. Lamer, N. Bharadwaj, K. Widener,
    K. Johnson, and E.E. Clothiaux, 2014: Scanning ARM Cloud Radars. Part II: Data
    Quality Control and Processing. J. Atmos. Oceanic Technol., 31, 583–598,
    https://doi.org/10.1175/JTECH-D-13-00045.1

    """

    if not isinstance(radar, Radar):
        raise ValueError("Please use a valid Py-ART Radar object.")

    # Range correct data and return the array from the Radar object
    data = range_correction(radar, field, height=height)

    # Pass each timestep into task list to calculate cloud threshhold
    # with a delayed dask process
    task = [dask.delayed(cloud_threshold)(row) for row in data]

    # Perform dask computation
    noise = dask.compute(*task)

    # Convert returned dask tuple into numpy array
    noise = np.array(noise, dtype=float)

    return noise


def cloud_threshold(data, n_avg=1.0, nffts=None):
    """
    Calculates the noise floor from a cloud threshold.

    Parameters
    ----------
    data : array
        Numpy array
    n_avg : float
        Number of points to average over
    nffts : int
        Number of heights to iterate over. If None will use the size of data.

    Returns
    -------
    result : numpy scalar float
        Returns the noise floor value for each time sample.

    References
    ----------
    Kollias, P., I. Jo, P. Borque, A. Tatarevic, K. Lamer, N. Bharadwaj, K. Widener,
    K. Johnson, and E.E. Clothiaux, 2014: Scanning ARM Cloud Radars. Part II: Data
    Quality Control and Processing. J. Atmos. Oceanic Technol., 31, 583–598,
    https://doi.org/10.1175/JTECH-D-13-00045.1

    """

    if nffts is None:
        nffts = data.size

    data = 10.0 ** (data / 10.0)
    data = np.sort(data)

    nthld = 10.0**-10.0
    dsum = 0.0
    sumSq = 0.0
    n = 0.0
    numNs = []
    sqrt_n_avg = np.sqrt(n_avg)
    for i in range(nffts):
        if data[i] > nthld:
            dsum += data[i]
            sumSq += data[i] ** 2.0
            n += 1.0
            a3 = dsum * dsum
            a1 = sqrt_n_avg * (n * sumSq - a3)
            if n > nffts / 4.0:
                if a1 <= a3:
                    sumNs = dsum
                    numNs = [n]
            else:
                sumNs = dsum
                numNs = [n]

    if len(numNs) > 0:
        n_mean = sumNs / numNs[0]
    else:
        n_mean = np.nan

    if n_mean == 0.0:
        result = np.nan
    else:
        result = 10.0 * np.log10(n_mean)

    return result


def range_correction(radar, field, height):
    """
    Corrects reflectivity for range to help get the
    correct noise floor values

    Parameters
    ----------
    radar : Radar
        Py-ART Radar object containing data.
    field : string
        Reflectivity field name to correct.
    height : string
        Height name to use in correction.

    Returns
    -------
    data : array
        Returns a range corrected array matching reflectivity field.

    """

    if not isinstance(radar, Radar):
        raise ValueError("Please use a valid Py-ART Radar object.")

    try:
        height_units = getattr(radar, height)["units"]
    except KeyError:
        warnings.warn(
            f"Height '{height} does not have units attribute. "
            "Assuming units are meters."
        )
        height_units = "m"

    height = getattr(radar, height)["data"]
    desired_unit = "m"
    if height_units is not desired_unit:
        if isinstance(height, np.ma.MaskedArray):
            height = height.filled(np.nan)
        unit_registry = pint.UnitRegistry()
        height = height * unit_registry.parse_expression(height_units)
        height = height.to(desired_unit)
        height = height.magnitude

    data = radar.fields[field]["data"]

    if isinstance(data, np.ma.MaskedArray) and not data.mask:
        data = data.data
    elif isinstance(data, np.ma.MaskedArray) and data.mask:
        data = data.filled(np.nan)

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", category=RuntimeWarning, message=".*divide by zero encountered.*"
        )
        data = data - 20.0 * np.log10(height / 1000.0)

    return data


def _first_mask(data, noise_threshold):
    mask = np.zeros_like(data, dtype=np.int16)
    mask[data > noise_threshold] = 1
    return mask
