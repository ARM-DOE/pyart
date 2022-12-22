"""
Functions for echo classification.

"""

from warnings import warn

import numpy as np

from ..config import get_field_name, get_fillvalue, get_metadata
from ._echo_class import _revised_conv_strat, steiner_class_buff


def steiner_conv_strat(
    grid,
    dx=None,
    dy=None,
    intense=42.0,
    work_level=3000.0,
    peak_relation="default",
    area_relation="medium",
    bkg_rad=11000.0,
    use_intense=True,
    fill_value=None,
    refl_field=None,
):
    """
    Partition reflectivity into convective-stratiform using the Steiner et
    al. (1995) algorithm.

    Parameters
    ----------
    grid : Grid
        Grid containing reflectivity field to partition.
    dx, dy : float, optional
        The x- and y-dimension resolutions in meters, respectively. If None
        the resolution is determined from the first two axes values.
    intense : float, optional
        The intensity value in dBZ. Grid points with a reflectivity
        value greater or equal to the intensity are automatically
        flagged as convective. See reference for more information.
    work_level : float, optional
        The working level (separation altitude) in meters. This is the height
        at which the partitioning will be done, and should minimize bright band
        contamination. See reference for more information.
    peak_relation : 'default' or 'sgp', optional
        The peakedness relation. See reference for more information.
    area_relation : 'small', 'medium', 'large', or 'sgp', optional
        The convective area relation. See reference for more information.
    bkg_rad : float, optional
        The background radius in meters. See reference for more information.
    use_intense : bool, optional
        True to use the intensity criteria.
    fill_value : float, optional
         Missing value used to signify bad data points. A value of None
         will use the default fill value as defined in the Py-ART
         configuration file.
    refl_field : str, optional
         Field in grid to use as the reflectivity during partitioning. None
         will use the default reflectivity field name from the Py-ART
         configuration file.

    Returns
    -------
    eclass : dict
        Steiner convective-stratiform classification dictionary.

    References
    ----------
    Steiner, M. R., R. A. Houze Jr., and S. E. Yuter, 1995: Climatological
    Characterization of Three-Dimensional Storm Structure from Operational
    Radar and Rain Gauge Data. J. Appl. Meteor., 34, 1978-2007.

    """
    # Get fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field parameters
    if refl_field is None:
        refl_field = get_field_name("reflectivity")

    # parse dx and dy
    if dx is None:
        dx = grid.x["data"][1] - grid.x["data"][0]
    if dy is None:
        dy = grid.y["data"][1] - grid.y["data"][0]

    # Get coordinates
    x = grid.x["data"]
    y = grid.y["data"]
    z = grid.z["data"]

    # Get reflectivity data
    ze = np.ma.copy(grid.fields[refl_field]["data"])
    ze = ze.filled(np.NaN)

    eclass = steiner_class_buff(
        ze,
        x,
        y,
        z,
        dx=dx,
        dy=dy,
        bkg_rad=bkg_rad,
        work_level=work_level,
        intense=intense,
        peak_relation=peak_relation,
        area_relation=area_relation,
        use_intense=use_intense,
    )

    return {
        "data": eclass.astype(np.int32),
        "standard_name": "echo_classification",
        "long_name": "Steiner echo classification",
        "valid_min": 0,
        "valid_max": 2,
        "comment_1": (
            "Convective-stratiform echo "
            "classification based on "
            "Steiner et al. (1995)"
        ),
        "comment_2": ("0 = Undefined, 1 = Stratiform, " "2 = Convective"),
    }


def conv_strat_yuter(
    grid,
    dx=None,
    dy=None,
    level_m=None,
    always_core_thres=42,
    bkg_rad_km=11,
    use_cosine=True,
    max_diff=5,
    zero_diff_cos_val=55,
    scalar_diff=1.5,
    use_addition=True,
    calc_thres=0.75,
    weak_echo_thres=5.0,
    min_dBZ_used=5.0,
    dB_averaging=True,
    remove_small_objects=True,
    min_km2_size=10,
    val_for_max_conv_rad=30,
    max_conv_rad_km=5.0,
    cs_core=3,
    nosfcecho=0,
    weakecho=3,
    sf=1,
    conv=2,
    refl_field=None,
    estimate_flag=True,
    estimate_offset=5,
):
    """
    Partition reflectivity into convective-stratiform using the Yuter et al. (2005)
    and Yuter and Houze (1997) algorithm.

    Parameters
    ----------
    grid : Grid
        Grid containing reflectivity field to partition.
    dx, dy : float, optional
        The x- and y-dimension resolutions in meters, respectively. If None
        the resolution is determined from the first two axes values parsed from grid object.
    level_m : float, optional
        Desired height in meters to classify with convective stratiform algorithm.
    always_core_thres : float, optional
        Threshold for points that are always convective. All values above the threshold are classifed as convective
        See Yuter et al. (2005) for more detail.
    bkg_rad_km : float, optional
        Radius to compute background reflectivity in kilometers. Default is 11 km. Recommended to be at least 3 x
        grid spacing
    use_cosine : bool, optional
        Boolean used to determine if a cosine scheme (see Yuter and Houze (1997)) should be used for identifying
        convective cores (True) or if a simpler scalar scheme (False) should be used.
    max_diff : float, optional
        Maximum difference between background average and reflectivity in order to be classified as convective.
        "a" value in Eqn. B1 in Yuter and Houze (1997)
    zero_diff_cos_val : float, optional
        Value where difference between background average and reflectivity is zero in the cosine function
        "b" value in Eqn. B1 in Yuter and Houze (1997)
    scalar_diff : float, optional
        If using a scalar difference scheme, this value is the multiplier or addition to the background average
    use_addition : bool, optional
        Determines if a multiplier (False) or addition (True) in the scalar difference scheme should be used
    calc_thres : float, optional
        Minimum percentage of points needed to be considered in background average calculation
    weak_echo_thres : float, optional
        Threshold for determining weak echo. All values below this threshold will be considered weak echo
    min_dBZ_used : float, optional
        Minimum dBZ value used for classification. All values below this threshold will be considered no surface echo
        See Yuter and Houze (1997) and Yuter et al. (2005) for more detail.
    dB_averaging : bool, optional
        True if using dBZ reflectivity values that need to be converted to linear Z before averaging. False for
        other non-dBZ values (i.e. snow rate)
    remove_small_objects : bool, optional
        Determines if small objects should be removed from convective core array. Default is True.
    min_km2_size : float, optional
        Minimum size of convective cores to be considered. Cores less than this size will be removed. Default is 10
        km^2.
    val_for_max_conv_rad : float, optional
        dBZ for maximum convective radius. Convective cores with values above this will have the maximum convective
        radius
    max_conv_rad_km : float, optional
        Maximum radius around convective cores to classify as convective. Default is 5 km
    cs_core : int, optional
        Value for points classified as convective cores
    nosfcecho : int, optional
        Value for points classified as no surface echo, based on min_dBZ_used
    weakecho : int, optional
        Value for points classified as weak echo, based on weak_echo_thres
    sf : int, optional
        Value for points classified as stratiform
    conv : int, optional
        Value for points classified as convective
    refl_field : str, optional
        Field in grid to use as the reflectivity during partitioning. None will use the default reflectivity
        field name from the Py-ART configuration file.
    estimate_flag : bool, optional
        Determines if over/underestimation should be applied. If true, the algorithm will also be run on the same field
        wih the estimate_offset added and the same field with the estimate_offset subtracted.
        Default is True (recommended)
    estimate_offset : float, optional
        Value used to offset the reflectivity values by for the over/underestimation application. Default value is 5
        dBZ.

    Returns
    -------
    convsf_dict : dict
        Convective-stratiform classification dictionary.

    References
    ----------
    Yuter, S. E., and R. A. Houze, Jr., 1997: Measurements of raindrop size
    distributions over the Pacific warm pool and implications for Z-R relations.
    J. Appl. Meteor., 36, 847-867.
    https://doi.org/10.1175/1520-0450(1997)036%3C0847:MORSDO%3E2.0.CO;2

    Yuter, S. E., R. A. Houze, Jr., E. A. Smith, T. T. Wilheit, and E. Zipser,
    2005: Physical characterization of tropical oceanic convection observed in
    KWAJEX. J. Appl. Meteor., 44, 385-415. https://doi.org/10.1175/JAM2206.1

    """

    # Maxmimum convective radius must be less than 5 km
    if max_conv_rad_km > 5:
        print("Max conv radius must be less than 5 km, exiting")
        raise

    # Parse field parameters
    if refl_field is None:
        refl_field = get_field_name("reflectivity")
        dB_averaging = True

    # parse dx and dy if None
    if dx is None:
        dx = grid.x["data"][1] - grid.x["data"][0]
    if dy is None:
        dy = grid.y["data"][1] - grid.y["data"][0]

    # add catch for background radius size
    if bkg_rad_km * 1000 < 2 * dx or bkg_rad_km * 1000 < 2 * dy:
        print(
            "Background radius for averaging must be at least 2 times dx and dy, exiting"
        )
        raise

    # Get coordinates
    z = grid.z["data"]

    # Get reflectivity data at desired level
    if level_m is None:
        try:
            ze = np.ma.copy(grid.fields[refl_field]["data"][0, :, :])
        except:
            ze = np.ma.copy(grid.fields[refl_field]["data"][:, :])
    else:
        zslice = np.argmin(np.abs(z - level_m))
        ze = np.ma.copy(grid.fields[refl_field]["data"][zslice, :, :])

    # run convective stratiform algorithm
    _, _, convsf_best = _revised_conv_strat(
        ze,
        dx,
        dy,
        always_core_thres=always_core_thres,
        bkg_rad_km=bkg_rad_km,
        use_cosine=use_cosine,
        max_diff=max_diff,
        zero_diff_cos_val=zero_diff_cos_val,
        scalar_diff=scalar_diff,
        use_addition=use_addition,
        calc_thres=calc_thres,
        weak_echo_thres=weak_echo_thres,
        min_dBZ_used=min_dBZ_used,
        dB_averaging=dB_averaging,
        remove_small_objects=remove_small_objects,
        min_km2_size=min_km2_size,
        val_for_max_conv_rad=val_for_max_conv_rad,
        max_conv_rad_km=max_conv_rad_km,
        cs_core=cs_core,
        nosfcecho=nosfcecho,
        weakecho=weakecho,
        sf=sf,
        conv=conv,
    )

    # put data into a dictionary to be added as a field
    convsf_dict = {
        "convsf": {
            "data": convsf_best,
            "standard_name": "convsf",
            "long_name": "Convective stratiform classification",
            "valid_min": 0,
            "valid_max": 3,
            "comment_1": (
                "Convective-stratiform echo "
                "classification based on "
                "Yuter and Houze (1997) and Yuter et al. (2005)"
            ),
            "comment_2": (
                "0 = No surface echo/Undefined, 1 = Stratiform, "
                "2 = Convective, 3 = weak echo"
            ),
        }
    }

    # If estimation is True, run the algorithm on the field with offset subtracted and the field with the offset added
    if estimate_flag:
        _, _, convsf_under = _revised_conv_strat(
            ze - estimate_offset,
            dx,
            dy,
            always_core_thres=always_core_thres,
            bkg_rad_km=bkg_rad_km,
            use_cosine=use_cosine,
            max_diff=max_diff,
            zero_diff_cos_val=zero_diff_cos_val,
            scalar_diff=scalar_diff,
            use_addition=use_addition,
            calc_thres=calc_thres,
            weak_echo_thres=weak_echo_thres,
            min_dBZ_used=min_dBZ_used,
            dB_averaging=dB_averaging,
            remove_small_objects=remove_small_objects,
            min_km2_size=min_km2_size,
            val_for_max_conv_rad=val_for_max_conv_rad,
            max_conv_rad_km=max_conv_rad_km,
            cs_core=cs_core,
            nosfcecho=nosfcecho,
            weakecho=weakecho,
            sf=sf,
            conv=conv,
        )

        _, _, convsf_over = _revised_conv_strat(
            ze + estimate_offset,
            dx,
            dy,
            always_core_thres=always_core_thres,
            bkg_rad_km=bkg_rad_km,
            use_cosine=use_cosine,
            max_diff=max_diff,
            zero_diff_cos_val=zero_diff_cos_val,
            scalar_diff=scalar_diff,
            use_addition=use_addition,
            calc_thres=calc_thres,
            weak_echo_thres=weak_echo_thres,
            min_dBZ_used=min_dBZ_used,
            dB_averaging=dB_averaging,
            remove_small_objects=remove_small_objects,
            min_km2_size=min_km2_size,
            val_for_max_conv_rad=val_for_max_conv_rad,
            max_conv_rad_km=max_conv_rad_km,
            cs_core=cs_core,
            nosfcecho=nosfcecho,
            weakecho=weakecho,
            sf=sf,
            conv=conv,
        )

        # save into dictionaries
        convsf_dict["convsf_under"] = {
            "data": convsf_under,
            "standard_name": "convsf_under",
            "long_name": "Convective stratiform classification (Underestimate)",
            "valid_min": 0,
            "valid_max": 3,
            "comment_1": (
                "Convective-stratiform echo "
                "classification based on "
                "Yuter and Houze (1997) and Yuter et al. (2005)"
            ),
            "comment_2": (
                "0 = Undefined, 1 = Stratiform, " "2 = Convective, 3 = weak echo"
            ),
        }

        convsf_dict["convsf_over"] = {
            "data": convsf_over,
            "standard_name": "convsf_under",
            "long_name": "Convective stratiform classification (Overestimate)",
            "valid_min": 0,
            "valid_max": 3,
            "comment_1": (
                "Convective-stratiform echo "
                "classification based on "
                "Yuter and Houze (1997)"
            ),
            "comment_2": (
                "0 = Undefined, 1 = Stratiform, " "2 = Convective, 3 = weak echo"
            ),
        }

    return convsf_dict


def hydroclass_semisupervised(
    radar,
    mass_centers=None,
    weights=np.array([1.0, 1.0, 1.0, 0.75, 0.5]),
    refl_field=None,
    zdr_field=None,
    rhv_field=None,
    kdp_field=None,
    temp_field=None,
    hydro_field=None,
):
    """
    Classifies precipitation echoes following the approach by Besic et al
    (2016).

    Parameters
    ----------
    radar : radar
        Radar object.
    mass_centers : ndarray 2D, optional
        The centroids for each variable and hydrometeor class in (nclasses,
        nvariables).
    weights : ndarray 1D, optional
        The weight given to each variable.
    refl_field, zdr_field, rhv_field, kdp_field, temp_field : str, optional
        Inputs. Field names within the radar object which represent the
        horizonal reflectivity, the differential reflectivity, the copolar
        correlation coefficient, the specific differential phase and the
        temperature field. A value of None for any of these parameters will
        use the default field name as defined in the Py-ART configuration
        file.
    hydro_field : str, optional
        Output. Field name which represents the hydrometeor class field.
        A value of None will use the default field name as defined in the
        Py-ART configuration file.

    Returns
    -------
    hydro : dict
        Hydrometeor classification field.

    References
    ----------
    Besic, N., Figueras i Ventura, J., Grazioli, J., Gabella, M., Germann, U.,
    and Berne, A.: Hydrometeor classification through statistical clustering
    of polarimetric radar measurements: a semi-supervised approach,
    Atmos. Meas. Tech., 9, 4425-4445, doi:10.5194/amt-9-4425-2016, 2016

    """
    lapse_rate = -6.5

    # select the centroids as a function of frequency band
    if mass_centers is None:
        # assign coefficients according to radar frequency
        if radar.instrument_parameters is not None:
            if "frequency" in radar.instrument_parameters:
                mass_centers = _get_mass_centers(
                    radar.instrument_parameters["frequency"]["data"][0]
                )
            else:
                mass_centers = _mass_centers_table()["C"]
                warn(
                    "Radar frequency unknown. "
                    "Default coefficients for C band will be applied."
                )
        else:
            mass_centers = _mass_centers_table()["C"]
            warn(
                "Radar instrument parameters is empty. So frequency is "
                "unknown. Default coefficients for C band will be applied."
            )

    # parse the field parameters
    if refl_field is None:
        refl_field = get_field_name("reflectivity")
    if zdr_field is None:
        zdr_field = get_field_name("differential_reflectivity")
    if rhv_field is None:
        rhv_field = get_field_name("cross_correlation_ratio")
    if kdp_field is None:
        kdp_field = get_field_name("specific_differential_phase")
    if temp_field is None:
        temp_field = get_field_name("temperature")
    if hydro_field is None:
        hydro_field = get_field_name("radar_echo_classification")

    # extract fields and parameters from radar
    radar.check_field_exists(refl_field)
    radar.check_field_exists(zdr_field)
    radar.check_field_exists(rhv_field)
    radar.check_field_exists(kdp_field)
    radar.check_field_exists(temp_field)

    refl = radar.fields[refl_field]["data"]
    zdr = radar.fields[zdr_field]["data"]
    rhohv = radar.fields[rhv_field]["data"]
    kdp = radar.fields[kdp_field]["data"]
    temp = radar.fields[temp_field]["data"]

    # convert temp in relative height respect to iso0
    relh = temp * (1000.0 / lapse_rate)

    # standardize data
    refl_std = _standardize(refl, "Zh")
    zdr_std = _standardize(zdr, "ZDR")
    kdp_std = _standardize(kdp, "KDP")
    rhohv_std = _standardize(rhohv, "RhoHV")
    relh_std = _standardize(relh, "relH")

    # standardize centroids
    mc_std = np.zeros(np.shape(mass_centers))
    mc_std[:, 0] = _standardize(mass_centers[:, 0], "Zh")
    mc_std[:, 1] = _standardize(mass_centers[:, 1], "ZDR")
    mc_std[:, 2] = _standardize(mass_centers[:, 2], "KDP")
    mc_std[:, 3] = _standardize(mass_centers[:, 3], "RhoHV")
    mc_std[:, 4] = _standardize(mass_centers[:, 4], "relH")

    # assign to class
    hydroclass_data, min_dist = _assign_to_class(
        refl_std, zdr_std, kdp_std, rhohv_std, relh_std, mc_std, weights=weights
    )

    # prepare output fields
    hydro = get_metadata(hydro_field)
    hydro["data"] = hydroclass_data

    return hydro


def _standardize(data, field_name, mx=None, mn=None):
    """
    Streches the radar data to -1 to 1 interval.

    Parameters
    ----------
    data : array
        Radar field.
    field_name : str
        Type of field (relH, Zh, ZDR, KDP or RhoHV).
    mx, mn : floats or None, optional
        Data limits for array values.

    Returns
    -------
    field_std : dict
        Standardized radar data.

    """
    if field_name == "relH":
        field_std = 2.0 / (1.0 + np.ma.exp(-0.005 * data)) - 1.0
        return field_std

    if (mx is None) or (mn is None):
        dlimits_dict = _data_limits_table()
        if field_name not in dlimits_dict:
            raise ValueError(
                "Field "
                + field_name
                + " unknown. "
                + "Valid field names for standardizing are: "
                + "relH, Zh, ZDR, KDP and RhoHV"
            )
        mx, mn = dlimits_dict[field_name]

    if field_name == "KDP":
        data[data < -0.5] = -0.5
        data = 10.0 * np.ma.log10(data + 0.6)
    elif field_name == "RhoHV":
        data = 10.0 * np.ma.log10(1.0 - data)

    mask = np.ma.getmaskarray(data)
    field_std = 2.0 * (data - mn) / (mx - mn) - 1.0
    field_std[data < mn] = -1.0
    field_std[data > mx] = 1.0
    field_std[mask] = np.ma.masked

    return field_std


def _assign_to_class(
    zh,
    zdr,
    kdp,
    rhohv,
    relh,
    mass_centers,
    weights=np.array([1.0, 1.0, 1.0, 0.75, 0.5]),
):
    """
    Assigns an hydrometeor class to a radar range bin computing
    the distance between the radar variables an a centroid.

    Parameters
    ----------
    zh, zdr, kdp, rhohv, relh : radar fields
        Variables used for assignment normalized to [-1, 1] values.
    mass_centers : matrix
        Centroids normalized to [-1, 1] values.
    weights : array, optional
        The weight given to each variable.

    Returns
    -------
    hydroclass : int array
        The index corresponding to the assigned class.
    mind_dist : float array
        The minimum distance to the centroids.

    """
    # prepare data
    nrays = zh.shape[0]
    nbins = zdr.shape[1]
    nclasses = mass_centers.shape[0]
    nvariables = mass_centers.shape[1]

    data = np.ma.array([zh, zdr, kdp, rhohv, relh])
    weights_mat = np.broadcast_to(
        weights.reshape(nvariables, 1, 1), (nvariables, nrays, nbins)
    )
    dist = np.ma.zeros((nclasses, nrays, nbins), dtype="float64")

    # compute distance: masked entries will not contribute to the distance
    for i in range(nclasses):
        centroids_class = mass_centers[i, :]
        centroids_class = np.broadcast_to(
            centroids_class.reshape(nvariables, 1, 1), (nvariables, nrays, nbins)
        )
        dist[i, :, :] = np.ma.sqrt(
            np.ma.sum(((centroids_class - data) ** 2.0) * weights_mat, axis=0)
        )

    # use very large fill_value so that masked entries will be sorted at the
    # end. There should not be any masked entry anyway
    class_vec = dist.argsort(axis=0, fill_value=10e40)

    # get minimum distance. Acts as a confidence value
    dist.sort(axis=0, fill_value=10e40)
    min_dist = dist[0, :, :]

    # Entries with non-valid reflectivity values are set to 0 (No class)
    mask = np.ma.getmaskarray(zh)
    hydroclass = class_vec[0, :, :] + 1
    hydroclass[mask] = 0

    return hydroclass, min_dist


def _get_mass_centers(freq):
    """
    Get mass centers for a particular frequency.

    Parameters
    ----------
    freq : float
        Radar frequency [Hz].

    Returns
    -------
    mass_centers : ndarray 2D
        The centroids for each variable and hydrometeor class in (nclasses,
        nvariables).

    """
    mass_centers_dict = _mass_centers_table()

    freq_band = get_freq_band(freq)
    if (freq_band is not None) and (freq_band in mass_centers_dict):
        return mass_centers_dict[freq_band]

    if freq < 4e9:
        freq_band_aux = "C"
    elif freq > 12e9:
        freq_band_aux = "X"

    mass_centers = mass_centers_dict[freq_band_aux]
    warn(
        "Radar frequency out of range. "
        + "Centroids only valid for C or X band. "
        + freq_band_aux
        + " band centroids will be applied"
    )

    return mass_centers


def _mass_centers_table():
    """
    Defines the mass centers look up table for each frequency band.

    Returns
    -------
    mass_centers_dict : dict
        A dictionary with the mass centers for each frequency band.

    """
    nclasses = 9
    nvariables = 5
    mass_centers_c = np.zeros((nclasses, nvariables))
    mass_centers_x = np.zeros((nclasses, nvariables))

    mass_centers_dict = dict()
    # C-band centroids derived for MeteoSwiss Albis radar
    #                       Zh        ZDR     kdp   RhoHV    delta_Z
    mass_centers_c[0, :] = [13.5829, 0.4063, 0.0497, 0.9868, 1330.3]  # DS
    mass_centers_c[1, :] = [02.8453, 0.2457, 0.0000, 0.9798, 0653.8]  # CR
    mass_centers_c[2, :] = [07.6597, 0.2180, 0.0019, 0.9799, -1426.5]  # LR
    mass_centers_c[3, :] = [31.6815, 0.3926, 0.0828, 0.9978, 0535.3]  # GR
    mass_centers_c[4, :] = [39.4703, 1.0734, 0.4919, 0.9876, -1036.3]  # RN
    mass_centers_c[5, :] = [04.8267, -0.5690, 0.0000, 0.9691, 0869.8]  # VI
    mass_centers_c[6, :] = [30.8613, 0.9819, 0.1998, 0.9845, -0066.1]  # WS
    mass_centers_c[7, :] = [52.3969, 2.1094, 2.4675, 0.9730, -1550.2]  # MH
    mass_centers_c[8, :] = [50.6186, -0.0649, 0.0946, 0.9904, 1179.9]  # IH/HDG

    mass_centers_dict.update({"C": mass_centers_c})

    # X-band centroids derived for MeteoSwiss DX50 radar
    #                       Zh        ZDR     kdp    RhoHV   delta_Z
    mass_centers_x[0, :] = [19.0770, 0.4139, 0.0099, 0.9841, 1061.7]  # DS
    mass_centers_x[1, :] = [03.9877, 0.5040, 0.0000, 0.9642, 0856.6]  # CR
    mass_centers_x[2, :] = [20.7982, 0.3177, 0.0004, 0.9858, -1375.1]  # LR
    mass_centers_x[3, :] = [34.7124, -0.3748, 0.0988, 0.9828, 1224.2]  # GR
    mass_centers_x[4, :] = [33.0134, 0.6614, 0.0819, 0.9802, -1169.8]  # RN
    mass_centers_x[5, :] = [08.2610, -0.4681, 0.0000, 0.9722, 1100.7]  # VI
    mass_centers_x[6, :] = [35.1801, 1.2830, 0.1322, 0.9162, -0159.8]  # WS
    mass_centers_x[7, :] = [52.4539, 2.3714, 1.1120, 0.9382, -1618.5]  # MH
    mass_centers_x[8, :] = [44.2216, -0.3419, 0.0687, 0.9683, 1272.7]  # IH/HDG

    mass_centers_dict.update({"X": mass_centers_x})

    return mass_centers_dict


def _data_limits_table():
    """
    Defines the data limits used in the standardization.

    Returns
    -------
    dlimits_dict : dict
        A dictionary with the limits for each variable.

    """
    dlimits_dict = dict()
    dlimits_dict.update({"Zh": (60.0, -10.0)})
    dlimits_dict.update({"ZDR": (5.0, -5.0)})
    dlimits_dict.update({"KDP": (7.0, -10.0)})
    dlimits_dict.update({"RhoHV": (-5.23, -50.0)})

    return dlimits_dict


def get_freq_band(freq):
    """
    Returns the frequency band name (S, C, X, ...).

    Parameters
    ----------
    freq : float
        Radar frequency [Hz].

    Returns
    -------
    freq_band : str
        Frequency band name.

    """
    if freq >= 2e9 and freq < 4e9:
        return "S"
    if freq >= 4e9 and freq < 8e9:
        return "C"
    if freq >= 8e9 and freq <= 12e9:
        return "X"

    warn("Unknown frequency band")

    return None
