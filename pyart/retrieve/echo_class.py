"""
Functions for echo classification.

"""

from warnings import warn

import numpy as np

# Local imports
from ..config import get_field_name, get_fillvalue, get_metadata
from ..core import Grid
from ._echo_class import _feature_detection, steiner_class_buff
from ._echo_class_wt import calc_scale_break, wavelet_reclass


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
    ze = ze.filled(np.nan)

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

    warn(
        "This function will be deprecated in Py-ART 2.0."
        " Please use feature_detection function.",
        DeprecationWarning,
    )

    feature_dict = feature_detection(
        grid,
        dx=dx,
        dy=dy,
        level_m=level_m,
        always_core_thres=always_core_thres,
        bkg_rad_km=bkg_rad_km,
        use_cosine=use_cosine,
        max_diff=max_diff,
        zero_diff_cos_val=zero_diff_cos_val,
        scalar_diff=scalar_diff,
        use_addition=use_addition,
        calc_thres=calc_thres,
        weak_echo_thres=weak_echo_thres,
        min_val_used=min_dBZ_used,
        dB_averaging=dB_averaging,
        remove_small_objects=remove_small_objects,
        min_km2_size=min_km2_size,
        binary_close=False,
        val_for_max_rad=val_for_max_conv_rad,
        max_rad_km=max_conv_rad_km,
        core_val=cs_core,
        nosfcecho=nosfcecho,
        weakecho=weakecho,
        bkgd_val=sf,
        feat_val=conv,
        field=refl_field,
        estimate_flag=estimate_flag,
        estimate_offset=5,
        overest_field=None,
        underest_field=None,
    )

    return feature_dict


def feature_detection(
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
    min_val_used=5.0,
    dB_averaging=True,
    remove_small_objects=True,
    min_km2_size=10,
    binary_close=False,
    val_for_max_rad=30,
    max_rad_km=5.0,
    core_val=3,
    nosfcecho=0,
    weakecho=3,
    bkgd_val=1,
    feat_val=2,
    field=None,
    estimate_flag=True,
    estimate_offset=5,
    overest_field=None,
    underest_field=None,
):
    """
    This function can be used to detect features in a field (e.g. reflectivity, rain rate, snow rate,
    etc.) described by Tomkins et al. (2023) and based on original convective-stratiform algorithms developed by
    Steiner et al. (1995), Yuter et al. (2005) and Yuter and Houze (1997) algorithm.

    Author: Laura Tomkins (@lauratomkins)

    Parameters
    ----------
    grid : Grid
        Grid containing reflectivity field to partition.
    dx, dy : float, optional
        The x- and y-dimension resolutions in meters, respectively. If None
        the resolution is determined from the first two axes values parsed from grid object.
    level_m : float, optional
        Desired height in meters to run feature detection algorithm.
    always_core_thres : float, optional
        Threshold for points that are always features. All values above the threshold are classified as features.
    bkg_rad_km : float, optional
        Radius to compute background reflectivity in kilometers. Default is 11 km. Recommended to be at least 3 x
        grid spacing
    use_cosine : bool, optional
        Boolean used to determine if a cosine scheme (see Yuter and Houze (1997)) should be used for identifying
        cores (True) or if a simpler scalar scheme (False) should be used.
    max_diff : float, optional
        Maximum difference between background average and reflectivity in order to be classified as features.
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
    min_val_used : float, optional
        Minimum value used for classification. All values below this threshold will be considered no surface echo
        See Yuter and Houze (1997) and Yuter et al. (2005) for more detail. Units based on input field
    dB_averaging : bool, optional
        True if using dBZ reflectivity values that need to be converted to linear Z before averaging. False for
        other non-dBZ values (i.e. snow rate)
    remove_small_objects : bool, optional
        Determines if small objects should be removed from core array. Default is True.
    min_km2_size : float, optional
        Minimum size of Cores to be considered. Cores less than this size will be removed. Default is 10 km^2.
    binary_close : bool, optional
        Determines if a binary closing should be performed on the cores. Default is False.
    val_for_max_rad : float, optional
        value used for maximum radius. Cores with values above this will have the maximum radius incorporated.
    max_rad_km : float, optional
        Maximum radius around cores to classify as feature. Default is 5 km
    core_val : int, optional
        Value for points classified as cores
    nosfcecho : int, optional
        Value for points classified as no surface echo, based on min_val_used
    weakecho : int, optional
        Value for points classified as weak echo, based on weak_echo_thres.
    bkgd_val : int, optional
        Value for points classified as background echo.
    feat_val : int, optional
        Value for points classified as features.
    field : str, optional
        Field in grid to find objects in. None will use the default reflectivity field name from the Py-ART
        configuration file.
    estimate_flag : bool, optional
        Determines if over/underestimation should be applied. If true, the algorithm will also be run on the same field
        wih the estimate_offset added and the same field with the estimate_offset subtracted.
        Default is True (recommended)
    estimate_offset : float, optional
        Value used to offset the field values by for the over/underestimation application. Default value is 5 dBZ.
    overest_field : str, optional
        Name of field in grid object used to calculate the overestimate if estimate_flag is True.
    underest_field : str, optional
        Name of field in grid object used to calculate the underestimate if estimate_flag is True.

    Returns
    -------
    feature_dict : dict
        Feature detection classification dictionary.

    References
    ----------
    Steiner, M. R., R. A. Houze Jr., and S. E. Yuter, 1995: Climatological
    Characterization of Three-Dimensional Storm Structure from Operational
    Radar and Rain Gauge Data. J. Appl. Meteor., 34, 1978-2007.

    Yuter, S. E., and R. A. Houze, Jr., 1997: Measurements of raindrop size
    distributions over the Pacific warm pool and implications for Z-R relations.
    J. Appl. Meteor., 36, 847-867.
    https://doi.org/10.1175/1520-0450(1997)036%3C0847:MORSDO%3E2.0.CO;2

    Yuter, S. E., R. A. Houze, Jr., E. A. Smith, T. T. Wilheit, and E. Zipser,
    2005: Physical characterization of tropical oceanic convection observed in
    KWAJEX. J. Appl. Meteor., 44, 385-415. https://doi.org/10.1175/JAM2206.1

    Tomkins, L. M., S. E. Yuter, and M. A. Miller, 2024: Objective identification
    of faint and strong features in radar observations of winter storms. in prep.

    """

    # Maxmimum radius must be less than 5 km
    if max_rad_km > 5:
        print("Max radius must be less than 5 km, exiting")
        raise

    # Parse field parameters
    if field is None:
        field = get_field_name("reflectivity")
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
            ze = np.ma.copy(grid.fields[field]["data"][0, :, :])
        except:
            ze = np.ma.copy(grid.fields[field]["data"][:, :])
    else:
        zslice = np.argmin(np.abs(z - level_m))
        ze = np.ma.copy(grid.fields[field]["data"][zslice, :, :])

    # run feature detection algorithm
    _, _, feature_best = _feature_detection(
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
        min_val_used=min_val_used,
        dB_averaging=dB_averaging,
        remove_small_objects=remove_small_objects,
        min_km2_size=min_km2_size,
        binary_close=binary_close,
        val_for_max_rad=val_for_max_rad,
        max_rad_km=max_rad_km,
        core_val=core_val,
        nosfcecho=nosfcecho,
        weakecho=weakecho,
        bkgd_val=bkgd_val,
        feat_val=feat_val,
    )

    # put data into a dictionary to be added as a field
    feature_dict = {
        "feature_detection": {
            "data": feature_best[None, :, :],
            "standard_name": "feature_detection",
            "long_name": "Feature Detection",
            "valid_min": 0,
            "valid_max": 3,
            "comment_1": (
                f"{nosfcecho} = No surface echo/Undefined, {bkgd_val} = Background echo, {feat_val} = Features, {weakecho} = weak echo"
            ),
        }
    }

    # If estimation is True, run the algorithm on the field with offset subtracted and the field with the offset added
    if estimate_flag:
        # get underestimate field
        if underest_field is None:
            under_field = ze - estimate_offset
        elif underest_field is not None:
            under_field = np.ma.copy(grid.fields[underest_field]["data"][0, :, :])

        # run algorithm to get feature detection underestimate
        _, _, feature_under = _feature_detection(
            under_field,
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
            min_val_used=min_val_used,
            dB_averaging=dB_averaging,
            remove_small_objects=remove_small_objects,
            min_km2_size=min_km2_size,
            binary_close=binary_close,
            val_for_max_rad=val_for_max_rad,
            max_rad_km=max_rad_km,
            core_val=core_val,
            nosfcecho=nosfcecho,
            weakecho=weakecho,
            bkgd_val=bkgd_val,
            feat_val=feat_val,
        )

        # get overestimate field
        if overest_field is None:
            over_field = ze + estimate_offset
        elif overest_field is not None:
            over_field = np.ma.copy(grid.fields[overest_field]["data"][0, :, :])

        # run algorithm to get feature detection underestimate
        _, _, feature_over = _feature_detection(
            over_field,
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
            min_val_used=min_val_used,
            dB_averaging=dB_averaging,
            remove_small_objects=remove_small_objects,
            min_km2_size=min_km2_size,
            binary_close=binary_close,
            val_for_max_rad=val_for_max_rad,
            max_rad_km=max_rad_km,
            core_val=core_val,
            nosfcecho=nosfcecho,
            weakecho=weakecho,
            bkgd_val=bkgd_val,
            feat_val=feat_val,
        )

        # save into dictionaries
        feature_dict["feature_under"] = {
            "data": feature_under[None, :, :],
            "standard_name": "feature_detection_under",
            "long_name": "Feature Detection Underestimate",
            "valid_min": 0,
            "valid_max": 3,
            "comment_1": (
                f"{nosfcecho} = No surface echo/Undefined, {bkgd_val} = Background echo, {feat_val} = Features, {weakecho} = weak echo"
            ),
        }

        feature_dict["feature_over"] = {
            "data": feature_over[None, :, :],
            "standard_name": "feature_detection_over",
            "long_name": "Feature Detection Overestimate",
            "valid_min": 0,
            "valid_max": 3,
            "comment_1": (
                f"{nosfcecho} = No surface echo/Undefined, {bkgd_val} = Background echo, {feat_val} = Features, {weakecho} = weak echo"
            ),
        }

    return feature_dict


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
    radar_freq=None,
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
    radar_freq : str, optional
        Radar frequency in Hertz (Hz) used for classification.
        This parameter will be ignored, if the radar object has frequency information.

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

    Notes
    -----
    The default hydrometeor classification is valid for C-band radars. For X-band radars,
    if frequency information is not present in the `radar.instrument_parameters`, the user-supplied
    `radar_freq` will be used with a warning. If both `radar.instrument_parameters` and
    `radar_freq` parameter are missing, the algorithm defaults to the C band.

    If the radar frequency information is missing from the radar object, you can add it in
    `radar.instrument_parameters`, as follows:
    .. code-block:: python

        radar.instrument_parameters["frequency"] = {
            "long_name": "Radar frequency",
            "units": "Hz",
            "data": [9.2e9]
        }
    """
    lapse_rate = -6.5

    # select the centroids as a function of frequency band
    if mass_centers is None:
        # assign coefficients according to radar frequency
        if radar.instrument_parameters and "frequency" in radar.instrument_parameters:
            frequency = radar.instrument_parameters["frequency"]["data"][0]
            mass_centers = _get_mass_centers(frequency)
            warn(f"Using radar frequency from instrument parameters: {frequency}")
        elif radar_freq is not None:
            mass_centers = _get_mass_centers(radar_freq)
            warn(
                f"Radar instrument parameters are empty. Using user-supplied radar frequency: {radar_freq}"
            )
        else:
            mass_centers = _mass_centers_table()["C"]
            warn(
                "Radar instrument parameters and radar_freq param are empty."
                "So frequency is unknown. Default coefficients for C band will be applied."
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


def conv_strat_raut(
    grid,
    refl_field,
    cappi_level=0,
    zr_a=200,
    zr_b=1.6,
    core_wt_threshold=5,
    conv_wt_threshold=1.5,
    conv_scale_km=25,
    min_reflectivity=5,
    conv_min_refl=25,
    conv_core_threshold=42,
    override_checks=False,
):
    """
    A computationally efficient method to classify radar echoes into convective cores, mixed convection,
    and stratiform regions for gridded radar reflectivity field.

    This function uses à trous wavelet transform (ATWT) for multiresolution (i.e. scale) analysis of radar field,
    focusing on precipitation structure over reflectivity thresholds for robust echo classification (Raut et al 2008, 2020).

    Parameters
    ----------
    grid : PyART Grid
        Grid object containing radar data.
    refl_field : str
        Field name for reflectivity data in the Py-ART grid object.
    zr_a : float, optional
        Coefficient 'a' in the Z-R relationship Z = a*R^b for reflectivity to rain rate conversion.
        The algorithm is not sensitive to precise values of 'zr_a' and 'zr_b'; however,
        they must be adjusted based on the type of radar used.
        Default is 200.
    zr_b : float, optional
        Coefficient 'b' in the Z-R relationship Z = a*R^b. Default is 1.6.
    core_wt_threshold : float, optional
        Threshold for wavelet components to separate convective cores from mix-intermediate type.
        Default is 5. Recommended values are between 4 and 6.
    conv_wt_threshold : float, optional
        Threshold for significant wavelet components to separate all convection from stratiform.
        Default is 1.5. Recommended values are between 1 and 2.
    conv_scale_km : float, optional
        Approximate scale break (in km) between convective and stratiform scales.
        Scale break may vary over different regions and seasons
        (Refere to Raut et al 2018 for more discussion on scale-breaks). Note that the
        algorithm is insensitive to small variations in the scale break due to the
        dyadic nature of the scaling. The actual scale break used in the calculation of wavelets
        is returned in the output dictionary by parameter `scale_break_used`.
        Default is 25 km. Recommended values are between 16 and 32 km.
    min_reflectivity : float, optional
        Minimum reflectivity threshold. Reflectivities below this value are not classified.
        Default is 5 dBZ. This value must be greater than or equal to '0'.
    conv_min_refl : float, optional
        Reflectivity values lower than this threshold will be always considered as non-convective.
        Default is 25 dBZ. Recommended values are between 25 and 30 dBZ.
    conv_core_threshold : float, optional
        Reflectivities above this threshold are classified as convective cores if wavelet components are significant (See: conv_wt_threshold).
        Default is 42 dBZ.
        Recommended value must be is greater than or equal to 40 dBZ. The algorithm is not sensitive to this value.
    override_checks : bool, optional
        If set to True, the function will bypass the sanity checks for above parameter values.
        This allows the user to use custom values for parameters, even if they fall outside
        the recommended ranges. The default is False.

    Returns
    -------

    dict :
    A dictionary structured as a Py-ART grid field, suitable for adding to a Py-ART Grid object. The dictionary
    contains the classification data and associated metadata. The classification categories are as follows:
        - 3: Convective Cores: associated with strong updrafts and active collision-coalescence.
        - 2: Mixed-Intermediate: capturing a wide range of convective activities, excluding the convective cores.
        - 1: Stratiform: remaining areas with more uniform and less intense precipitation.
        - 0: Unclassified: for reflectivity below the minimum threshold.


    References
    ----------
    Raut, B. A., Karekar, R. N., & Puranik, D. M. (2008). Wavelet-based technique to extract convective clouds from
    infrared satellite images. IEEE Geosci. Remote Sens. Lett., 5(3), 328-330.

    Raut, B. A., Seed, A. W., Reeder, M. J., & Jakob, C. (2018). A multiplicative cascade model for high‐resolution
    space‐time downscaling of rainfall. J. Geophys. Res. Atmos., 123(4), 2050-2067.

    Raut, B. A., Louf, V., Gayatri, K., Murugavel, P., Konwar, M., & Prabhakaran, T. (2020). A multiresolution technique
    for the classification of precipitation echoes in radar data. IEEE Trans. Geosci. Remote Sens., 58(8), 5409-5415.
    """

    # Check if the grid is a Py-ART Grid object
    if not isinstance(grid, Grid):
        raise TypeError("The 'grid' is not a Py-ART Grid object.")

    # Check if dx and dy are the same, and warn if not
    dx = grid.x["data"][1] - grid.x["data"][0]
    dy = grid.y["data"][1] - grid.y["data"][0]
    if dx != dy:
        warn(
            "Warning: Grid resolution `dx` and `dy` should be comparable for correct results.",
            UserWarning,
        )

    # Compute scale break (dyadic) here to paas it on as parameter to user dictionary
    scale_break = calc_scale_break(res_meters=dx, conv_scale_km=conv_scale_km)

    # From dyadic scale to km
    scale_break_km = (2 ** (scale_break - 1)) * dx / 1000

    # Sanity checks for parameters if override_checks is False
    if not override_checks:
        conv_core_threshold = max(
            40, conv_core_threshold
        )  # Ensure conv_core_threshold is at least 40 dBZ
        core_wt_threshold = max(
            4, min(core_wt_threshold, 6)
        )  # core_wt_threshold should be between 4 and 6
        conv_wt_threshold = max(
            1, min(conv_wt_threshold, 2)
        )  # conv_wt_threshold should be between 1 and 2
        conv_scale_km = max(
            16, min(conv_scale_km, 32)
        )  # conv_scale_km should be between 15 and 30 km
        min_reflectivity = max(
            0, min_reflectivity
        )  # min_reflectivity should be non-negative
        conv_min_refl = max(
            25, min(conv_min_refl, 30)
        )  # conv_min_refl should be between 25 and 30 dBZ

    # Call the actual wavelet_relass function to obtain radar echo classificatino
    reclass = wavelet_reclass(
        grid,
        refl_field,
        cappi_level,
        zr_a,
        zr_b,
        core_wt_threshold=core_wt_threshold,
        conv_wt_threshold=conv_wt_threshold,
        scale_break=scale_break,
        min_reflectivity=min_reflectivity,
        conv_min_refl=conv_min_refl,
        conv_core_threshold=conv_core_threshold,
    )

    reclass = np.expand_dims(reclass, axis=0)

    # put data into a dictionary to be added as a field
    reclass_dict = {
        "wt_reclass": {
            "data": reclass,
            "standard_name": "wavelet_echo_class",
            "long_name": "Wavelet-based multiresolution radar echo classification",
            "valid_min": 0,
            "valid_max": 3,
            "classification_description": "0: Unclassified, 1: Stratiform, 2: Mixed-Intermediate, 3: Convective Cores",
            "parameters": {
                "refl_field": refl_field,
                "cappi_level": cappi_level,
                "zr_a": zr_a,
                "zr_b": zr_b,
                "core_wt_threshold": core_wt_threshold,
                "conv_wt_threshold": conv_wt_threshold,
                "conv_scale_km": conv_scale_km,
                "scale_break_used": int(scale_break_km),
                "min_reflectivity": min_reflectivity,
                "conv_min_refl": conv_min_refl,
                "conv_core_threshold": conv_core_threshold,
            },
        }
    }

    return reclass_dict
