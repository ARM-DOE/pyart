import numpy as np
import scipy.ndimage


def _steiner_conv_strat(
    refl,
    x,
    y,
    dx,
    dy,
    intense=42,
    peak_relation=0,
    area_relation=1,
    bkg_rad=11000,
    use_intense=True,
):
    """
    We perform the Steiner et al. (1995) algorithm for echo classification
    using only the reflectivity field in order to classify each grid point
    as either convective, stratiform or undefined. Grid points are
    classified as follows,

    0 = Undefined
    1 = Stratiform
    2 = Convective
    """

    def convective_radius(ze_bkg, area_relation):
        """
        Given a mean background reflectivity value, we determine via a step
        function what the corresponding convective radius would be.

        Higher background reflectivitives are expected to have larger
        convective influence on surrounding areas, so a larger convective
        radius would be prescribed.
        """
        if area_relation == 0:
            if ze_bkg < 30:
                conv_rad = 1000.0
            elif (ze_bkg >= 30) & (ze_bkg < 35.0):
                conv_rad = 2000.0
            elif (ze_bkg >= 35.0) & (ze_bkg < 40.0):
                conv_rad = 3000.0
            elif (ze_bkg >= 40.0) & (ze_bkg < 45.0):
                conv_rad = 4000.0
            else:
                conv_rad = 5000.0

        if area_relation == 1:
            if ze_bkg < 25:
                conv_rad = 1000.0
            elif (ze_bkg >= 25) & (ze_bkg < 30.0):
                conv_rad = 2000.0
            elif (ze_bkg >= 30.0) & (ze_bkg < 35.0):
                conv_rad = 3000.0
            elif (ze_bkg >= 35.0) & (ze_bkg < 40.0):
                conv_rad = 4000.0
            else:
                conv_rad = 5000.0

        if area_relation == 2:
            if ze_bkg < 20:
                conv_rad = 1000.0
            elif (ze_bkg >= 20) & (ze_bkg < 25.0):
                conv_rad = 2000.0
            elif (ze_bkg >= 25.0) & (ze_bkg < 30.0):
                conv_rad = 3000.0
            elif (ze_bkg >= 30.0) & (ze_bkg < 35.0):
                conv_rad = 4000.0
            else:
                conv_rad = 5000.0

        if area_relation == 3:
            if ze_bkg < 40:
                conv_rad = 0.0
            elif (ze_bkg >= 40) & (ze_bkg < 45.0):
                conv_rad = 1000.0
            elif (ze_bkg >= 45.0) & (ze_bkg < 50.0):
                conv_rad = 2000.0
            elif (ze_bkg >= 50.0) & (ze_bkg < 55.0):
                conv_rad = 6000.0
            else:
                conv_rad = 8000.0

        return conv_rad

    def peakedness(ze_bkg, peak_relation):
        """
        Given a background reflectivity value, we determine what the necessary
        peakedness (or difference) has to be between a grid point's
        reflectivity and the background reflectivity in order for that grid
        point to be labeled convective.
        """
        if peak_relation == 0:
            if ze_bkg < 0.0:
                peak = 10.0
            elif (ze_bkg >= 0.0) and (ze_bkg < 42.43):
                peak = 10.0 - ze_bkg**2 / 180.0
            else:
                peak = 0.0

        elif peak_relation == 1:
            if ze_bkg < 0.0:
                peak = 14.0
            elif (ze_bkg >= 0.0) and (ze_bkg < 42.43):
                peak = 14.0 - ze_bkg**2 / 180.0
            else:
                peak = 4.0

        return peak

    sclass = np.zeros(refl.shape, dtype=int)
    ny, nx = refl.shape

    for i in range(0, nx):
        # Get stencil of x grid points within the background radius
        imin = np.max(np.array([1, (i - bkg_rad / dx)], dtype=int))
        imax = np.min(np.array([nx, (i + bkg_rad / dx)], dtype=int))

        for j in range(0, ny):
            # First make sure that the current grid point has not already been
            # classified. This can happen when grid points within the
            # convective radius of a previous grid point have also been
            # classified.
            if ~np.isnan(refl[j, i]) & (sclass[j, i] == 0):
                # Get stencil of y grid points within the background radius
                jmin = np.max(np.array([1, (j - bkg_rad / dy)], dtype=int))
                jmax = np.min(np.array([ny, (j + bkg_rad / dy)], dtype=int))

                n = 0
                sum_ze = 0

                # Calculate the mean background reflectivity for the current
                # grid point, which will be used to determine the convective
                # radius and the required peakedness.

                for r in range(imin, imax):
                    for m in range(jmin, jmax):
                        if not np.isnan(refl[m, r]):
                            rad = np.sqrt((x[r] - x[i]) ** 2 + (y[m] - y[j]) ** 2)

                            # The mean background reflectivity will first be
                            # computed in linear units, i.e. mm^6/m^3, then
                            # converted to decibel units.
                            if rad <= bkg_rad:
                                n += 1
                                sum_ze += 10.0 ** (refl[m, r] / 10.0)

                if n == 0:
                    ze_bkg = np.inf
                else:
                    ze_bkg = 10.0 * np.log10(sum_ze / n)

                # Now get the corresponding convective radius knowing the mean
                # background reflectivity.
                conv_rad = convective_radius(ze_bkg, area_relation)

                # Now we want to investigate the points surrounding the current
                # grid point that are within the convective radius, and whether
                # they too are convective, stratiform or undefined.

                # Get stencil of x and y grid points within the convective
                # radius.
                lmin = np.max(np.array([1, int(i - conv_rad / dx)], dtype=int))
                lmax = np.min(np.array([nx, int(i + conv_rad / dx)], dtype=int))
                mmin = np.max(np.array([1, int(j - conv_rad / dy)], dtype=int))
                mmax = np.min(np.array([ny, int(j + conv_rad / dy)], dtype=int))

                if use_intense and (refl[j, i] >= intense):
                    sclass[j, i] = 2

                    for r in range(lmin, lmax):
                        for m in range(mmin, mmax):
                            if not np.isnan(refl[m, r]):
                                rad = np.sqrt((x[r] - x[i]) ** 2 + (y[m] - y[j]) ** 2)

                                if rad <= conv_rad:
                                    sclass[m, r] = 2

                else:
                    peak = peakedness(ze_bkg, peak_relation)

                    if refl[j, i] - ze_bkg >= peak:
                        sclass[j, i] = 2

                        for r in range(imin, imax):
                            for m in range(jmin, jmax):
                                if not np.isnan(refl[m, r]):
                                    rad = np.sqrt(
                                        (x[r] - x[i]) ** 2 + (y[m] - y[j]) ** 2
                                    )

                                    if rad <= conv_rad:
                                        sclass[m, r] = 2

                    else:
                        # If by now the current grid point has not been
                        # classified as convective by either the intensity
                        # criteria or the peakedness criteria, then it must be
                        # stratiform.
                        sclass[j, i] = 1

    return sclass


def steiner_class_buff(
    ze,
    x,
    y,
    z,
    dx,
    dy,
    bkg_rad,
    work_level,
    intense,
    peak_relation,
    area_relation,
    use_intense,
):
    zslice = np.argmin(np.abs(z - work_level))
    refl = ze[zslice, :, :]

    area_rel = {"small": 0, "medium": 1, "large": 2, "sgp": 3}
    peak_rel = {"default": 0, "sgp": 1}

    sclass = _steiner_conv_strat(
        refl,
        x,
        y,
        dx,
        dy,
        intense=intense,
        peak_relation=peak_rel[peak_relation],
        area_relation=area_rel[area_relation],
        bkg_rad=11000,
        use_intense=True,
    )

    return sclass


def _feature_detection(
    field,
    dx,
    dy,
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
):
    """
    These functions are used to detect features in fields based on how distinct they are from the background
    average. Methodology described in Tomkins et al. (2023), originally based on convective-stratiform algorithms
    developed by Steiner et al. (1995), Yuter and Houze (1997), and Yuter et al. (2005) Grid points are
    classified as follows,

    nosfcecho = No Surface Echo/ Undefined
    bkgd_val = Background echo (e.g. Stratiform)
    feat_val = Feature (e.g. Convective)
    weakecho = Weak Echo

    field : array
        array of values to find features
    x, y : array
        x and y coordinates of field array, respectively
    dx, dy : float
        The x- and y-dimension resolutions in meters, respectively.
    always_core_thres : float, optional
        Threshold for points that are always features. All values above the threshold are classified as features.
    bkg_rad_km : float, optional
        Radius to compute background field in kilometers. Default is 11 km. Recommended to be at least 3 x grid spacing
    use_cosine : bool, optional
        Boolean used to determine if a cosine scheme (see Yuter and Houze (1997)) should be used for identifying
        cores (True) or if a simpler scalar scheme (False) should be used.
    max_diff : float, optional
        Maximum difference between background average and grid value in order to be classified as features.
        "a" value in Eqn. B1 in Yuter and Houze (1997)
    zero_diff_cos_val : float, optional
        Value where difference between background average and grid value is zero in the cosine function
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
        True if using dBZ values that need to be converted to linear Z before averaging. False for other types of values
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

    Returns
    -------
    field_bkg : array
        Array of background values
    core_array : array
        Array of initial cores (identified features without radii applied)
    feature_array : array
        Array of feature detection with radii applied
    """

    # Set up mask arrays for background average and prepare for mask arrays
    # calculate maximum diameter from max. radius (input)
    max_diameter = int(np.floor((max_rad_km / (dx / 1000)) * 2))
    # if diameter is even, make odd
    if max_diameter % 2 == 0:
        max_diameter = max_diameter + 1
    # find center point
    center_mask_x = int(np.floor(max_diameter / 2))

    # prepare background mask array for computing background average
    # calculate number of pixels for background array given requested background radius and dx
    bkg_diameter_pix = int(np.floor((bkg_rad_km / (dx / 1000)) * 2))
    # set diameter to odd if even
    if bkg_diameter_pix % 2 == 0:
        bkg_diameter_pix = bkg_diameter_pix + 1
    # find center point
    bkg_center = int(np.floor(bkg_diameter_pix / 2))
    # create background array
    bkg_mask_array = np.ones((bkg_diameter_pix, bkg_diameter_pix), dtype=float)
    # mask outside circular region
    bkg_mask_array = create_radial_mask(
        bkg_mask_array,
        min_rad_km=0,
        max_rad_km=bkg_rad_km,
        x_pixsize=dx / 1000,
        y_pixsize=dy / 1000,
        center_x=bkg_center,
        center_y=bkg_center,
        circular=True,
    )

    # Feature detection
    # start by making field a masked array
    field = np.ma.masked_invalid(field)
    # Compute background radius
    field_bkg = calc_bkg_intensity(field, bkg_mask_array, dB_averaging, calc_thres)

    # Get core array from cosine scheme, or scalar scheme
    if use_cosine:
        core_array = core_cos_scheme(
            field, field_bkg, max_diff, zero_diff_cos_val, always_core_thres, core_val
        )
    else:
        core_array = core_scalar_scheme(
            field,
            field_bkg,
            scalar_diff,
            always_core_thres,
            core_val,
            use_addition=use_addition,
        )

    # Assign radii based on background field
    radius_array_km = assign_feature_radius_km(
        field_bkg, val_for_max_rad=val_for_max_rad, max_rad=max_rad_km
    )

    # remove small objects in core array
    if remove_small_objects:
        # calculate minimum pixel size given dx and dy
        min_pix_size = min_km2_size / ((dx / 1000) * (dy / 1000))
        # label connected objects in core array
        cc_labels, _ = scipy.ndimage.label(core_array)
        # mask labels where core array is masked
        cc_labels = np.ma.masked_where(core_array.mask, cc_labels)

        # loop through each unique label
        for lab in np.unique(cc_labels):
            # calculate number of pixels for each label
            size_lab = np.count_nonzero(cc_labels == lab)
            # if number of pixels is less than minimum, then remove core
            if size_lab < min_pix_size:
                core_array[cc_labels == lab] = 0

    # perform binary closing
    if binary_close:
        # binary closing - returns binary array
        close_core = scipy.ndimage.binary_closing(core_array).astype(int)
        # set values to core values
        core_array = close_core * core_val

    # Incorporate radius using binary dilation
    # Create empty array for assignment
    temp_assignment = np.zeros_like(core_array)

    # Loop through radii
    for radius in np.arange(1, max_rad_km + 1):
        # create mask array for radius incorporation
        radius_mask_array = create_radius_mask(
            max_diameter, radius, dx / 1000, dy / 1000, center_mask_x
        )
        # find location of radius
        temp = radius_array_km == radius
        # get cores for given radius
        temp_core = np.ma.masked_where(~temp, core_array)
        # dilate cores
        temp_dilated = scipy.ndimage.binary_dilation(
            temp_core.filled(0), radius_mask_array
        )
        # add to assignment array
        temp_assignment = temp_assignment + temp_dilated

    # add dilated cores to original array
    core_copy = np.ma.copy(core_array)
    core_copy[temp_assignment >= 1] = core_val

    # Now do feature detection
    feature_array = np.zeros_like(field)
    feature_array = classify_feature_array(
        field,
        feature_array,
        core_copy,
        nosfcecho,
        feat_val,
        bkgd_val,
        weakecho,
        core_val,
        min_val_used,
        weak_echo_thres,
    )
    # mask where field is masked
    feature_array = np.ma.masked_where(field.mask, feature_array)

    return field_bkg, core_array, feature_array


# functions


def create_radial_mask(
    mask_array,
    min_rad_km,
    max_rad_km,
    x_pixsize,
    y_pixsize,
    center_x,
    center_y,
    circular=True,
):
    """
    Computes a radial distance mask, everything with distance between minradiuskm
    and maxradiuskm is assigned 1, everything else is assigned 0. This version can
    handle rectangular arrays and pixels as well as square ones.

    Parameters
    ----------
    mask_array : array
        Array to mask
    min_rad_km, max_rad_km : float
        The minimum and maximum radius of the non-masked region in kilometers.
    x_pixsize, y_pixsize : float
        The pixel size in the x- and y-dimension in kilometers, respectively
    center_x, center_y : int
        The center pixel in the x- and y-dimension, respectively
    circular : bool
        True returns circular mask, False returns a rectangular mask.

    Returns
    -------
    mask_array : array
        Rectangular array masked by a radial distance.
    """

    xsize, ysize = mask_array.shape

    for j in np.arange(0, ysize, 1):
        for i in np.arange(0, xsize, 1):
            # compute range to pixel
            if circular:
                x_range_sq = ((center_x - i) * x_pixsize) ** 2
                y_range_sq = ((center_y - j) * y_pixsize) ** 2
                circ_range = np.sqrt(x_range_sq + y_range_sq)
            # if circular is False, use square mask
            else:
                x_range = abs(int(np.floor(center_x - i) * x_pixsize))
                y_range = abs(int(np.floor(center_y - j) * y_pixsize))

                if x_range > y_range:
                    circ_range = x_range
                else:
                    circ_range = y_range
            # if range is within min and max, set to True
            if (circ_range <= max_rad_km) and (circ_range >= min_rad_km):
                mask_array[j, i] = 1
            else:
                mask_array[j, i] = 0

    return mask_array


def calc_bkg_intensity(field, bkg_mask_array, dB_averaging, calc_thres=None):
    """
    Calculate the background of the given field. The footprint used to
    calculate the average for each pixel is given by bkg_mask_array

    Parameters
    ----------
    field : array
        Field array to compute average
    bkg_mask_array : array
        Array of radial points to use for average
    dB_averaging : bool
        If True, converts dBZ to linear Z before averaging
    calc_thres : float
        Minimum percentage of points needed to be considered in background average calculation

    Returns
    -------
    field_bkg : array
        Array of average values
    """

    # if dBaverage is true, convert reflectivity to linear Z
    if dB_averaging:
        field = 10 ** (field / 10)

    # calculate background field with circular footprint
    field_bkg = scipy.ndimage.generic_filter(
        field.filled(np.nan),
        function=np.nanmean,
        mode="constant",
        footprint=bkg_mask_array.astype(bool),
        cval=np.nan,
    )

    # if calc_thres is not none, then calculate the number of points used to calculate average
    if calc_thres is not None:
        # count valid points
        field_count = scipy.ndimage.generic_filter(
            field.filled(0),
            function=np.count_nonzero,
            mode="constant",
            footprint=bkg_mask_array.astype(bool),
            cval=0,
        )
        # find threshold number of points
        val = calc_thres * np.count_nonzero(bkg_mask_array)
        # mask out values
        field_bkg = np.ma.masked_where(field_count < val, field_bkg)

    # mask where original field is invalid
    field_bkg = np.ma.masked_where(field.mask, field_bkg)

    # if dBaveraging is true, convert background reflectivity to dBZ
    if dB_averaging:
        field_bkg = 10 * np.log10(field_bkg)
        # mask where original field is invalid
        field_bkg = np.ma.masked_where(field.mask, field_bkg)

    return field_bkg


def core_cos_scheme(
    field, field_bkg, max_diff, zero_diff_cos_val, always_core_thres, CS_CORE
):
    """
    Function for assigning cores based on a cosine function

    Parameters
    ----------
    field : array
        Field values
    field_bkg : array
        Background average of field values
    max_diff : float
        Maximum difference between field and field_bkg needed for feature detection
    zero_diff_cos_val : float
        Value where the cosine function returns a zero difference
    always_core_thres : float
        All values above this threshold considered to be features
    CS_CORE : int
        Value assigned to features

    Returns
    -------
    core_array : array
        Array of booleans if point is feature (1) or not (0)
    """

    # initialize entire array to not a core
    core_array = np.zeros_like(field)

    # calculate zeDiff for entire array
    zDiff = max_diff * np.cos((np.pi * field_bkg) / (2 * zero_diff_cos_val))
    zDiff[zDiff < 0] = 0  # where difference less than zero, set to zero
    zDiff[field_bkg < 0] = max_diff  # where background less than zero, set to max. diff

    # set core values
    # where field >= always_core_thres and where difference exceeds minimum
    core_elements = np.logical_or(
        (field >= always_core_thres), (field - field_bkg) >= zDiff
    )
    core_elements = core_elements.filled(0)
    core_array[core_elements] = CS_CORE

    # mask by field array
    core_array = np.ma.masked_where(field.mask, core_array)

    return core_array


def core_scalar_scheme(
    field, field_bkg, max_diff, always_core_thres, CORE, use_addition=False
):
    """
    Function for assigning cores based on a scalar difference

    Parameters
    ----------
    field : array
        Field values
    field_bkg : array
        Background average of field values
    max_diff : float
        Maximum difference between field and field_bkg needed for feature detection
    always_core_thres : float
        All values above this threshold considered to be a feature
    CORE : int
        Value assigned to features
    use_addition : bool
        Boolean to determine if scalar should be added (True) or multiplied (False)

    Returns
    -------
    core_array : array
        Array of booleans if point is feature (1) or not (0)
    """

    # initialize entire array to not a core
    core_array = np.zeros_like(field)

    # calculate zDiff for entire array
    # if addition, add difference. Else, multiply difference
    if use_addition:
        zDiff = (max_diff + field_bkg) - field_bkg
    else:
        zDiff = (max_diff * field_bkg) - field_bkg

    zDiff[zDiff < 0] = 0  # where difference less than zero, set to zero
    zDiff[field_bkg < 0] = 0  # where background less than zero, set to zero

    # set core values
    # where field >= always_core_thres and where difference exceeds minimum
    core_elements = np.logical_or(
        (field >= always_core_thres), (field - field_bkg) >= zDiff
    )
    core_elements = core_elements.filled(0)
    core_array[core_elements] = CORE

    # mask by field array
    core_array = np.ma.masked_where(field.mask, core_array)

    return core_array


def create_radius_mask(max_diameter, radius_km, x_spacing, y_spacing, center_mask_x):
    """
    Creates a circular mask based on input diameter

    Parameters
    ----------
    max_diameter : int
        maximum diameter in kilometers
    radius_km : int
        radius in kilometers
    x_spacing, y_spacing : float
        x- and y-dimension pixel size in meters, respectively
    center_mask_x : int
        index of center point

    Returns
    -------
    feature_mask_array : array
        array masked based on distance of diameter to incorporate
    """

    feature_mask_array = np.zeros((max_diameter, max_diameter))
    feature_mask_array = create_radial_mask(
        feature_mask_array,
        0,
        radius_km,
        x_spacing,
        y_spacing,
        center_mask_x,
        center_mask_x,
        True,
    )

    return feature_mask_array


def assign_feature_radius_km(field_bkg, val_for_max_rad, max_rad=5):
    """
    Assigns the radius in kilometers based on the background values

    Parameters
    ----------
    field_bkg : array
        array of background field values
    val_for_max_rad : float
        field value for maximum radius (5 km)
    max_rad : float, optional
        maximum radius in kilometers

    Returns
    -------
    radius_array_km : array
        array of radii based on background values and val for max. radius
    """

    radius_array_km = np.ones_like(field_bkg)

    radius_array_km[field_bkg >= (val_for_max_rad - 15)] = max_rad - 3
    radius_array_km[field_bkg >= (val_for_max_rad - 10)] = max_rad - 2
    radius_array_km[field_bkg >= (val_for_max_rad - 5)] = max_rad - 1
    radius_array_km[field_bkg >= val_for_max_rad] = max_rad

    return radius_array_km


def classify_feature_array(
    field,
    feature_array,
    core_array,
    NOSFCECHO,
    FEAT_VAL,
    BKGD_VAL,
    WEAKECHO,
    CORE,
    MINDBZUSE,
    WEAKECHOTHRES,
):
    """
    Does an initial feature detection

    Parameters
    ----------
    field : array
        Array of field values
    feature_array : array
        Array with feature detection
    core_array : array
        Array with cores
    NOSFCECHO : int
        Value to assign points classified as no surface echo
    FEAT_VAL : int
        Value to assign points classified as features
    BKGD_VAL : int
        Value to assign points classified as background echo
    WEAKECHO : int
        Value to assign points classfied as weak echo
    CORE : int
        Value assigned to cores in core_array
    MINDBZUSE : float
        Minimum dBZ value to consider in classification, all values below this will be set to NOSFCECHO
    WEAKECHOTHRES : float
        dBZ threshold for weak echo classification, all values below this will be set to WEAKECHO

    Returns
    -------
    feature_array : array
        array of classifications
    """

    # assuming order so that each point is only assigned one time, no overlapping assignment
    # initially, assign every point to stratiform
    feature_array[:] = BKGD_VAL
    # where field is masked, set to no surface echo
    feature_array[field.mask] = NOSFCECHO
    # assign cores to FEAT
    feature_array[core_array == CORE] = FEAT_VAL
    # assign field values less than weakechothres to weak echo
    feature_array[field < WEAKECHOTHRES] = WEAKECHO
    # assign field values less than minimum to no surface echo
    feature_array[field < MINDBZUSE] = NOSFCECHO

    return feature_array
