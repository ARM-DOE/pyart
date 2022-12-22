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


def _revised_conv_strat(
    refl,
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
):
    """
    We perform the Yuter and Houze (1997) algorithm for echo classification
    using only the reflectivity field in order to classify each grid point
    as either convective, stratiform or undefined. Grid points are
    classified as follows,

    nosfcecho = No Surface Echo/ Undefined
    sf = Stratiform
    conv = Convective
    weakecho = Weak Echo

    refl : array
        array of reflectivity values
    x, y : array
        x and y coordinates of reflectivity array, respectively
    dx, dy : float
        The x- and y-dimension resolutions in meters, respectively.
    always_core_thres : float, optional
        Threshold for points that are always convective. All values above the threshold are classifed as convective
    bkg_rad_km : float, optional
        Radius to compute background reflectivity in kilometers. Default is 11 km. Recommended to be at least 3 x
        grid spacing
    use_cosine : bool, optional
        Boolean used to determine if cosine scheme should be used for identifying convective cores (True) or a scalar
        scheme (False)
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
    dB_averaging : bool, optional
        True if using dBZ values that need to be converted to linear Z before averaging. False for other types of values
    remove_small_objects : bool, optional
        Determines if small objects should be removed from convective core array. Default is True.
    min_km2_size : float, optional
        Minimum size of convective cores to be considered. Cores less than this size will be removed. Default is 10
        km^2.
    val_for_max_conv_rad : float, optional
        dBZ for maximum convective radius. Convective cores with values above this will have the maximum convective
        radius
    max_conv_rad_km : float, optional
        Maximum radius around convective cores to classify as convective. Default is 5 km.
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

    Returns
    -------
    refl_bkg : array
        Array of background values
    conv_core_array : array
        Array of initial convective cores (identified convective elements without convective radii applied)
    conv_strat_array : array
        Array of convective stratiform classifcation with convective radii applied
    """

    # Set up mask arrays for background average and
    # prepare for convective mask arrays
    # calculate maximum convective diameter from max. convective radius (input)
    max_conv_diameter = int(np.floor((max_conv_rad_km / (dx / 1000)) * 2))
    # if diameter is even, make odd
    if max_conv_diameter % 2 == 0:
        max_conv_diameter = max_conv_diameter + 1
    # find center point
    center_conv_mask_x = int(np.floor(max_conv_diameter / 2))

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

    # Convective stratiform detection
    # start by making reflectivity a masked array
    refl = np.ma.masked_invalid(refl)
    # Compute background radius
    refl_bkg = calc_bkg_intensity(refl, bkg_mask_array, dB_averaging, calc_thres)
    # mask reflectivity field
    refl = np.ma.masked_where(refl_bkg.mask, refl)

    # Get convective core array from cosine scheme, or scalar scheme
    if use_cosine:
        conv_core_array = convcore_cos_scheme(
            refl, refl_bkg, max_diff, zero_diff_cos_val, always_core_thres, cs_core
        )
    else:
        conv_core_array = convcore_scalar_scheme(
            refl,
            refl_bkg,
            scalar_diff,
            always_core_thres,
            cs_core,
            use_addition=use_addition,
        )

    # Assign convective radii based on background reflectivity
    conv_radius_km = assign_conv_radius_km(
        refl_bkg,
        val_for_max_conv_rad=val_for_max_conv_rad,
        max_conv_rad=max_conv_rad_km,
    )

    # remove small objects in convective core array
    if remove_small_objects:
        # calculate minimum pixel size given dx and dy
        min_pix_size = min_km2_size / ((dx / 1000) * (dy / 1000))
        # label connected objects in convective core array
        cc_labels, _ = scipy.ndimage.label(conv_core_array)
        # mask labels where convective core array is masked
        cc_labels = np.ma.masked_where(conv_core_array.mask, cc_labels)

        # loop through each unique label
        for lab in np.unique(cc_labels):
            # calculate number of pixels for each label
            size_lab = np.count_nonzero(cc_labels == lab)
            # if number of pixels is less than minimum, then remove core
            if size_lab < min_pix_size:
                conv_core_array[cc_labels == lab] = 0

    # Incorporate convective radius using binary dilation
    # Create empty array for assignment
    temp_assignment = np.zeros_like(conv_core_array)

    # Loop through radii
    for radius in np.arange(1, max_conv_rad_km + 1):
        # create mask array for radius incorporation
        conv_mask_array = create_conv_radius_mask(
            max_conv_diameter, radius, dx / 1000, dy / 1000, center_conv_mask_x
        )
        # find location of radius
        temp = conv_radius_km == radius
        # get cores for given radius
        temp_core = np.ma.masked_where(~temp, conv_core_array)
        # dilate cores
        temp_dilated = scipy.ndimage.binary_dilation(
            temp_core.filled(0), conv_mask_array
        )
        # add to assignment array
        temp_assignment = temp_assignment + temp_dilated

    # add dilated cores to original array
    conv_core_copy = np.ma.copy(conv_core_array)
    conv_core_copy[temp_assignment >= 1] = cs_core

    # Now do convective stratiform classification
    conv_strat_array = np.zeros_like(refl)
    conv_strat_array = classify_conv_strat_array(
        refl,
        conv_strat_array,
        conv_core_copy,
        nosfcecho,
        conv,
        sf,
        weakecho,
        cs_core,
        min_dBZ_used,
        weak_echo_thres,
    )
    # mask where reflectivity is masked
    conv_strat_array = np.ma.masked_where(refl.mask, conv_strat_array)

    return refl_bkg, conv_core_array, conv_strat_array


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


def calc_bkg_intensity(refl, bkg_mask_array, dB_averaging, calc_thres=None):
    """
    Calculate the background of the given refl array. The footprint used to
    calculate the average for each pixel is given by bkg_mask_array

    Parameters
    ----------
    refl : array
        Reflectivity array to compute average
    bkg_mask_array : array
        Array of radial points to use for average
    dB_averaging : bool
        If True, converts dBZ to linear Z before averaging
    calc_thres : float
        Minimum percentage of points needed to be considered in background average calculation

    Returns
    -------
    refl_bkg : array
        Array of average values
    """

    # if dBaverage is true, convert reflectivity to linear Z
    if dB_averaging:
        refl = 10 ** (refl / 10)

    # calculate background reflectivity with circular footprint
    refl_bkg = scipy.ndimage.generic_filter(
        refl.filled(np.nan),
        function=np.nanmean,
        mode="constant",
        footprint=bkg_mask_array.astype(bool),
        cval=np.nan,
    )

    # if calc_thres is not none, then calculate the number of points used to calculate average
    if calc_thres is not None:
        # count valid points
        refl_count = scipy.ndimage.generic_filter(
            refl.filled(0),
            function=np.count_nonzero,
            mode="constant",
            footprint=bkg_mask_array.astype(bool),
            cval=0,
        )
        # find threshold number of points
        val = calc_thres * np.count_nonzero(bkg_mask_array)
        # mask out values
        refl_bkg = np.ma.masked_where(refl_count < val, refl_bkg)

    # mask where original reflectivity is invalid
    refl_bkg = np.ma.masked_where(refl.mask, refl_bkg)

    # if dBaveraging is true, convert background reflectivity to dBZ
    if dB_averaging:
        refl_bkg = 10 * np.log10(refl_bkg)
        # mask where original reflectivity is invalid
        refl_bkg = np.ma.masked_where(refl.mask, refl_bkg)

    return refl_bkg


def convcore_cos_scheme(
    refl, refl_bkg, max_diff, zero_diff_cos_val, always_core_thres, CS_CORE
):
    """
    Function for assigning convective cores based on a cosine function

    Parameters
    ----------
    refl : array
        Reflectivity values
    refl_bkg : array
        Background average of reflectivity values
    max_diff : float
        Maximum difference between refl and refl_bkg needed for convective classification
    zero_diff_cos_val : float
        Value where the cosine function returns a zero difference
    always_core_thres : float
        All values above this threshold considered to be convective
    CS_CORE : int
        Value assigned to convective pixels

    Returns
    -------
    conv_core_array : array
        Array of booleans if point is convective (1) or not (0)
    """

    # initialize entire array to not a convective core
    conv_core_array = np.zeros_like(refl)

    # calculate zeDiff for entire array
    zDiff = max_diff * np.cos((np.pi * refl_bkg) / (2 * zero_diff_cos_val))
    zDiff[zDiff < 0] = 0  # where difference less than zero, set to zero
    zDiff[refl_bkg < 0] = max_diff  # where background less than zero, set to max. diff

    # set values
    conv_core_array[
        refl >= always_core_thres
    ] = CS_CORE  # where Z is greater than always_core_thres, set to core
    conv_core_array[
        (refl - refl_bkg) >= zDiff
    ] = CS_CORE  # where difference exceeeds minimum, set to core

    return conv_core_array


def convcore_scalar_scheme(
    refl, refl_bkg, max_diff, always_core_thres, CS_CORE, use_addition=False
):
    """
    Function for assigning convective cores based on a scalar difference

    Parameters
    ----------
    refl : array
        Reflectivity values
    refl_bkg : array
        Background average of reflectivity values
    max_diff : float
        Maximum difference between refl and refl_bkg needed for convective classification
    always_core_thres : float
        All values above this threshold considered to be convective
    CS_CORE : int
        Value assigned to convective pixels
    use_addition : bool
        Boolean to determine if scalar should be added (True) or multiplied (False)

    Returns
    -------
    conv_core_array : array
        Array of booleans if point is convective (1) or not (0)
    """

    # initialize entire array to not a convective core
    conv_core_array = np.zeros_like(refl)

    # calculate zDiff for entire array
    # if addition, add difference. Else, multiply difference
    if use_addition:
        zDiff = max_diff + refl_bkg
    else:
        zDiff = max_diff * refl_bkg

    zDiff[zDiff < 0] = 0  # where difference less than zero, set to zero
    zDiff[refl_bkg < 0] = 0  # where background less than zero, set to zero

    # set values
    conv_core_array[
        refl >= always_core_thres
    ] = CS_CORE  # where Z is greater than always_core_thres, set to core
    conv_core_array[
        refl >= zDiff
    ] = CS_CORE  # where difference exceeeds minimum, set to core

    return conv_core_array


def create_conv_radius_mask(
    max_conv_diameter, radius_km, x_spacing, y_spacing, center_conv_mask_x
):
    """
    Does and initial convective stratiform classification

    Parameters
    ----------
    max_conv_diameter : int
        maximum convective diameter in kilometers
    radius_km : int
        convective radius in kilometers
    x_spacing, y_spacing : float
        x- and y-dimension pixel size in meters, respectively
    center_conv_mask_x : int
        index of center point

    Returns
    -------
    conv_mask_array : array
        array masked based on distance of convective diameter
    """

    conv_mask_array = np.zeros((max_conv_diameter, max_conv_diameter))
    conv_mask_array = create_radial_mask(
        conv_mask_array,
        0,
        radius_km,
        x_spacing,
        y_spacing,
        center_conv_mask_x,
        center_conv_mask_x,
        True,
    )

    return conv_mask_array


def assign_conv_radius_km(refl_bkg, val_for_max_conv_rad, max_conv_rad=5):
    """
    Assigns the convective radius in kilometers based on the background values

    Parameters
    ----------
    refl_bkg : array
        array of background reflectivity values
    val_for_max_conv_rad : float
        reflectivity value for maximum convective radius (5 km)
    max_conv_rad : float, optional
        maximum convective radius in kilometers

    Returns
    -------
    convRadiuskm : array
        array of convective radii based on background values and val for max. conv radius
    """

    convRadiuskm = np.ones_like(refl_bkg)

    convRadiuskm[refl_bkg >= (val_for_max_conv_rad - 15)] = max_conv_rad - 3
    convRadiuskm[refl_bkg >= (val_for_max_conv_rad - 10)] = max_conv_rad - 2
    convRadiuskm[refl_bkg >= (val_for_max_conv_rad - 5)] = max_conv_rad - 1
    convRadiuskm[refl_bkg >= val_for_max_conv_rad] = max_conv_rad

    return convRadiuskm


def classify_conv_strat_array(
    refl,
    conv_strat_array,
    conv_core_array,
    NOSFCECHO,
    CONV,
    SF,
    WEAKECHO,
    CS_CORE,
    MINDBZUSE,
    WEAKECHOTHRES,
):
    """
    Does and initial convective stratiform classification

    Parameters
    ----------
    refl : array
        Array of reflectivity values
    conv_strat_array : array
        Array with convective stratiform classifications
    conv_core_array : array
        Array with convective cores
    NOSFCECHO : int
        Value to assign points classified as no surface echo
    CONV : int
        Value to assign points classified as convective
    SF : int
        Value to assign points classified as stratiform
    WEAKECHO : int
        Value to assign points classfied as weak echo
    CS_CORE : int
        Value assigned to convective cores in conv_core_array
    MINDBZUSE : float
        Minimum dBZ value to consider in classification, all values below this will be set to NOSFCECHO
    WEAKECHOTHRES : float
        dBZ threshold for weak echo classification, all values below this will be set to WEAKECHO

    Returns
    -------
    conv_strat_array : array
        array of classifications
    """

    # assuming order so that each point is only assigned one time, no overlapping assignment
    # initially, assign every point to stratiform
    conv_strat_array[:] = SF
    # where reflectivity is masked, set to no surface echo
    conv_strat_array[refl.mask] = NOSFCECHO
    # assign convective cores to CONV
    conv_strat_array[conv_core_array == CS_CORE] = CONV
    # assign reflectivity less than weakechothres to weak echo
    conv_strat_array[refl < WEAKECHOTHRES] = WEAKECHO
    # assign reflectivity less than minimum to no surface echo
    conv_strat_array[refl < MINDBZUSE] = NOSFCECHO

    return conv_strat_array
