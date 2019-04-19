import numpy as np


def _steiner_conv_strat(refl, x, y, dx, dy, intense=42, peak_relation=0,
                        area_relation=1, bkg_rad=11000, use_intense=True):
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
                conv_rad = 1000.
            elif (ze_bkg >= 30) & (ze_bkg < 35.):
                conv_rad = 2000.
            elif (ze_bkg >= 35.) & (ze_bkg < 40.):
                conv_rad = 3000.
            elif (ze_bkg >= 40.) & (ze_bkg < 45.):
                conv_rad = 4000.
            else:
                conv_rad = 5000.

        if area_relation == 1:
            if ze_bkg < 25:
                conv_rad = 1000.
            elif (ze_bkg >= 25) & (ze_bkg < 30.):
                conv_rad = 2000.
            elif (ze_bkg >= 30.) & (ze_bkg < 35.):
                conv_rad = 3000.
            elif (ze_bkg >= 35.) & (ze_bkg < 40.):
                conv_rad = 4000.
            else:
                conv_rad = 5000.

        if area_relation == 2:
            if ze_bkg < 20:
                conv_rad = 1000.
            elif (ze_bkg >= 20) & (ze_bkg < 25.):
                conv_rad = 2000.
            elif (ze_bkg >= 25.) & (ze_bkg < 30.):
                conv_rad = 3000.
            elif (ze_bkg >= 30.) & (ze_bkg < 35.):
                conv_rad = 4000.
            else:
                conv_rad = 5000.

        if area_relation == 3:
            if ze_bkg < 40:
                conv_rad = 0.
            elif (ze_bkg >= 40) & (ze_bkg < 45.):
                conv_rad = 1000.
            elif (ze_bkg >= 45.) & (ze_bkg < 50.):
                conv_rad = 2000.
            elif (ze_bkg >= 50.) & (ze_bkg < 55.):
                conv_rad = 6000.
            else:
                conv_rad = 8000.

        return conv_rad

    def peakedness(ze_bkg, peak_relation):
        """
        Given a background reflectivity value, we determine what the necessary
        peakedness (or difference) has to be between a grid point's
        reflectivity and the background reflectivity in order for that grid
        point to be labeled convective.
        """
        if peak_relation == 0:
            if ze_bkg < 0.:
                peak = 10.
            elif (ze_bkg >= 0.) and (ze_bkg < 42.43):
                peak = 10. - ze_bkg ** 2 / 180.
            else:
                peak = 0.

        elif peak_relation == 1:
            if ze_bkg < 0.:
                peak = 14.
            elif (ze_bkg >= 0.) and (ze_bkg < 42.43):
                peak = 14. - ze_bkg ** 2 / 180.
            else:
                peak = 4.

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

                for l in range(imin, imax):
                    for m in range(jmin, jmax):
                        if not np.isnan(refl[m, l]):
                            rad = np.sqrt(
                                (x[l] - x[i]) ** 2 + (y[m] - y[j]) ** 2)

                        # The mean background reflectivity will first be
                        # computed in linear units, i.e. mm^6/m^3, then
                        # converted to decibel units.
                            if rad <= bkg_rad:
                                n += 1
                                sum_ze += 10. ** (refl[m, l] / 10.)

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
                lmin = np.max(
                    np.array([1, int(i - conv_rad / dx)], dtype=int))
                lmax = np.min(
                    np.array([nx, int(i + conv_rad / dx)], dtype=int))
                mmin = np.max(
                    np.array([1, int(j - conv_rad / dy)], dtype=int))
                mmax = np.min(
                    np.array([ny, int(j + conv_rad / dy)], dtype=int))

                if use_intense and (refl[j, i] >= intense):
                    sclass[j, i] = 2

                    for l in range(lmin, lmax):
                        for m in range(mmin, mmax):
                            if not np.isnan(refl[m, l]):
                                rad = np.sqrt(
                                    (x[l] - x[i]) ** 2
                                    + (y[m] - y[j]) ** 2)

                                if rad <= conv_rad:
                                    sclass[m, l] = 2

                else:
                    peak = peakedness(ze_bkg, peak_relation)

                    if refl[j, i] - ze_bkg >= peak:
                        sclass[j, i] = 2

                        for l in range(imin, imax):
                            for m in range(jmin, jmax):
                                if not np.isnan(refl[m, l]):
                                    rad = np.sqrt(
                                        (x[l] - x[i]) ** 2
                                        + (y[m] - y[j]) ** 2)

                                    if rad <= conv_rad:
                                        sclass[m, l] = 2

                    else:
                        # If by now the current grid point has not been
                        # classified as convective by either the intensity
                        # criteria or the peakedness criteria, then it must be
                        # stratiform.
                        sclass[j, i] = 1

    return sclass


def steiner_class_buff(ze, x, y, z, dx, dy, bkg_rad,
                       work_level, intense, peak_relation,
                       area_relation, use_intense):

    zslice = np.argmin(np.abs(z - work_level))
    refl = ze[zslice, :, :]

    area_rel = {"small": 0, "medium": 1, "large": 2, "sgp": 3}
    peak_rel = {"default": 0, "sgp": 1}

    sclass = _steiner_conv_strat(refl, x, y, dx, dy, intense=intense,
                                 peak_relation=peak_rel[peak_relation],
                                 area_relation=area_rel[area_relation],
                                 bkg_rad=11000, use_intense=True)

    return sclass
