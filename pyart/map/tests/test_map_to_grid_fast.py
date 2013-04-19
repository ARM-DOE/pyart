#!/usr/bin/env python
# script to test map_to_grid function, creates figure_map_to_grid_fast.png

import numpy as np
import matplotlib.pyplot as plt
import pyart


def dms_to_d(dms):
    return dms[0] + (dms[1] + dms[2] / 60.0) / 60.0


def my_qrf(xg, yg, zg):
    hz_dis = np.sqrt(xg ** 2 + yg ** 2)
    hyp = np.sqrt(hz_dis ** 2 + zg**2)+1.0
    theta = np.arccos(hz_dis / hyp) * 180.0 / np.pi
    qrf = (zg / 20.0) + np.sqrt(yg ** 2 + xg ** 2) * \
        np.tan(1.0 * 1.5 * np.pi / 180.0) + 500.0 + \
        zg * 1000.0 * theta ** 2/(45.0 ** 2) / 17000.
    return qrf


def test_map_to_grid_fast():

    # read in the data
    radar = pyart.io.read_netcdf('sgpcsaprsurcmacI7.c0.20110520.110100.nc')

    # mask out last 10 gates along each ray
    radar.fields['reflectivity_horizontal']['data'][:, -10:] = np.ma.masked

    # set origin
    cf_lat = dms_to_d([36.0, 36.0, 18.35])
    cf_lon = -1.0 * dms_to_d([97.0, 29.0, 10.69])

    # perform gridding
    grids = pyart.map.map_to_grid(
        (radar,),
        grid_shape=(241, 241, 2),
        grid_limits=((-117000.0, 123000.0), (-100000.0, 140000.0),
                     (1500, 2000)),
        grid_origin=(cf_lat, cf_lon),
        fields=['reflectivity_horizontal', 'copol_coeff'],
        toa=20000.0,
        qrf_func=my_qrf,
        copy_field_data=True)

    # create plot
    refl = grids['reflectivity_horizontal']
    refl = np.ma.masked_equal(refl, -9999.0)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(refl[1], origin='lower')
    fig.savefig('figure_map_to_grid_fast.png')

if __name__ == "__main__":
    test_map_to_grid_fast()
