#! /usr/bin/env python
# Example script for plotting six panels of rain data from a radar file
# usage : six_panel_rain_plot filename outdir

import sys
import copy

import matplotlib
from pylab import *
import netCDF4

from pyart.graph import radar_display
from pyart.io import radar, py4dd, py_mdv


if __name__ == "__main__":

    # parse the command line arguements
    filename = sys.argv[1]
    outdir = sys.argv[2]

    print "plotting ", filename, " to ", outdir

    # read in the data
    if ".mdv" in filename:
        my_object = py_mdv.read_mdv(filename, debug=True)
        myradar = radar.Radar(my_object)
    elif ".nc" in filename:
        my_object = netCDF4.Dataset(filename)
        myradar = radar.Radar(my_object)
        my_object.close()
    else:
        py4dd.RSL_radar_verbose_on()
        my_object = py4dd.RSL_anyformat_to_radar(filename)
        myradar = radar.Radar(my_object)

    #calc R
    R = 300.0 * (myradar.fields['specific_attenuation']['data']) ** 0.89
    rainrate = copy.deepcopy(myradar.fields['diff_phase'])
    rainrate['data'] = R
    rainrate['valid_min'] = 0.0
    rainrate['valid_max'] = 400.0
    rainrate['standard_name'] = 'rainfall_rate'
    rainrate['long_name'] = 'rainfall_rate'
    rainrate['least_significant_digit'] = 1
    rainrate['units'] = 'mm/hr'
    myradar.fields.update({'rain_rate_A': rainrate})
    my_display = radar_display.radar_display(myradar)

    f = figure(figsize=[15, 18])
    tilt = 0

    subplot(3, 2, 1)
    my_display.plot_ppi('norm_coherent_power', tilt)
    my_display.append_x()
    my_display.append_y()
    gca().set_title(my_display.generate_title('norm_coherent_power', tilt))
    my_display.add_cb()

    subplot(3, 2, 2)
    my_display.plot_ppi('recalculated_diff_phase', tilt)
    gca().set_title(my_display.generate_title('recalculated_diff_phase', tilt))
    my_display.append_x()
    my_display.add_cb()

    subplot(3, 2, 3)
    my_display.plot_ppi('specific_attenuation', tilt)
    gca().set_title(my_display.generate_title('specific_attenuation', tilt))
    my_display.append_x()
    my_display.add_cb()

    subplot(3, 2, 4)
    my_display.plot_ppi('reflectivity_horizontal', tilt, vmin=-0, vmax=60.0)
    gca().set_title(my_display.generate_title('reflectivity_horizontal', tilt))
    my_display.append_x()
    my_display.add_cb()

    subplot(3, 2, 5)
    my_display.plot_ppi('corrected_reflectivity_horizontal', tilt, vmin=-0,
                        vmax=60.0)
    gca().set_title(my_display.generate_title(
        'corrected_reflectivity_horizontal', tilt))
    my_display.append_x()
    my_display.add_cb()

    subplot(3, 2, 6)
    my_display.plot_ppi('rain_rate_A', tilt, vmax=150, cmap='gist_stern')
    gca().set_title(my_display.generate_title('rain_rate_A', tilt))
    my_display.append_x()
    my_display.add_cb()

    figname = my_display.generate_filename('six_panel', tilt)
    savefig(outdir + '/' + figname.replace(' ', '_'))
    close(f)
