""" Unit Tests for Py-ART's graph/radar_display.py module. """
# execute this script to create figure_plot_radar_display_*.png files.

# TODO use matplotlib's @image_comparison decorator to compare to file
# in baseline_images directory. Current this test only determines if files can
# be created, not that they are correct.

import matplotlib.pyplot as plt
import pyart


def test_radar_display_rhi(outfile=None):
    radar = pyart.io.read_netcdf(pyart.testing.NETCDF_RHI_FILE)
    radar.metadata['instrument_name'] = ''
    del radar.fields['reflectivity_horizontal']['valid_min']
    del radar.fields['reflectivity_horizontal']['valid_max']
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_rhi('reflectivity_horizontal', 0, ax=ax,
                     axislabels=('X-axis', 'Y-axis'))
    if outfile:
        fig.savefig(outfile)


def test_radar_display_generate_filename():
    radar = pyart.io.read_netcdf(pyart.testing.NETCDF_PPI_FILE)
    display = pyart.graph.RadarDisplay(radar, shift=(0.1, 0.0))
    filename = display.generate_filename('test', 0)
    assert filename == 'xsapr-sgp_test_00_20110520105416.png'


def test_radar_display_ppi(outfile=None):
    radar = pyart.io.read_netcdf(pyart.testing.NETCDF_PPI_FILE)
    display = pyart.graph.RadarDisplay(radar, shift=(0.1, 0.0))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    display.plot_ppi('reflectivity_horizontal', 0, colorbar_flag=False,
                     title="Fancy PPI Plot", axislabels=('X-axis', 'Y-axis'))
    display.plot_colorbar()
    display.plot_range_rings([10, 20, 30, 40], ax=ax)
    display.plot_labels(['tree'], [(36.68, -97.62)], symbols='k+', ax=ax)
    display.plot_cross_hair(2)
    display.set_limits(ylim=[-50, 50], xlim=[-50, 50])
    if outfile:
        fig.savefig(outfile)


if __name__ == "__main__":
    test_radar_display_rhi('figure_radar_display_rhi.png')
    test_radar_display_ppi('figure_radar_display_ppi.png')
