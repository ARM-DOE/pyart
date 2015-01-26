""" Unit Tests for Py-ART's correct/region_dealias.py module. """

# python test_
# to recreate a dealias_plot.png file showing the before and after doppler
# velocities

from __future__ import print_function

import pyart
import numpy as np
from numpy.testing import assert_allclose

REF_DATA = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,
            12.5, 13.5, 12.5, 11.5, 10.5, 9.5, 8.5, 7.5, 6.5, 5.5, 4.5, 3.5,
            2.5, 1.5, 0.5]


def test_dealias_region_based():
    radar, dealias_vel = perform_dealias()
    assert_allclose(dealias_vel['data'][13, :27], REF_DATA)
    assert np.ma.is_masked(dealias_vel['data'][13]) is False


def test_dealias_region_based_masked_wrap_around():
    radar = pyart.testing.make_velocity_aliased_radar()
    vdata = radar.fields['velocity']['data']
    radar.fields['velocity']['data'] = np.ma.masked_array(vdata)
    radar.fields['velocity']['data'][13, -4:] = [-7.5, 8.5, 0, 0]
    radar.fields['velocity']['data'][0, 25] = np.ma.masked
    radar.fields['velocity']['data'][-1, 30] = np.ma.masked
    radar.fields['velocity']['data'][180, -1] = np.ma.masked
    radar.fields['velocity']['data'][180, 0] = np.ma.masked
    dealias_vel = pyart.correct.dealias_region_based(
        radar, rays_wrap_around=True)
    assert_allclose(dealias_vel['data'][13, :27], REF_DATA)
    assert np.ma.is_masked(dealias_vel['data'][13]) is False
    assert np.ma.is_masked(dealias_vel['data'][0, 25])
    assert not np.ma.is_masked(dealias_vel['data'][0, 24])
    assert np.ma.is_masked(dealias_vel['data'][-1, 30])
    assert not np.ma.is_masked(dealias_vel['data'][-1, 29])
    assert np.ma.is_masked(dealias_vel['data'][180, -1])
    assert not np.ma.is_masked(dealias_vel['data'][180, -2])
    assert np.ma.is_masked(dealias_vel['data'][180, 0])
    assert not np.ma.is_masked(dealias_vel['data'][180, 1])


def test_dealias_region_based_no_wrap_around():
    radar = pyart.testing.make_velocity_aliased_radar()
    vdata = radar.fields['velocity']['data']
    radar.fields['velocity']['data'] = np.ma.masked_array(vdata)
    radar.fields['velocity']['data'][13, -4:] = [-7.5, 8.5, 0, 0]
    radar.fields['velocity']['data'][0, 25] = np.ma.masked
    radar.fields['velocity']['data'][0, 25] = np.ma.masked
    radar.fields['velocity']['data'][-1, 30] = np.ma.masked
    dealias_vel = pyart.correct.dealias_region_based(
        radar, rays_wrap_around=False)
    assert_allclose(dealias_vel['data'][13, :27], REF_DATA)
    assert np.ma.is_masked(dealias_vel['data'][13]) is False
    assert np.ma.is_masked(dealias_vel['data'][0, 25])
    assert not np.ma.is_masked(dealias_vel['data'][0, 24])
    assert np.ma.is_masked(dealias_vel['data'][-1, 30])
    assert not np.ma.is_masked(dealias_vel['data'][-1, 29])


def perform_dealias(**kwargs):
    """ Perform velocity dealiasing on reference data. """
    radar = pyart.testing.make_velocity_aliased_radar()
    # speckling that will not be not be dealiased.
    radar.fields['velocity']['data'][13, -4:] = [-7.5, 8.5, 0, 0]
    dealias_vel = pyart.correct.dealias_region_based(
        radar, **kwargs)
    return radar, dealias_vel


def main():

    radar, dealias_vel = perform_dealias()
    radar.fields['dealiased_velocity'] = dealias_vel

    # print out results
    print("ray 13 velocitites before dealias:")
    print(radar.fields['velocity']['data'][13])
    print("ray 13 velocities after dealias:")
    print(radar.fields['dealiased_velocity']['data'][13])

    # create plot
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=[5, 10])
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    rd = pyart.graph.RadarDisplay(radar)
    rd.plot_ppi('velocity', 0, ax=ax1, colorbar_flag=False,
                title='', vmin=-10, vmax=10)
    rd.plot_ppi('dealiased_velocity', 0, ax=ax2, colorbar_flag=False,
                title='', vmin=-10, vmax=20)
    plt.show()
    fig.savefig('dealias_plot.png')


if __name__ == "__main__":
    main()
