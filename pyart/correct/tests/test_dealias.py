""" Unit Tests for Py-ART's correct/dealias.py module. """

# python test_dealias.py
# to recreate a dealias_plot.png file showing the before and after doppler
# velocities

import datetime
import warnings

import netCDF4
import pyart
import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal
from numpy.testing.decorators import skipif
from numpy.testing import assert_raises


@skipif(not pyart.correct.dealias._FOURDD_AVAILABLE)
def test_dealias_sounding():
    radar, dealias_vel = perform_dealias()
    assert_allclose(
        dealias_vel['data'][13, :27],
        [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,
         12.5, 13.5, 12.5, 11.5, 10.5, 9.5, 8.5, 7.5, 6.5, 5.5, 4.5, 3.5,
         2.5, 1.5, 0.5])
    assert dealias_vel['data'][13, 46] is np.ma.masked
    assert dealias_vel['data'][13, 47] is np.ma.masked


@skipif(not pyart.correct.dealias._FOURDD_AVAILABLE)
def test_set_limits():

    radar, dealias_vel = perform_dealias(set_limits=True)
    assert 'valid_min' in dealias_vel
    assert_almost_equal(dealias_vel['valid_min'], -30.0)
    assert 'valid_max' in dealias_vel
    assert_almost_equal(dealias_vel['valid_max'], 30.0)

    radar, dealias_vel = perform_dealias(set_limits=False)
    assert 'valid_min' not in dealias_vel
    assert 'valid_max' not in dealias_vel


@skipif(not pyart.correct.dealias._FOURDD_AVAILABLE)
def test_dealias_sounding_keep_original():
    radar, dealias_vel = perform_dealias(True)
    assert_allclose(
        dealias_vel['data'][13, :27],
        [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,
         12.5, 13.5, 12.5, 11.5, 10.5, 9.5, 8.5, 7.5, 6.5, 5.5, 4.5, 3.5,
         2.5, 1.5, 0.5])
    assert dealias_vel['data'][13, 46] == -7.5
    assert dealias_vel['data'][13, 47] == 8.5


def perform_dealias(keep_original=False, **kwargs):
    """ Perform velocity dealiasing on reference data. """
    radar = pyart.testing.make_velocity_aliased_radar()
    # speckling that will not be dealiased.
    radar.fields['velocity']['data'][13, -4:] = [-7.5, 8.5, 0, 0]
    height = np.linspace(150, 250, 10).astype('float32')
    speed = np.ones((10), dtype='float32') * 0.5
    direction = np.ones((10), dtype='float32') * 5.
    profile = pyart.core.HorizontalWindProfile(height, speed, direction)
    dealias_vel = pyart.correct.dealias_fourdd(
        radar, sonde_profile=profile, keep_original=keep_original,
        **kwargs)
    return radar, dealias_vel


@skipif(not pyart.correct.dealias._FOURDD_AVAILABLE)
def test_dealias_last_radar_and_sounding():
    radar = pyart.testing.make_velocity_aliased_radar()
    last_radar = pyart.testing.make_velocity_aliased_radar(False)
    height = np.linspace(150, 250, 10).astype('float32')
    speed = np.ones((10), dtype='float32') * 0.5
    direction = np.ones((10), dtype='float32') * 5.
    profile = pyart.core.HorizontalWindProfile(height, speed, direction)
    dealias_vel = pyart.correct.dealias_fourdd(
        radar, sonde_profile=profile, last_radar=last_radar,
        last_vel_field='velocity')
    assert_allclose(
        dealias_vel['data'][13, :27],
        [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,
         12.5, 13.5, 12.5, 11.5, 10.5, 9.5, 8.5, 7.5, 6.5, 5.5, 4.5, 3.5,
         2.5, 1.5, 0.5])
    return


@skipif(not pyart.correct.dealias._FOURDD_AVAILABLE)
def test_dealias_last_radar_and_sounding():
    radar = pyart.testing.make_velocity_aliased_radar()
    last_radar = pyart.testing.make_velocity_aliased_radar(False)
    dealias_vel = pyart.correct.dealias_fourdd(
        radar, last_radar=last_radar, last_vel_field='velocity')
    assert_allclose(
        dealias_vel['data'][13, :27],
        [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,
         12.5, 13.5, 12.5, 11.5, 10.5, 9.5, 8.5, 7.5, 6.5, 5.5, 4.5, 3.5,
         2.5, 1.5, 0.5])
    return


@skipif(not pyart.correct.dealias._FOURDD_AVAILABLE)
def test_error_raising():

    # ValueError when no sounding or last_radar provided
    radar = pyart.testing.make_velocity_aliased_radar()
    assert_raises(ValueError, pyart.correct.dealias_fourdd, radar)

    return


@skipif(not pyart.correct.dealias._FOURDD_AVAILABLE)
def test_segmentation_fault_error():
    # See GitHub issue #571

    # ValueError when sounding is > 999
    radar = pyart.testing.make_velocity_aliased_radar()
    height = np.arange(2000)
    speed = np.zeros((2000, ))
    direction = np.zeros((2000, ))
    profile = pyart.core.HorizontalWindProfile(height, speed, direction)
    assert_raises(
        ValueError, pyart.correct.dealias_fourdd, radar,
        sonde_profile=profile)
    return


# Remove this test when the sounding_ parameters of dealias_fourdd are
# Deprecated
@skipif(not pyart.correct.dealias._FOURDD_AVAILABLE)
def test_dealias_deprecated_arguments():

    radar = pyart.testing.make_velocity_aliased_radar()
    # speckling that will not be dealiased.
    radar.fields['velocity']['data'][13, -4:] = [-7.5, 8.5, 0, 0]
    height = np.linspace(150, 250, 10).astype('float32')
    speed = np.ones((10), dtype='float32') * 0.5
    direction = np.ones((10), dtype='float32') * 5.
    profile = pyart.core.HorizontalWindProfile(height, speed, direction)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=DeprecationWarning)
        dealias_vel = pyart.correct.dealias_fourdd(
            radar, sounding_heights=height, sounding_wind_speeds=speed,
            sounding_wind_direction=direction)
    assert_allclose(
        dealias_vel['data'][13, :27],
        [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,
         12.5, 13.5, 12.5, 11.5, 10.5, 9.5, 8.5, 7.5, 6.5, 5.5, 4.5, 3.5,
         2.5, 1.5, 0.5])
    assert dealias_vel['data'][13, 46] is np.ma.masked
    assert dealias_vel['data'][13, 47] is np.ma.masked


# Remove this test when the find_time_in_interp_sonde function is Deprecated
@skipif(not pyart.correct.dealias._FOURDD_AVAILABLE)
def test_find_time_in_interp_sounde():
    target = datetime.datetime(2011, 5, 10, 11, 30, 8)
    interp_sounde = netCDF4.Dataset(pyart.testing.INTERP_SOUNDE_FILE)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=DeprecationWarning)
        t = pyart.correct.dealias.find_time_in_interp_sonde(
            interp_sounde, target)
    height, speed, direction = t

    assert height.shape == (316,)
    assert speed.shape == (316,)
    assert direction.shape == (316,)

    assert height.dtype == 'float32'
    assert speed.dtype == 'float32'
    assert direction.dtype == 'float32'

    assert_almost_equal(height[100], 2.32, 2)
    assert_almost_equal(speed[100], 15.54, 2) == 15.54
    assert_almost_equal(direction[100], 231.8, 2) == 231.8


if __name__ == "__main__":

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
    fig.savefig('dealias_plot.png')
