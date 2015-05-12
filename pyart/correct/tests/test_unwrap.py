""" Unit Tests for Py-ART's correct/unwrap.py module. """

# python test_dealias.py
# to recreate a dealias_plot.png file showing the before and after doppler
# velocities

from __future__ import print_function

import pyart
import numpy as np
from numpy.testing import assert_allclose
from numpy.testing import assert_raises

REF_DATA = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,
            12.5, 13.5, 12.5, 11.5, 10.5, 9.5, 8.5, 7.5, 6.5, 5.5, 4.5, 3.5,
            2.5, 1.5, 0.5]


def test_dealias_unwrap_phase_ray():
    radar, dealias_vel = perform_dealias('ray')
    assert_allclose(dealias_vel['data'][13, :27], REF_DATA)
    assert np.ma.is_masked(dealias_vel['data'][13]) is False


def test_dealias_unwrap_phase_sweep():
    radar, dealias_vel = perform_dealias('sweep')
    assert_allclose(dealias_vel['data'][13, :27], REF_DATA)
    assert np.ma.is_masked(dealias_vel['data'][13]) is False


def test_dealias_unwrap_phase_volume():
    radar, dealias_vel = perform_dealias('volume')
    assert_allclose(dealias_vel['data'][13, :27], REF_DATA)
    assert np.ma.is_masked(dealias_vel['data'][13]) is False


def test_dealias_unwrap_phase_no_gatefilter():
    radar, dealias_vel = perform_dealias('sweep', gatefilter=False)
    assert_allclose(dealias_vel['data'][13, :27], REF_DATA)
    assert np.ma.is_masked(dealias_vel['data'][13]) is False


def test_dealias_unwrap_phase_explicit_gatefilter():
    radar = pyart.testing.make_velocity_aliased_radar()
    gf = pyart.correct.GateFilter(radar)
    radar.fields['velocity']['data'][13, -4:] = [-7.5, 8.5, 0, 0]
    dealias_vel = pyart.correct.dealias_unwrap_phase(radar, gatefilter=gf)
    assert_allclose(dealias_vel['data'][13, :27], REF_DATA)
    assert np.ma.is_masked(dealias_vel['data'][13]) is False


def test_dealias_unwrap_phase_masked_field_sweep():
    radar = pyart.testing.make_velocity_aliased_radar()
    vdata = radar.fields['velocity']['data']
    radar.fields['velocity']['data'] = np.ma.masked_array(vdata)
    radar.fields['velocity']['data'][13, -4:] = [-7.5, 8.5, 0, 0]
    radar.fields['velocity']['data'][180, 25] = np.ma.masked
    dealias_vel = pyart.correct.dealias_unwrap_phase(
        radar, unwrap_unit='volume')
    assert_allclose(dealias_vel['data'][13, :27], REF_DATA)
    assert np.ma.is_masked(dealias_vel['data'][13]) is False
    assert np.ma.is_masked(dealias_vel['data'][180, 25])
    assert not np.ma.is_masked(dealias_vel['data'][180, 24])


def test_dealias_unwrap_phase_masked_field_volume():
    radar = pyart.testing.make_velocity_aliased_radar()
    vdata = radar.fields['velocity']['data']
    radar.fields['velocity']['data'] = np.ma.masked_array(vdata)
    radar.fields['velocity']['data'][13, -4:] = [-7.5, 8.5, 0, 0]
    radar.fields['velocity']['data'][180, 25] = np.ma.masked
    dealias_vel = pyart.correct.dealias_unwrap_phase(
        radar, unwrap_unit='volume')
    assert_allclose(dealias_vel['data'][13, :27], REF_DATA)
    assert np.ma.is_masked(dealias_vel['data'][13]) is False
    assert np.ma.is_masked(dealias_vel['data'][180, 25])
    assert not np.ma.is_masked(dealias_vel['data'][180, 24])


def test_dealias_unwrap_phase_masked_field_ray():
    radar = pyart.testing.make_velocity_aliased_radar()
    vdata = radar.fields['velocity']['data']
    radar.fields['velocity']['data'] = np.ma.masked_array(vdata)
    radar.fields['velocity']['data'][13, -4:] = [-7.5, 8.5, 0, 0]
    radar.fields['velocity']['data'][180, 25] = np.ma.masked
    dealias_vel = pyart.correct.dealias_unwrap_phase(
        radar, unwrap_unit='ray')
    assert_allclose(dealias_vel['data'][13, :27], REF_DATA)
    assert np.ma.is_masked(dealias_vel['data'][13]) is False
    assert np.ma.is_masked(dealias_vel['data'][180, 25])
    assert not np.ma.is_masked(dealias_vel['data'][180, 24])


def test_dealias_unwrap_phase_rhi_sweep():
    radar = pyart.testing.make_velocity_aliased_rhi_radar()
    radar.fields['velocity']['data'][13, -4:] = [-7.5, 8.5, 0, 0]
    dealias_vel = pyart.correct.dealias_unwrap_phase(radar)
    assert_allclose(dealias_vel['data'][13, :27], REF_DATA)
    assert np.ma.is_masked(dealias_vel['data'][13]) is False


def test_dealias_unwrap_phase_rhi_volume():
    radar = pyart.testing.make_velocity_aliased_rhi_radar()
    radar.fields['velocity']['data'][13, -4:] = [-7.5, 8.5, 0, 0]
    dealias_vel = pyart.correct.dealias_unwrap_phase(
        radar, unwrap_unit='volume')
    assert_allclose(dealias_vel['data'][13, :27], REF_DATA)
    assert np.ma.is_masked(dealias_vel['data'][13]) is False


def test_dealias_unwrap_phase_raises():

    # invalid unwrap_unit
    radar = pyart.testing.make_velocity_aliased_radar()
    assert_raises(ValueError, pyart.correct.dealias_unwrap_phase, radar,
                  unwrap_unit='fuzz')

    # no explicit nyquist
    radar = pyart.testing.make_velocity_aliased_radar()
    radar.instrument_parameters = None
    assert_raises(LookupError, pyart.correct.dealias_unwrap_phase, radar)

    # non-sequential
    radar = pyart.testing.make_velocity_aliased_radar()
    radar.azimuth['data'][10] = 190.
    assert_raises(ValueError, pyart.correct.dealias_unwrap_phase, radar)

    # non-aligned sweeps
    radar = pyart.testing.make_empty_ppi_radar(1, 10, 2)
    radar.fields['velocity'] = {'data': np.zeros((20, 1))}
    assert_raises(ValueError, pyart.correct.dealias_unwrap_phase, radar,
                  nyquist_vel=10, unwrap_unit='volume')

    # non-cubic radar
    radar = pyart.testing.make_empty_ppi_radar(1, 10, 2)
    radar.fields['velocity'] = {'data': np.zeros((20, 1))}
    radar.azimuth['data'][-10:] = range(10)
    radar.sweep_end_ray_index['data'][-1] = 18
    assert_raises(ValueError, pyart.correct.dealias_unwrap_phase, radar,
                  nyquist_vel=10, unwrap_unit='volume')

    # invalid scan type
    radar = pyart.testing.make_velocity_aliased_radar()
    radar.scan_type = 'fuzz'
    assert_raises(ValueError, pyart.correct.dealias_unwrap_phase, radar)


def perform_dealias(unwrap_unit='sweep', **kwargs):
    """ Perform velocity dealiasing on reference data. """
    radar = pyart.testing.make_velocity_aliased_radar()
    # speckling that will not be not be dealiased.
    radar.fields['velocity']['data'][13, -4:] = [-7.5, 8.5, 0, 0]
    dealias_vel = pyart.correct.dealias_unwrap_phase(
        radar, unwrap_unit=unwrap_unit, **kwargs)
    return radar, dealias_vel


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
