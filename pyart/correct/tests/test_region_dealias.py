""" Unit Tests for Py-ART's correct/region_dealias.py module. """

# python test_
# to recreate a dealias_plot.png file showing the before and after doppler
# velocities

import warnings

import pyart
import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal

REF_DATA = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,
            12.5, 13.5, 12.5, 11.5, 10.5, 9.5, 8.5, 7.5, 6.5, 5.5, 4.5, 3.5,
            2.5, 1.5, 0.5]


def test_dealias_region_based():
    radar, dealias_vel = perform_dealias()
    assert_allclose(dealias_vel['data'][13, :27], REF_DATA)
    assert np.ma.is_masked(dealias_vel['data'][13]) is False


def test_dealias_interval_limits():
    interval_limits = [-10, 0, 10]
    radar, dealias_vel = perform_dealias(interval_limits=interval_limits)
    assert_allclose(dealias_vel['data'][13, :27], REF_DATA)
    assert np.ma.is_masked(dealias_vel['data'][13]) is False


def test_all_masked():
    radar = pyart.testing.make_velocity_aliased_radar()
    radar.fields['velocity']['data'] = np.ma.array(
        radar.fields['velocity']['data'], mask=True)
    dealias_vel = pyart.correct.dealias_region_based(radar)
    assert np.all(np.ma.getmaskarray(dealias_vel['data']))
    assert 'valid_min' not in dealias_vel
    assert 'valid_max' not in dealias_vel


def test_no_edges():
    # Dealiasing should work when there are no edges shared between
    # regions found
    radar = pyart.testing.make_velocity_aliased_radar()
    radar.fields['velocity']['data'] = np.ma.array(
        radar.fields['velocity']['data'], mask=True)
    radar.fields['velocity']['data'][0, 40] = 0.0
    radar.fields['velocity']['data'][180, 40] = 0.0
    dealias_vel = pyart.correct.dealias_region_based(radar)
    assert_almost_equal(dealias_vel['data'][0, 40], 0.0)
    assert_almost_equal(dealias_vel['data'][180, 40], 0.0)


def test_dealias_reference_velocity_field():
    radar = pyart.testing.make_velocity_aliased_radar()
    ref_data = np.ones_like(radar.fields['velocity']['data']) * 20.
    radar.add_field_like('velocity', 'ref_velocity', ref_data)
    dealias_vel = pyart.correct.dealias_region_based(
        radar, ref_vel_field='ref_velocity')
    offset_ref_data = np.array(REF_DATA) + 20
    assert_allclose(dealias_vel['data'][13, :27], offset_ref_data)


def test_keep_original():
    radar = pyart.testing.make_velocity_aliased_radar()
    radar.fields['velocity']['data'][180, 5] = 88
    gf = pyart.filters.GateFilter(radar)
    gf.exclude_above('velocity', 40)

    dealias_vel = pyart.correct.dealias_region_based(
        radar, gatefilter=gf, keep_original=False)
    assert np.ma.is_masked(dealias_vel['data'][180, 5]) is True

    dealias_vel = pyart.correct.dealias_region_based(
        radar, gatefilter=gf, keep_original=True)
    assert_almost_equal(dealias_vel['data'][180, 5], 88)
    assert np.ma.is_masked(dealias_vel['data'][180, 5]) is False


def test_set_limits():

    radar, dealias_vel = perform_dealias(set_limits=True)
    assert 'valid_min' in dealias_vel
    assert_almost_equal(dealias_vel['valid_min'], -30.0)
    assert 'valid_max' in dealias_vel
    assert_almost_equal(dealias_vel['valid_max'], 30.0)

    radar, dealias_vel = perform_dealias(set_limits=False)
    assert 'valid_min' not in dealias_vel
    assert 'valid_max' not in dealias_vel


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


def test_dealias_outside():
    # unwrap when the velocity field contains gates with values outside of
    # the interval limit, but a warnings should be raised
    radar = pyart.testing.make_velocity_aliased_radar()
    radar.fields['velocity']['data'][250, 25] = 20.1
    radar.fields['velocity']['data'][270, 25] = -20.1
    with warnings.catch_warnings(record=True) as w:
        dealias_vel = pyart.correct.dealias_region_based(radar)
        assert len(w) == 1
    assert_allclose(dealias_vel['data'][13, :27], REF_DATA)
    assert np.ma.is_masked(dealias_vel['data'][13]) is False
    assert_almost_equal(dealias_vel['data'][250, 25], 0.10, 2)
    assert_almost_equal(dealias_vel['data'][270, 25], -0.10, 2)


def test_fillvalue(**kwargs):
    radar = pyart.testing.make_velocity_aliased_radar()
    radar.fields['velocity']['_FillValue'] = 8888.0

    # ensure that at least one gate is excluded and should be masked
    radar.fields['velocity']['data'][100, 25] = 44
    gf = pyart.filters.GateFilter(radar)
    gf.exclude_above('velocity', 40)

    dealias_vel = pyart.correct.dealias_region_based(radar, gatefilter=gf)
    assert isinstance(dealias_vel['data'], np.ma.MaskedArray)
    assert np.ma.is_masked(dealias_vel['data'][100, 25])

    assert '_FillValue' in dealias_vel
    assert dealias_vel['_FillValue'] == 8888.0
    assert dealias_vel['data'].fill_value == 8888.0


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
