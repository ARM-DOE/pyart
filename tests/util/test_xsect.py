""" Unit tests for the xsect.py module. """

import pyart

from numpy.testing import assert_almost_equal


def test_cross_section_ppi():
    radar = pyart.testing.make_target_radar()
    xsect = pyart.util.cross_section_ppi(radar, [45, 60])

    assert xsect.nsweeps == 2
    assert xsect.nrays == 2
    assert len(xsect.time['data']) == 2
    assert_almost_equal(xsect.azimuth['data'][0], 45)
    assert_almost_equal(xsect.azimuth['data'][1], 60)
    assert_almost_equal(xsect.elevation['data'][0], 0.75, 2)
    assert_almost_equal(xsect.elevation['data'][1], 0.75, 2)
    assert xsect.fields['reflectivity']['data'].shape == (2, 50)
    assert len(xsect.sweep_number['data']) == 2
    assert len(xsect.sweep_mode['data']) == 2
    assert len(xsect.fixed_angle['data']) == 2
    assert xsect.scan_type == 'rhi'
    assert_almost_equal(xsect.sweep_start_ray_index['data'][0], 0)
    assert_almost_equal(xsect.sweep_start_ray_index['data'][1], 1)
    assert_almost_equal(xsect.sweep_end_ray_index['data'][0], 0)
    assert_almost_equal(xsect.sweep_end_ray_index['data'][1], 1)
