""" Unit tests for Py-ART's correct/filters.py module. """

import numpy as np
from numpy.testing import assert_raises

import pyart

radar = pyart.testing.make_empty_ppi_radar(10, 36, 1)

# a simple field
fdata = np.tile(np.arange(10.), 36).reshape(36, 10)
radar.add_field('test_field', {'data': fdata})

# more
fdata2 = np.ma.masked_array(fdata, copy=True)
fdata2[2, 2] = np.ma.masked
fdata2[3, 3] = np.NAN
fdata2[4, 4] = np.PINF
fdata2[5, 5] = np.NINF
radar.add_field('test_field2', {'data': fdata2})


def test_gatefilter_init():
    gfilter = pyart.correct.GateFilter(radar)
    assert np.all(gfilter.gate_excluded == np.False_)


def test_gatefilter_exclude_below():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_below('test_field', 5)
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, -1] is np.False_
    gfilter.exclude_below('test_field', 99)
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, -1] is np.True_


def test_gatefilter_exclude_above():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_above('test_field', 5)
    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[0, -1] is np.True_
    gfilter.exclude_above('test_field', -5)
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, -1] is np.True_


def test_gatefilter_exclude_inside():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_inside('test_field', 2, 5)
    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[0, 3] is np.True_
    assert gfilter.gate_excluded[0, -1] is np.False_
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_inside('test_field', 5, 2)
    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[0, 3] is np.True_
    assert gfilter.gate_excluded[0, -1] is np.False_


def test_gatefilter_exclude_outside():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_outside('test_field', 2, 5)
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, 3] is np.False_
    assert gfilter.gate_excluded[0, -1] is np.True_
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_outside('test_field', 5, 2)
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, 3] is np.False_
    assert gfilter.gate_excluded[0, -1] is np.True_


def test_gatefilter_exclude_equal():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_equal('test_field', 2)
    assert gfilter.gate_excluded[0, 2] is np.True_
    assert gfilter.gate_excluded[0, 3] is np.False_


def test_gatefilter_exclude_not_equal():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_not_equal('test_field', 2)
    assert gfilter.gate_excluded[0, 2] is np.False_
    assert gfilter.gate_excluded[0, 3] is np.True_


def test_gatefilter_exclude_all():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_all()
    assert np.all(gfilter.gate_excluded == np.True_)


def test_gatefilter_exclude_none():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_none()
    assert np.all(gfilter.gate_excluded == np.False_)


def test_gatefilter_exclude_masked():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_masked('test_field2')
    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[2, 2] is np.True_
    assert gfilter.gate_excluded[3, 3] is np.False_
    assert gfilter.gate_excluded[4, 4] is np.False_
    assert gfilter.gate_excluded[5, 5] is np.False_


def test_gatefilter_exclude_invalid():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_invalid('test_field2')
    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[2, 2] is np.True_
    assert gfilter.gate_excluded[3, 3] is np.True_
    assert gfilter.gate_excluded[4, 4] is np.True_
    assert gfilter.gate_excluded[5, 5] is np.True_


def test_gatefilter_ops():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_below('test_field', 0.5, op='or')
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, 1] is np.False_
    assert gfilter.gate_excluded[0, 9] is np.False_
    gfilter.exclude_above('test_field', 8.5, op='or')
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, 1] is np.False_
    assert gfilter.gate_excluded[0, 9] is np.True_

    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_below('test_field', 0.5, op='or')
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, 1] is np.False_
    assert gfilter.gate_excluded[0, 9] is np.False_
    gfilter.exclude_above('test_field', 8.5, op='and')
    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[0, 1] is np.False_
    assert gfilter.gate_excluded[0, 9] is np.False_

    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_below('test_field', 0.5, op='or')
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, 1] is np.False_
    assert gfilter.gate_excluded[0, 9] is np.False_
    gfilter.exclude_above('test_field', 8.5, op='new')
    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[0, 1] is np.False_
    assert gfilter.gate_excluded[0, 9] is np.True_


def test_gatefilter_raises():
    gfilter = pyart.correct.GateFilter(radar)
    assert_raises(ValueError, gfilter.exclude_below, 'test_field', 0.5,
                  op='fuzz')
    assert_raises(ValueError, gfilter.exclude_below, 'test_field', 0.5,
                  exclude_masked='fuzz')


#################
# include tests #
#################


def test_gatefilter_gate_included_attribute():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_below('test_field', 0.5, op='or')
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, 1] is np.False_
    assert gfilter.gate_included[0, 0] is np.False_
    assert gfilter.gate_included[0, 1] is np.True_


def test_gatefilter_include_below():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_below('test_field', 5)
    assert gfilter.gate_included[0, 0] is np.True_
    assert gfilter.gate_included[0, -1] is np.False_
    gfilter.include_below('test_field', 99)
    assert gfilter.gate_included[0, 0] is np.True_
    assert gfilter.gate_included[0, -1] is np.True_


def test_gatefilter_include_above():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_above('test_field', 5)
    assert gfilter.gate_included[0, 0] is np.False_
    assert gfilter.gate_included[0, -1] is np.True_
    gfilter.include_above('test_field', -5)
    assert gfilter.gate_included[0, 0] is np.True_
    assert gfilter.gate_included[0, -1] is np.True_


def test_gatefilter_include_inside():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_inside('test_field', 2, 5)
    assert gfilter.gate_included[0, 0] is np.False_
    assert gfilter.gate_included[0, 3] is np.True_
    assert gfilter.gate_included[0, -1] is np.False_
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_inside('test_field', 5, 2)
    assert gfilter.gate_included[0, 0] is np.False_
    assert gfilter.gate_included[0, 3] is np.True_
    assert gfilter.gate_included[0, -1] is np.False_


def test_gatefilter_include_outside():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_outside('test_field', 2, 5)
    assert gfilter.gate_included[0, 0] is np.True_
    assert gfilter.gate_included[0, 3] is np.False_
    assert gfilter.gate_included[0, -1] is np.True_
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_outside('test_field', 5, 2)
    assert gfilter.gate_included[0, 0] is np.True_
    assert gfilter.gate_included[0, 3] is np.False_
    assert gfilter.gate_included[0, -1] is np.True_


def test_gatefilter_include_equal():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_equal('test_field', 2)
    assert gfilter.gate_included[0, 2] is np.True_
    assert gfilter.gate_included[0, 3] is np.False_


def test_gatefilter_include_not_equal():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_not_equal('test_field', 2)
    assert gfilter.gate_included[0, 2] is np.False_
    assert gfilter.gate_included[0, 3] is np.True_


def test_gatefilter_include_all():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_all()
    assert np.all(gfilter.gate_included == np.True_)


def test_gatefilter_include_none():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_none()
    assert np.all(gfilter.gate_included == np.False_)


def test_gatefilter_include_masked():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_not_masked('test_field2')
    assert gfilter.gate_included[0, 0] is np.True_
    assert gfilter.gate_included[2, 2] is np.False_
    assert gfilter.gate_included[3, 3] is np.True_
    assert gfilter.gate_included[4, 4] is np.True_
    assert gfilter.gate_included[5, 5] is np.True_


def test_gatefilter_include_valid():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_valid('test_field2')
    assert gfilter.gate_included[0, 0] is np.True_
    assert gfilter.gate_included[2, 2] is np.False_
    assert gfilter.gate_included[3, 3] is np.False_
    assert gfilter.gate_included[4, 4] is np.False_
    assert gfilter.gate_included[5, 5] is np.False_
