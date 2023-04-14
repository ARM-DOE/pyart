""" Unit tests for Py-ART's correct/filters.py module. """

import netCDF4
import numpy as np
import pytest

import pyart

radar = pyart.testing.make_empty_ppi_radar(10, 36, 1)

# a simple field
fdata = np.tile(np.arange(10.0), 36).reshape(36, 10)
radar.add_field("test_field", {"data": fdata})

# more
fdata2 = np.ma.masked_array(fdata, copy=True)
fdata2[2, 2] = np.ma.masked
fdata2[3, 3] = np.NAN
fdata2[4, 4] = np.PINF
fdata2[5, 5] = np.NINF
radar.add_field("test_field2", {"data": fdata2})


def test_gatefilter_init():
    gfilter = pyart.correct.GateFilter(radar)
    assert np.all(gfilter.gate_excluded == np.False_)


def test_gatefilter_copy():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter2 = gfilter.copy()
    gfilter.exclude_below("test_field", 5)
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, -1] is np.False_

    assert gfilter2.gate_excluded[0, 0] is np.False_
    assert gfilter2.gate_excluded[0, -1] is np.False_


def test_exclude_transition():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_transition()
    assert np.all(gfilter.gate_excluded == np.False_)

    trans_radar = pyart.testing.make_empty_ppi_radar(12, 5, 1)
    trans_radar.antenna_transition = {"data": np.array([1, 0, 0, 0, 1])}
    gfilter = pyart.correct.GateFilter(trans_radar)
    gfilter.exclude_transition()

    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, 1] is np.True_
    assert gfilter.gate_excluded[-1, 0] is np.True_
    assert gfilter.gate_excluded[-1, 1] is np.True_
    assert gfilter.gate_excluded[1, 0] is np.False_
    assert gfilter.gate_excluded[1, 1] is np.False_


def test_include_not_transition():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_not_transition()
    assert np.all(gfilter.gate_excluded == np.False_)

    trans_radar = pyart.testing.make_empty_ppi_radar(12, 5, 1)
    trans_radar.antenna_transition = {"data": np.array([1, 0, 0, 0, 1])}
    gfilter = pyart.correct.GateFilter(trans_radar, exclude_based=False)
    gfilter.include_not_transition()

    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, 1] is np.True_
    assert gfilter.gate_excluded[-1, 0] is np.True_
    assert gfilter.gate_excluded[-1, 1] is np.True_
    assert gfilter.gate_excluded[1, 0] is np.False_
    assert gfilter.gate_excluded[1, 1] is np.False_


def test_gatefilter_exclude_below():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_below("test_field", 5)
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, -1] is np.False_

    assert gfilter.gate_excluded[0, 5] is np.False_
    gfilter.exclude_below("test_field", 5, inclusive=True)
    assert gfilter.gate_excluded[0, 5] is np.True_

    gfilter.exclude_below("test_field", 99)
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, -1] is np.True_


def test_gatefilter_exclude_above():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_above("test_field", 5)
    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[0, -1] is np.True_

    assert gfilter.gate_excluded[0, 5] is np.False_
    gfilter.exclude_above("test_field", 5, inclusive=True)
    assert gfilter.gate_excluded[0, 5] is np.True_

    gfilter.exclude_above("test_field", -5)
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, -1] is np.True_


def test_gatefilter_exclude_inside():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_inside("test_field", 2, 5, inclusive=False)
    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[0, 3] is np.True_
    assert gfilter.gate_excluded[0, -1] is np.False_

    assert gfilter.gate_excluded[0, 2] is np.False_
    gfilter.exclude_inside("test_field", 2, 5, inclusive=True)
    assert gfilter.gate_excluded[0, 2] is np.True_

    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_inside("test_field", 5, 2)
    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[0, 3] is np.True_
    assert gfilter.gate_excluded[0, -1] is np.False_


def test_gatefilter_exclude_outside():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_outside("test_field", 2, 5)
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, 3] is np.False_
    assert gfilter.gate_excluded[0, -1] is np.True_

    assert gfilter.gate_excluded[0, 2] is np.False_
    gfilter.exclude_outside("test_field", 2, 5, inclusive=True)
    assert gfilter.gate_excluded[0, 2] is np.True_

    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_outside("test_field", 5, 2)
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, 3] is np.False_
    assert gfilter.gate_excluded[0, -1] is np.True_


def test_gatefilter_exclude_equal():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_equal("test_field", 2)
    assert gfilter.gate_excluded[0, 2] is np.True_
    assert gfilter.gate_excluded[0, 3] is np.False_


def test_gatefilter_exclude_not_equal():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_not_equal("test_field", 2)
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
    gfilter.exclude_masked("test_field2")
    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[2, 2] is np.True_
    assert gfilter.gate_excluded[3, 3] is np.False_
    assert gfilter.gate_excluded[4, 4] is np.False_
    assert gfilter.gate_excluded[5, 5] is np.False_


def test_gatefilter_exclude_invalid():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_invalid("test_field2")
    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[2, 2] is np.True_
    assert gfilter.gate_excluded[3, 3] is np.True_
    assert gfilter.gate_excluded[4, 4] is np.True_
    assert gfilter.gate_excluded[5, 5] is np.True_


def test_gatefilter_exclude_gates():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.include_all()
    gates = gfilter.gate_excluded
    gates[2, 0] = np.True_
    gates[2, 2] = np.True_
    gfilter.exclude_gates(gates)
    # exclude when included
    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[2, 0] is np.True_
    # exclude when already excluded
    gfilter.exclude_gates(gates)
    assert gfilter.gate_excluded[2, 0] is np.True_
    gates[2, 0] = np.False_
    gates[0, 2] = np.True_
    # exclude with op='and'
    gfilter.exclude_gates(gates, op="and")
    assert gfilter.gate_excluded[2, 0] is np.False_
    assert gfilter.gate_excluded[0, 2] is np.False_
    assert gfilter.gate_excluded[2, 2] is np.True_


def test_gatefilter_ops():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_below("test_field", 0.5, op="or")
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, 1] is np.False_
    assert gfilter.gate_excluded[0, 9] is np.False_
    gfilter.exclude_above("test_field", 8.5, op="or")
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, 1] is np.False_
    assert gfilter.gate_excluded[0, 9] is np.True_

    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_below("test_field", 0.5, op="or")
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, 1] is np.False_
    assert gfilter.gate_excluded[0, 9] is np.False_
    gfilter.exclude_above("test_field", 8.5, op="and")
    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[0, 1] is np.False_
    assert gfilter.gate_excluded[0, 9] is np.False_

    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_below("test_field", 0.5, op="or")
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, 1] is np.False_
    assert gfilter.gate_excluded[0, 9] is np.False_
    gfilter.exclude_above("test_field", 8.5, op="new")
    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[0, 1] is np.False_
    assert gfilter.gate_excluded[0, 9] is np.True_


def test_gatefilter_raises():
    gfilter = pyart.correct.GateFilter(radar)
    pytest.raises(ValueError, gfilter.exclude_below, "test_field", 0.5, op="fuzz")
    pytest.raises(
        ValueError, gfilter.exclude_below, "test_field", 0.5, exclude_masked="fuzz"
    )
    pytest.raises(ValueError, gfilter.exclude_gates, np.arange(2))
    pytest.raises(ValueError, gfilter.include_gates, np.arange(2))


#################
# include tests #
#################


def test_gatefilter_gate_included_attribute():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_below("test_field", 0.5, op="or")
    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, 1] is np.False_
    assert gfilter.gate_included[0, 0] is np.False_
    assert gfilter.gate_included[0, 1] is np.True_


def test_gatefilter_include_below():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_below("test_field", 5)
    assert gfilter.gate_included[0, 0] is np.True_
    assert gfilter.gate_included[0, -1] is np.False_

    assert gfilter.gate_included[0, 5] is np.False_
    gfilter.include_below("test_field", 5, inclusive=True)
    assert gfilter.gate_included[0, 5] is np.True_

    gfilter.include_below("test_field", 99)
    assert gfilter.gate_included[0, 0] is np.True_
    assert gfilter.gate_included[0, -1] is np.True_


def test_gatefilter_include_above():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_above("test_field", 5)
    assert gfilter.gate_included[0, 0] is np.False_
    assert gfilter.gate_included[0, -1] is np.True_

    assert gfilter.gate_included[0, 5] is np.False_
    gfilter.include_above("test_field", 5, inclusive=True)
    assert gfilter.gate_included[0, 5] is np.True_

    gfilter.include_above("test_field", -5)
    assert gfilter.gate_included[0, 0] is np.True_
    assert gfilter.gate_included[0, -1] is np.True_


def test_gatefilter_include_inside():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_inside("test_field", 2, 5, inclusive=False)
    assert gfilter.gate_included[0, 0] is np.False_
    assert gfilter.gate_included[0, 3] is np.True_
    assert gfilter.gate_included[0, -1] is np.False_

    assert gfilter.gate_included[0, 2] is np.False_
    gfilter.include_inside("test_field", 2, 5, inclusive=True)
    assert gfilter.gate_included[0, 2] is np.True_

    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_inside("test_field", 5, 2)
    assert gfilter.gate_included[0, 0] is np.False_
    assert gfilter.gate_included[0, 3] is np.True_
    assert gfilter.gate_included[0, -1] is np.False_


def test_gatefilter_include_outside():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_outside("test_field", 2, 5)
    assert gfilter.gate_included[0, 0] is np.True_
    assert gfilter.gate_included[0, 3] is np.False_
    assert gfilter.gate_included[0, -1] is np.True_

    assert gfilter.gate_included[0, 2] is np.False_
    gfilter.include_outside("test_field", 2, 5, inclusive=True)
    assert gfilter.gate_included[0, 2] is np.True_

    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_outside("test_field", 5, 2)
    assert gfilter.gate_included[0, 0] is np.True_
    assert gfilter.gate_included[0, 3] is np.False_
    assert gfilter.gate_included[0, -1] is np.True_


def test_gatefilter_include_equal():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_equal("test_field", 2)
    assert gfilter.gate_included[0, 2] is np.True_
    assert gfilter.gate_included[0, 3] is np.False_


def test_gatefilter_include_not_equal():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_not_equal("test_field", 2)
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
    gfilter.include_not_masked("test_field2")
    assert gfilter.gate_included[0, 0] is np.True_
    assert gfilter.gate_included[2, 2] is np.False_
    assert gfilter.gate_included[3, 3] is np.True_
    assert gfilter.gate_included[4, 4] is np.True_
    assert gfilter.gate_included[5, 5] is np.True_


def test_gatefilter_include_valid():
    gfilter = pyart.correct.GateFilter(radar, exclude_based=False)
    gfilter.include_valid("test_field2")
    assert gfilter.gate_included[0, 0] is np.True_
    assert gfilter.gate_included[2, 2] is np.False_
    assert gfilter.gate_included[3, 3] is np.False_
    assert gfilter.gate_included[4, 4] is np.False_
    assert gfilter.gate_included[5, 5] is np.False_


def test_gatefilter_include_gates():
    gfilter = pyart.correct.GateFilter(radar)
    gfilter.exclude_all()
    gates = gfilter.gate_included
    gates[2, 0] = np.True_
    gates[2, 2] = np.True_
    gfilter.include_gates(gates)
    # include when excluded
    assert gfilter.gate_included[0, 0] is np.False_
    assert gfilter.gate_included[2, 0] is np.True_
    # include when already included
    gfilter.include_gates(gates)
    assert gfilter.gate_included[2, 0] is np.True_
    gates[2, 0] = np.False_
    gates[0, 2] = np.True_
    # include with op='or'
    gfilter.include_gates(gates, op="or")
    assert gfilter.gate_included[2, 0] is np.False_
    assert gfilter.gate_included[0, 2] is np.False_
    assert gfilter.gate_included[2, 2] is np.True_


radar2 = pyart.io.read(pyart.testing.NEXRAD_ARCHIVE_MSG31_FILE)


def test_moment_based_gate_filter():
    fdata = np.full((7200, 1832), 0.5)
    radar2.add_field("normalized_coherent_power", {"data": fdata})

    gfilter = pyart.correct.moment_based_gate_filter(
        radar2, min_ncp=0.4, min_refl=-33.0, max_refl=-31.0, min_rhv=0.207
    )
    gfilter2 = pyart.correct.moment_based_gate_filter(
        radar2, min_ncp=0.6, min_refl=-31.0, max_refl=-28.0, min_rhv=0.209
    )

    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[0, -1] is np.True_

    assert gfilter2.gate_excluded[0, 0] is np.True_
    assert gfilter2.gate_excluded[0, -1] is np.True_


sonde = netCDF4.Dataset(pyart.testing.SONDE_FILE)
radar3 = pyart.io.read(pyart.testing.NEXRAD_ARCHIVE_MSG31_COMPRESSED_FILE)
z_dict, temp_dict = pyart.retrieve.map_profile_to_gates(
    sonde.variables["tdry"][:], sonde.variables["alt"][:], radar3, toa=16500
)
radar3.add_field("temperature", temp_dict, replace_existing=True)
radar3.add_field("height", z_dict, replace_existing=True)


def test_moment_and_texture_based_gate_filter():
    gfilter = pyart.filters.moment_and_texture_based_gate_filter(
        radar3,
        phi_field="differential_phase",
        max_textphi=20.0,
        max_textrhv=0.3,
        max_textzdr=2.85,
        max_textrefl=8.0,
        min_rhv=0.6,
    )

    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, -1] is np.True_

    gfilter2 = pyart.filters.moment_and_texture_based_gate_filter(
        radar3,
        phi_field="differential_phase",
        max_textphi=280.0,
        max_textrhv=0.9,
        max_textzdr=10,
        max_textrefl=33.0,
        min_rhv=0.2,
    )

    assert gfilter2.gate_excluded[0, 0] is np.False_
    assert gfilter2.gate_excluded[0, -1] is np.True_


def test_temp_based_gate_filter():
    gfilter = pyart.filters.temp_based_gate_filter(radar3, min_temp=17.0)

    assert gfilter.gate_excluded[0, 0] is np.True_
    assert gfilter.gate_excluded[0, -1] is np.True_

    gfilter2 = pyart.filters.temp_based_gate_filter(
        radar3, thickness=None, min_temp=17.0
    )

    assert gfilter2.gate_excluded[0, 27] is np.False_
    assert gfilter2.gate_excluded[0, -1] is np.True_


def test_iso0_based_gate_filter():
    gfilter = pyart.filters.iso0_based_gate_filter(
        radar3, iso0_field="height", max_h_iso0=None, thickness=None, beamwidth=1.0
    )

    assert gfilter.gate_excluded[0, 0] is np.False_
    assert gfilter.gate_excluded[0, -1] is np.False_

    gfilter2 = pyart.filters.iso0_based_gate_filter(
        radar3, iso0_field="height", max_h_iso0=1000.0, thickness=400.0, beamwidth=None
    )

    assert gfilter2.gate_excluded[0, 0] is np.False_
    assert gfilter2.gate_excluded[0, -1] is np.True_
