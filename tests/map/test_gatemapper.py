from copy import deepcopy

import numpy as np

import pyart


def test_gatemapper():
    # Make fake radar with target
    old_radar = pyart.testing.make_target_radar()
    new_radar = deepcopy(old_radar)
    new_radar.latitude["data"] = old_radar.latitude["data"] + 0.001
    new_radar.longitude["data"] = old_radar.longitude["data"] + 0.001
    old_radar.fields["reflectivity_copy"] = old_radar.fields["reflectivity"]
    gate_mapper = pyart.map.GateMapper(old_radar, new_radar)
    mapped_radar = gate_mapper.mapped_radar(["reflectivity"])

    # Test point outside of 1 min tolerance
    assert gate_mapper[20, 20] == (None, None)
    assert gate_mapper[40, 40] == (40, 33)
    assert (
        mapped_radar.fields["reflectivity"]["data"][40, 33]
        == old_radar.fields["reflectivity"]["data"][40, 40]
    )

    # Check case where source radar has field destination doesn't
    mapped_radar = gate_mapper.mapped_radar(["reflectivity_copy"])
    assert (
        mapped_radar.fields["reflectivity_copy"]["data"][40, 33]
        == old_radar.fields["reflectivity_copy"]["data"][40, 40]
    )


def test_gatemapper_gatefilter():
    # Make fake radar with target
    old_radar = pyart.testing.make_target_radar()
    new_radar = deepcopy(old_radar)
    new_radar.latitude["data"] = old_radar.latitude["data"] + 0.001
    new_radar.longitude["data"] = old_radar.longitude["data"] + 0.001
    gatefilter = pyart.filters.GateFilter(old_radar)
    gatefilter.exclude_below("reflectivity", 40)
    gate_mapper = pyart.map.GateMapper(new_radar, old_radar, gatefilter_src=gatefilter)
    mapped_radar = gate_mapper.mapped_radar(["reflectivity"])
    assert gate_mapper[4, 4] == (26, 11)
    assert np.ma.is_masked(mapped_radar.fields["reflectivity"]["data"][26, 11])
