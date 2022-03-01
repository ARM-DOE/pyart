import pyart

from copy import deepcopy


def test_gatemapper():
    # Make fake radar with target
    old_radar = pyart.testing.make_target_radar()
    new_radar = deepcopy(old_radar)
    new_radar.latitude["data"] = old_radar.latitude["data"] + 0.001
    new_radar.longitude["data"] = old_radar.longitude["data"] + 0.001
    gate_mapper = pyart.map.GateMapper(new_radar, old_radar)
    mapped_radar = gate_mapper.mapped_radar(['reflectivity'])
    assert gate_mapper[20, 20] == (25, 27)
    assert mapped_radar.fields[
        'reflectivity']['data'][25, 27] == old_radar.fields[
                'reflectivity']['data'][20, 20]
