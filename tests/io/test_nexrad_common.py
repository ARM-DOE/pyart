"""Unit Tests for Py-ART's io/nexrad_common.py module."""

import pytest

from pyart.io import nexrad_common


def test_get_nexrad_location_known_station():
    lat, lon, elev = nexrad_common.get_nexrad_location("KTLX")
    assert lat == pytest.approx(35.33306, abs=1e-3)
    assert lon == pytest.approx(-97.2775, abs=1e-3)
    # bundled table stores elevation in feet (1213 ft); function must
    # return meters
    assert elev == pytest.approx(1213 * 0.3048, abs=1e-6)


def test_get_nexrad_location_repeated_calls_are_stable():
    # Regression test: get_nexrad_location used to mutate the shared
    # NEXRAD_LOCATIONS dict in place (loc["elev"] = loc["elev"] * 0.3048,
    # where loc is a reference into the module-level table, not a copy).
    # A second call for the same station then applied the feet->meters
    # conversion again on top of the already-converted value, silently
    # corrupting the elevation.
    _, _, elev_first = nexrad_common.get_nexrad_location("KTLX")
    _, _, elev_second = nexrad_common.get_nexrad_location("KTLX")
    _, _, elev_third = nexrad_common.get_nexrad_location("KTLX")
    assert elev_first == elev_second == elev_third
    # and the underlying table itself must be untouched (still in feet)
    assert nexrad_common.NEXRAD_LOCATIONS["KTLX"]["elev"] == 1213
