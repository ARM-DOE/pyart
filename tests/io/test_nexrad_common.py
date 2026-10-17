"""Unit Tests for Py-ART's io/nexrad_common.py module."""

import pytest

from pyart.io import nexrad_common

KTLX_ELEV_FT = 1213.0


@pytest.fixture
def ktlx_in_feet():
    """Put KTLX back into its as-shipped state and restore it afterwards.

    get_nexrad_location deliberately converts the shared NEXRAD_LOCATIONS
    entry in place, so without this any test that reads the table would
    depend on whether something else looked the station up first.
    """
    entry = nexrad_common.NEXRAD_LOCATIONS["KTLX"]
    saved = dict(entry)
    if entry.get("units") == "m":
        entry["elev"] = entry["elev"] / 0.3048
        del entry["units"]
    yield entry
    entry.clear()
    entry.update(saved)


def test_get_nexrad_location_known_station(ktlx_in_feet):
    lat, lon, elev = nexrad_common.get_nexrad_location("KTLX")
    assert lat == pytest.approx(35.33306, abs=1e-3)
    assert lon == pytest.approx(-97.2775, abs=1e-3)
    # table ships feet, the function returns meters
    assert elev == pytest.approx(KTLX_ELEV_FT * 0.3048, abs=1e-6)


def test_get_nexrad_location_repeated_calls_are_stable(ktlx_in_feet):
    # Regression test: the feet->meters conversion used to be re-applied on
    # every call, so a second lookup for the same station came back wrong
    # (369.72 m, then 112.69 m).
    _, _, first = nexrad_common.get_nexrad_location("KTLX")
    _, _, second = nexrad_common.get_nexrad_location("KTLX")
    _, _, third = nexrad_common.get_nexrad_location("KTLX")
    assert first == second == third


def test_get_nexrad_location_tags_table_units(ktlx_in_feet):
    # Converting the shared table is the existing behaviour and is kept on
    # purpose; the units tag is what makes it idempotent.
    entry = ktlx_in_feet
    assert "units" not in entry
    nexrad_common.get_nexrad_location("KTLX")
    assert entry["units"] == "m"
    assert entry["elev"] == pytest.approx(KTLX_ELEV_FT * 0.3048, abs=1e-6)
    nexrad_common.get_nexrad_location("KTLX")
    assert entry["elev"] == pytest.approx(KTLX_ELEV_FT * 0.3048, abs=1e-6)


def test_get_nexrad_location_is_case_insensitive(ktlx_in_feet):
    lower = nexrad_common.get_nexrad_location("ktlx")
    upper = nexrad_common.get_nexrad_location("KTLX")
    assert lower == upper


def test_get_nexrad_location_unknown_station_raises(ktlx_in_feet):
    with pytest.raises(KeyError):
        nexrad_common.get_nexrad_location("XXXX")
