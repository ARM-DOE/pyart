import numpy as np
from open_radar_data import DATASETS

import pyart
from pyart.retrieve import create_cappi


def test_create_cappi():
    # Load radar data
    file = DATASETS.fetch("RAW_NA_000_125_20080411190016")
    radar = pyart.io.read(file)

    # Create CAPPI at 10000 meters for the 'reflectivity' field
    cappi = create_cappi(radar, fields=["reflectivity"], height=10000)

    # Retrieve the 'reflectivity' field from the generated CAPPI
    reflectivity_cappi = cappi.fields["reflectivity"]

    # Test 1: Check the shape of the reflectivity CAPPI data
    expected_shape = (360, 992)  # As per the sample data provided
    assert (
        reflectivity_cappi["data"].shape == expected_shape
    ), "Shape mismatch in CAPPI data"

    # Test 2: Check the units of the reflectivity CAPPI
    assert (
        reflectivity_cappi["units"] == "dBZ"
    ), "Incorrect units for CAPPI reflectivity"

    # Test 3: Check that the elevation data is correctly set to zero
    assert np.all(
        cappi.elevation["data"] == 0
    ), "Elevation data should be all zeros in CAPPI"

    # Test 4: Verify the fill value
    assert (
        reflectivity_cappi["_FillValue"] == -9999.0
    ), "Incorrect fill value in CAPPI reflectivity"

    # Test 5: Check the long name and comment
    assert (
        reflectivity_cappi["long_name"] == "CAPPI reflectivity at 10000 meters"
    ), "Incorrect long name"
    assert (
        reflectivity_cappi["comment"]
        == "CAPPI reflectivity calculated at a height of 10000 meters"
    ), "Incorrect comment"
