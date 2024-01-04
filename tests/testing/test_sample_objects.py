""" Unit Tests for Py-ART's testing/sample_objects.py module. """

import numpy as np

from pyart.testing.sample_objects import make_gaussian_storm_grid


def test_gaussian_storm_grid_results_correct():
    """
    Test for the make_gaussian_storm_grid function.

    Checks grid shape, limits, field data, and masking.
    These test are focusing on statistical properties of the storm and not on comparing exact storm values.
    """
    grid_len = 32
    min_value = 5
    max_value = 45
    mask_margin = 3

    # Create grid
    gaussian_storm_2d = make_gaussian_storm_grid()

    # Test Shape
    assert gaussian_storm_2d.fields["reflectivity"]["data"].shape == (
        2,
        grid_len,
        grid_len,
    ), "Grid shape mismatch"

    # Test Data
    assert (
        gaussian_storm_2d.fields["reflectivity"]["data"] is not None
    ), "No data in reflectivity field"

    # Test Masking
    mask = gaussian_storm_2d.fields["reflectivity"]["data"].mask
    assert np.all(mask[:, :mask_margin, :]), "Masking at the boundary is incorrect"
    assert np.all(mask[:, -mask_margin:, :]), "Masking at the boundary is incorrect"
    assert np.all(mask[:, :, :mask_margin]), "Masking at the boundary is incorrect"
    assert np.all(mask[:, :, -mask_margin:]), "Masking at the boundary is incorrect"

    storm_data = gaussian_storm_2d.fields["reflectivity"]["data"]

    # Test for Max and Min
    assert np.isclose(
        np.max(storm_data), max_value
    ), "Maximum value does not match expected"
    assert np.isclose(
        np.min(storm_data[~storm_data.mask]), min_value
    ), "Minimum value does not match expected"

    # Test Mean and SD
    expected_mean = 8.666844653650797
    expected_std = 7.863066829145
    assert np.isclose(
        np.mean(storm_data), expected_mean, atol=5
    ), "Mean value out of expected range"
    assert np.isclose(
        np.std(storm_data), expected_std, atol=5
    ), "Standard deviation out of expected range"

    # Test Central Value
    assert storm_data[0, 15, 15] == max_value, "Maximum value is not at the center"
