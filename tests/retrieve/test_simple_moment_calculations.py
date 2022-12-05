""" Unit Tests for Py-ART's retrieve/simple_moment_calculation.py module. """

import numpy as np
from numpy.testing import assert_allclose

import pyart

# Setup a test radar object for various retrievals
test_radar = pyart.testing.make_empty_ppi_radar(100, 360, 5)
test_radar.range["data"] = test_radar.range["data"] * 100.0
range_grid = (
    np.meshgrid(test_radar.range["data"], np.ma.ones(test_radar.time["data"].shape))[0]
    + 1.0
)
foo_field = {"data": np.zeros([360 * 5, 100]) + 20.0 * np.log10(range_grid / 1000.0)}
test_radar.add_field("reflectivity", foo_field)


def test_calculate_snr_from_reflectivity():
    """Test the calculate_snr_from_reflectivity function from pyart.retrieve"""
    snr = pyart.retrieve.calculate_snr_from_reflectivity(test_radar, toa=500)
    assert snr["data"].mean() < 1e-6


def test_compute_noisedBZ():
    """Test the compute_noisedBZ function from pyart.retrieve"""
    noise = pyart.retrieve.compute_noisedBZ(
        nrays=test_radar.nrays,
        noisedBZ_val=1.0,
        _range=test_radar.range["data"],
        ref_dist=1.0,
        noise_field="noisedBZ_hh",
    )
    test_radar.add_field_like
    assert noise["data"].min() > 1.0
    assert noise["data"].max() < 45.0


def test_compute_snr():
    noise = pyart.retrieve.compute_noisedBZ(
        nrays=test_radar.nrays,
        noisedBZ_val=1.0,
        _range=test_radar.range["data"],
        ref_dist=1.0,
        noise_field="noisedBZ_hh",
    )
    # Add the noise field to the radar
    test_radar.add_field("noisedBZ_hh", noise)

    # Calculate snr from reflectivity and noise
    snr = pyart.retrieve.compute_snr(
        test_radar, refl_field="reflectivity", noise_field="noisedBZ_hh"
    )
    assert snr["data"].min() > -1.0


def test_calculate_velocity_texture():
    radar = pyart.testing.make_empty_ppi_radar(10, 36, 1)

    # a zero field
    fdata = np.tile(np.arange(10.0), 36).reshape(36, 10)
    fdata3 = np.zeros(fdata.shape)
    radar.add_field("zero_field", {"data": fdata3})

    vel_field = "zero_field"
    texture_field = pyart.retrieve.calculate_velocity_texture(
        radar, vel_field, wind_size=2, nyq=10
    )
    assert np.all(texture_field["data"] == 0)
    texture_field = pyart.retrieve.calculate_velocity_texture(
        radar, vel_field, wind_size=3, nyq=10
    )
    assert np.all(texture_field["data"] == 0)
    texture_field = pyart.retrieve.calculate_velocity_texture(
        radar, vel_field, wind_size=4, nyq=10
    )
    assert np.all(texture_field["data"] == 0)

    # Test none parameters
    radar2 = pyart.io.read(pyart.testing.NEXRAD_ARCHIVE_MSG31_FILE)
    texture_field = pyart.retrieve.calculate_velocity_texture(
        radar2, vel_field=None, wind_size=4, nyq=None
    )
    assert_allclose(
        texture_field["data"][-1][-5:],
        [0.000363, 0.000363, 0.000363, 0.000363, 0.000363],
        atol=1e-03,
    )
