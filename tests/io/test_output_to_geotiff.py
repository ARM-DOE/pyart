""" Unit Tests for Py-ART's output_to_geotiff.py module. """

import warnings

import numpy as np
import pytest

import pyart

# TODO : inspect the output file to verify their contents, currently only the
# fact that something was written is confirmed


def test__get_rgb_values_nan():
    data = np.arange(10, dtype="float32")
    data[5] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        rarr, barr, garr, aarr = pyart.io.output_to_geotiff._get_rgb_values(
            data, 0, 10, None, "jet", False, 1
        )
    assert np.isnan(rarr[5])
    assert np.isnan(barr[5])
    assert np.isnan(garr[5])


def test_raise_missingoptionaldepedency():
    backup = bool(pyart.io.output_to_geotiff.IMPORT_FLAG)
    pyart.io.output_to_geotiff.IMPORT_FLAG = False
    grid = make_tiny_grid()
    pytest.raises(
        pyart.exceptions.MissingOptionalDependency,
        pyart.io.write_grid_geotiff,
        grid,
        "test.foo",
        "reflectivity",
    )
    pyart.io.output_to_geotiff.IMPORT_FLAG = backup


def make_tiny_grid():
    """Make a tiny grid."""
    grid_shape = (2, 10, 8)
    grid_limits = ((0, 500), (-400000, 400000), (-300000, 300000))
    grid = pyart.testing.make_empty_grid(grid_shape, grid_limits)
    fdata = np.zeros((2, 10, 8), dtype="float32")
    fdata[:, 2:-2, 1:-1] = 10.0
    fdata[:, 3:-3, 2:-2] = 20.0
    fdata[:, 4:-4, 3:-3] = 30.0
    fdata[1] += 5
    rdic = {"data": fdata, "long_name": "reflectivity", "units": "dBz"}
    grid.fields = {"reflectivity": rdic}
    return grid


def make_tiny_grid_with_mask():
    """Make a tiny grid."""
    grid_shape = (2, 10, 8)
    grid_limits = ((0, 500), (-400000, 400000), (-300000, 300000))
    grid = pyart.testing.make_empty_grid(grid_shape, grid_limits)
    fdata = np.zeros((2, 10, 8), dtype="float32")
    fdata[:, 2:-2, 1:-1] = 10.0
    fdata[:, 3:-3, 2:-2] = 20.0
    fdata[:, 4:-4, 3:-3] = 30.0
    fdata[1] += 5
    fdata[fdata == 0] = np.nan
    rdic = {"data": fdata, "long_name": "reflectivity", "units": "dBz"}
    grid.fields = {"reflectivity": rdic}
    return grid


@pytest.mark.skipif(
    not pyart.io.output_to_geotiff.IMPORT_FLAG, reason="GDAL is not installed."
)
def test_write_grid_geotiff_tif_single_channel():
    grid = make_tiny_grid()
    with pyart.testing.InTemporaryDirectory():
        pyart.io.write_grid_geotiff(grid, "test.foo", "reflectivity")
        # check that something was written to the file
        with open("test.foo", "rb") as f:
            assert len(f.read()) > 0


@pytest.mark.skipif(
    not pyart.io.output_to_geotiff.IMPORT_FLAG, reason="GDAL is not installed."
)
def test_write_grid_geotiff_tif_rgb_with_missing():
    grid = make_tiny_grid()
    with pyart.testing.InTemporaryDirectory():
        pyart.io.write_grid_geotiff(grid, "test.foo", "reflectivity", rgb=True)
        # check that something was written to the file
        with open("test.foo", "rb") as f:
            assert len(f.read()) > 0


@pytest.mark.skipif(
    not pyart.io.output_to_geotiff.IMPORT_FLAG, reason="GDAL is not installed."
)
def test_write_grid_geotiff_tif_level():
    grid = make_tiny_grid()
    with pyart.testing.InTemporaryDirectory():
        pyart.io.write_grid_geotiff(grid, "test.foo", "reflectivity", level=1)
        # check that something was written to the file
        with open("test.foo", "rb") as f:
            assert len(f.read()) > 0


@pytest.mark.skipif(
    not pyart.io.output_to_geotiff.IMPORT_FLAG, reason="GDAL is not installed."
)
def test_write_grid_geotiff_tif_warp():
    grid = make_tiny_grid()
    with pyart.testing.InTemporaryDirectory():
        pyart.io.write_grid_geotiff(grid, "test.foo", "reflectivity", warp=True)
        # check that something was written to the file
        with open("test.foo", "rb") as f:
            assert len(f.read()) > 0


@pytest.mark.skipif(
    not pyart.io.output_to_geotiff.IMPORT_FLAG, reason="GDAL is not installed."
)
def test_write_grid_geotiff_tif_missing_suffix():
    grid = make_tiny_grid()
    with pyart.testing.InTemporaryDirectory():
        pyart.io.write_grid_geotiff(grid, "test", "reflectivity")
        # suffix missing, should be added to filename
        # check that something was written to the file
        with open("test.tif", "rb") as f:
            assert len(f.read()) > 0


@pytest.mark.skipif(
    not pyart.io.output_to_geotiff.IMPORT_FLAG, reason="GDAL is not installed."
)
def test_write_grid_geotiff_sld():
    grid = make_tiny_grid()
    with pyart.testing.InTemporaryDirectory():
        pyart.io.write_grid_geotiff(grid, "test.foo", "reflectivity", sld=True)
        # check that something was written to the file
        with open("test.foo", "rb") as f:
            assert len(f.read()) > 0


@pytest.mark.skipif(
    not pyart.io.output_to_geotiff.IMPORT_FLAG, reason="GDAL is not installed."
)
def test_write_grid_geotiff_missing_field():
    grid = make_tiny_grid()
    pytest.raises(KeyError, pyart.io.write_grid_geotiff, grid, "test.foo", "foobar")


@pytest.mark.skipif(
    not pyart.io.output_to_geotiff.IMPORT_FLAG, reason="GDAL is not installed."
)
def test_write_grid_geotiff_transparent_background():
    grid = make_tiny_grid_with_mask()
    with pyart.testing.InTemporaryDirectory():
        pyart.io.write_grid_geotiff(
            grid,
            "transparent.foo",
            "reflectivity",
            rgb=True,
            cmap="HomeyerRainbow",
            vmin=0,
            vmax=40,
            transparent_bg=True,
            opacity=1,
        )
        with open("transparent.foo", "rb") as f:
            assert len(f.read()) > 0


@pytest.mark.skipif(
    not pyart.io.output_to_geotiff.IMPORT_FLAG, reason="GDAL is not installed."
)
def test_write_grid_geotiff_opacity():
    grid = make_tiny_grid_with_mask()
    with pyart.testing.InTemporaryDirectory():
        pyart.io.write_grid_geotiff(
            grid,
            "opacity.foo",
            "reflectivity",
            rgb=True,
            cmap="HomeyerRainbow",
            vmin=0,
            vmax=40,
            transparent_bg=False,
            opacity=0.25,
        )
        with open("opacity.foo", "rb") as f:
            assert len(f.read()) > 0
