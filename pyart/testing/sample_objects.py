"""
Functions for creating sample Radar and Grid objects.

"""

import numpy as np
import scipy

from ..config import get_metadata
from ..core.grid import Grid
from ..core.radar import Radar
from ..exceptions import MissingOptionalDependency
from .sample_files import _EXAMPLE_RAYS_FILE

try:
    import xarray as xr

    from ..core.radar_spectra import RadarSpectra

    _XARRAY_AVAILABLE = True
except ImportError:
    _XARRAY_AVAILABLE = False


def make_empty_ppi_radar(ngates, rays_per_sweep, nsweeps):
    """
    Return an Radar object, representing a PPI scan.

    Parameters
    ----------
    ngates : int
        Number of gates per ray.
    rays_per_sweep : int
        Number of rays in each PPI sweep.
    nsweeps : int
        Number of sweeps.

    Returns
    -------
    radar : Radar
        Radar object with no fields, other parameters are set to default
        values.

    """
    nrays = rays_per_sweep * nsweeps

    time = get_metadata("time")
    _range = get_metadata("range")
    latitude = get_metadata("latitude")
    longitude = get_metadata("longitude")
    altitude = get_metadata("altitude")
    sweep_number = get_metadata("sweep_number")
    sweep_mode = get_metadata("sweep_mode")
    fixed_angle = get_metadata("fixed_angle")
    sweep_start_ray_index = get_metadata("sweep_start_ray_index")
    sweep_end_ray_index = get_metadata("sweep_end_ray_index")
    azimuth = get_metadata("azimuth")
    elevation = get_metadata("elevation")

    fields = {}
    scan_type = "ppi"
    metadata = {"instrument_name": "fake_radar"}

    time["data"] = np.arange(nrays, dtype="float64")
    time["units"] = "seconds since 1989-01-01T00:00:01Z"
    _range["data"] = np.linspace(0, 1000, ngates).astype("float32")

    latitude["data"] = np.array([36.5], dtype="float64")
    longitude["data"] = np.array([-97.5], dtype="float64")
    altitude["data"] = np.array([200], dtype="float64")

    sweep_number["data"] = np.arange(nsweeps, dtype="int32")
    sweep_mode["data"] = np.array(["azimuth_surveillance"] * nsweeps)
    fixed_angle["data"] = np.array([0.75] * nsweeps, dtype="float32")
    sweep_start_ray_index["data"] = np.arange(0, nrays, rays_per_sweep, dtype="int32")
    sweep_end_ray_index["data"] = np.arange(
        rays_per_sweep - 1, nrays, rays_per_sweep, dtype="int32"
    )

    azimuth["data"] = np.arange(nrays, dtype="float32")
    elevation["data"] = np.array([0.75] * nrays, dtype="float32")

    return Radar(
        time,
        _range,
        fields,
        metadata,
        scan_type,
        latitude,
        longitude,
        altitude,
        sweep_number,
        sweep_mode,
        fixed_angle,
        sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth,
        elevation,
        instrument_parameters=None,
    )


def make_empty_rhi_radar(ngates, rays_per_sweep, nsweeps):
    """
    Return an Radar object, representing a RHI scan.

    Parameters
    ----------
    ngates : int
        Number of gates per ray.
    rays_per_sweep : int
        Number of rays in each PPI sweep.
    nsweeps : int
        Number of sweeps.

    Returns
    -------
    radar : Radar
        Radar object with no fields, other parameters are set to default
        values.

    """
    radar = make_empty_ppi_radar(ngates, rays_per_sweep, nsweeps)
    radar.scan_type = "rhi"
    nrays = rays_per_sweep * nsweeps
    radar.sweep_mode["data"] = np.array(["rhi"] * nsweeps)
    radar.elevation["data"] = np.arange(nrays, dtype="float32")
    radar.azimuth["data"] = np.array([0.75] * nrays, dtype="float32")
    return radar


def make_target_radar():
    """
    Return a PPI radar with a target like reflectivity field.
    """
    radar = make_empty_ppi_radar(50, 360, 1)
    fields = {"reflectivity": get_metadata("reflectivity")}
    fdata = np.zeros((360, 50), dtype="float32")
    fdata[:, 0:10] = 0.0
    fdata[:, 10:20] = 10.0
    fdata[:, 20:30] = 20.0
    fdata[:, 30:40] = 30.0
    fdata[:, 40:50] = 40.0
    fields["reflectivity"]["data"] = fdata
    radar.fields = fields
    return radar


def make_velocity_aliased_radar(alias=True):
    """
    Return a PPI radar with a target like reflectivity field.

    Set alias to False to return a de-aliased radar.
    """
    radar = make_empty_ppi_radar(50, 360, 1)
    radar.range["meters_between_gates"] = 1.0
    radar.range["meters_to_center_of_first_gate"] = 1.0
    radar.instrument_parameters = {"nyquist_velocity": {"data": np.array([10.0] * 360)}}

    fields = {
        "reflectivity": get_metadata("reflectivity"),
        "velocity": get_metadata("velocity"),
    }

    # fake reflectivity data, all zero reflectivity
    fdata = np.zeros((360, 50), dtype="float32")
    fields["reflectivity"]["data"] = fdata

    # fake velocity data, all zeros except a wind burst on at ~13 degrees.
    # burst is partially aliased.
    vdata = np.zeros((360 * 1, 50), dtype="float32")

    for i, idx in enumerate(range(13, -1, -1)):
        vdata[i, idx : idx + i + 1] = np.arange(0.5, 0.5 + i * 1.0 + 0.001)
    vdata[:14, 14:27] = vdata[:14, 12::-1]  # left/right flip
    vdata[14:27] = vdata[12::-1, :]  # top/bottom flip
    aliased = np.where(vdata > 10.0)
    if alias:
        vdata[aliased] += -20.0
    fields["velocity"]["data"] = vdata

    radar.fields = fields
    return radar


def make_velocity_aliased_rhi_radar(alias=True):
    """
    Return a RHI radar with a target like reflectivity field.

    Set alias to False to return a de-aliased radar.
    """
    radar = make_empty_rhi_radar(50, 180, 1)
    radar.range["meters_between_gates"] = 1.0
    radar.range["meters_to_center_of_first_gate"] = 1.0
    radar.instrument_parameters = {"nyquist_velocity": {"data": np.array([10.0] * 360)}}

    fields = {
        "reflectivity": get_metadata("reflectivity"),
        "velocity": get_metadata("velocity"),
    }

    # fake reflectivity data, all zero reflectivity
    fdata = np.zeros((180, 50), dtype="float32")
    fields["reflectivity"]["data"] = fdata

    # fake velocity data, all zeros except a wind burst on at ~13 degrees.
    # burst is partially aliased.
    vdata = np.zeros((180 * 1, 50), dtype="float32")

    for i, idx in enumerate(range(13, -1, -1)):
        vdata[i, idx : idx + i + 1] = np.arange(0.5, 0.5 + i * 1.0 + 0.001)
    vdata[:14, 14:27] = vdata[:14, 12::-1]  # left/right flip
    vdata[14:27] = vdata[12::-1, :]  # top/bottom flip
    aliased = np.where(vdata > 10.0)
    if alias:
        vdata[aliased] += -20.0
    fields["velocity"]["data"] = vdata

    radar.fields = fields
    return radar


def make_single_ray_radar():
    """
    Return a PPI radar with a single ray taken from a ARM C-SAPR Radar.

    Radar object returned has 'reflectivity_horizontal',
    'norm_coherent_power', 'copol_coeff', 'dp_phase_shift', 'diff_phase', and
    'differential_reflectivity' fields with no metadata but a 'data' key.
    This radar is used for unit tests in correct modules.

    """
    radar = make_empty_ppi_radar(983, 1, 1)
    radar.range["data"] = 117.8784 + np.arange(983) * 119.91698
    f = np.load(_EXAMPLE_RAYS_FILE)
    for field_name in f:
        radar.fields[field_name] = {"data": f[field_name]}
    f.close()
    return radar


def make_empty_grid(grid_shape, grid_limits):
    """
    Make an empty grid object without any fields or metadata.

    Parameters
    ----------
    grid_shape : 3-tuple of floats
        Number of points in the grid (z, y, x).
    grid_limits : 3-tuple of 2-tuples
        Minimum and maximum grid location (inclusive) in meters for the
        z, y, x coordinates.

    Returns
    -------
    grid : Grid
        Empty Grid object, centered near the ARM SGP site (Oklahoma).

    """
    time = get_metadata("grid_time")
    time["data"] = np.array([0.0])
    time["units"] = "seconds since 2000-01-01T00:00:00Z"

    # grid coordinate dictionaries
    nz, ny, nx = grid_shape
    (z0, z1), (y0, y1), (x0, x1) = grid_limits

    x = get_metadata("x")
    x["data"] = np.linspace(x0, x1, nx)

    y = get_metadata("y")
    y["data"] = np.linspace(y0, y1, ny)

    z = get_metadata("z")
    z["data"] = np.linspace(z0, z1, nz)

    origin_altitude = get_metadata("origin_altitude")
    origin_altitude["data"] = np.array([300.0])

    origin_latitude = get_metadata("origin_latitude")
    origin_latitude["data"] = np.array([36.74])

    origin_longitude = get_metadata("origin_longitude")
    origin_longitude["data"] = np.array([-98.1])

    fields = {}
    metadata = {}

    radar_latitude = get_metadata("radar_latitude")
    radar_latitude["data"] = np.array([36.74])

    radar_longitude = get_metadata("radar_longitude")
    radar_longitude["data"] = np.array([-98.1])

    radar_altitude = get_metadata("radar_altitude")
    radar_altitude["data"] = np.array([300.0])

    radar_time = get_metadata("radar_time")
    radar_time["data"] = np.array([0.0])
    radar_time["units"] = "seconds since 2000-01-01T00:00:00Z"

    radar_name = get_metadata("radar_name")
    radar_name["data"] = np.array(["ExampleRadar"])

    return Grid(
        time,
        fields,
        metadata,
        origin_latitude,
        origin_longitude,
        origin_altitude,
        x,
        y,
        z,
        radar_latitude=radar_latitude,
        radar_longitude=radar_longitude,
        radar_altitude=radar_altitude,
        radar_time=radar_time,
        radar_name=radar_name,
    )


def make_target_grid():
    """
    Make a sample Grid with a rectangular target.
    """
    grid_shape = (2, 400, 320)
    grid_limits = ((0, 500), (-400000, 400000), (-300000, 300000))
    grid = make_empty_grid(grid_shape, grid_limits)

    fdata = np.zeros((2, 400, 320), dtype="float32")
    fdata[:, 50:-50, 40:-40] = 10.0
    fdata[:, 100:-100, 80:-80] = 20.0
    fdata[:, 150:-150, 120:-120] = 30.0
    fdata[1] += 5
    rdic = {"data": fdata, "long_name": "reflectivity", "units": "dBz"}
    grid.fields = {"reflectivity": rdic}
    return grid


def make_storm_grid():
    """
    Make a sample Grid with a rectangular storm target.
    """
    grid_shape = (2, 50, 40)
    grid_limits = ((0, 500), (-400000, 400000), (-300000, 300000))
    grid = make_empty_grid(grid_shape, grid_limits)

    fdata = np.ma.zeros((2, 50, 40), dtype="float32")
    fdata[:] = np.ma.masked
    fdata[:, 5:-5, 4:-4] = 20.0
    fdata[:, 10:-10, 8:-8] = 30.0
    fdata[:, 15:-15, 12:-12] = 60.0
    fdata[1] += 5
    rdic = {"data": fdata, "long_name": "reflectivity", "units": "dBz"}
    grid.fields = {"reflectivity": rdic}
    return grid


def make_normal_storm(sigma, mu):
    """
    Make a sample Grid with a gaussian storm target.
    """
    test_grid = make_empty_grid([1, 101, 101], [(1, 1), (-50, 50), (-50, 50)])
    x = test_grid.x["data"]
    y = test_grid.y["data"]
    z = test_grid.z["data"]
    zg, yg, xg = np.meshgrid(z, y, x, indexing="ij")
    r = np.sqrt((xg - mu[0]) ** 2 + (yg - mu[1]) ** 2)
    term1 = 1.0 / (sigma * np.sqrt(2.0 * np.pi))
    term2 = -1.0 * (r**2 / (2.0 * sigma**2))
    data = term1 * np.exp(term2)
    rdic = {"data": data, "long_name": "reflectivity", "units": "dBz"}
    test_grid.fields.update({"reflectivity": rdic})
    return test_grid


def make_empty_spectra_radar(nrays, ngates, npulses_max):
    """
    Return a Spectra Radar object.

    Parameters
    ----------
    nrays : int
        Number of rays in the object.
    ngates : int
        Number of gates per ray.
    npulses_max : int
        Number of pulses in each gate.

    Returns
    -------
    radar : RadarSpectra
        Radar spectra object with an empty spectra field, other parameters are
        set to default values.

    """
    if not _XARRAY_AVAILABLE:
        raise MissingOptionalDependency(
            "Xarray is required to use make_empty_spectra_radar "
            "but is not installed!"
        )

    time_dict = get_metadata("time")
    _range_dict = get_metadata("range")
    latitude_dict = get_metadata("latitude")
    longitude_dict = get_metadata("longitude")
    altitude_dict = get_metadata("altitude")
    sweep_number_dict = get_metadata("sweep_number")
    sweep_mode_dict = get_metadata("sweep_mode")
    fixed_angle_dict = get_metadata("fixed_angle")
    sweep_start_ray_index_dict = get_metadata("sweep_start_ray_index")
    sweep_end_ray_index_dict = get_metadata("sweep_end_ray_index")
    azimuth_dict = get_metadata("azimuth")
    elevation_dict = get_metadata("elevation")

    fields = {}
    scan_type = "vpt"
    c = 299792458
    wavelength = c / 34.830 * 1e-9
    metadata = {"instrument_name": "fake_spectra_radar", "wavelength": wavelength}
    time_dict["units"] = "seconds since 1989-01-01T00:00:01Z"
    time = xr.DataArray(np.arange(nrays, dtype="float32"), attrs=time_dict, dims="time")
    _range = xr.DataArray(
        np.linspace(0, 1000, ngates).astype("float32"), attrs=_range_dict, dims="range"
    )
    latitude = xr.DataArray(np.array([36.5], dtype="float64"), attrs=latitude_dict)
    longitude = xr.DataArray(np.array([-97.5], dtype="float64"), attrs=longitude_dict)
    altitude = xr.DataArray(np.array([200], dtype="float64"), attrs=altitude_dict)

    sweep_number = xr.DataArray(np.arange(1, dtype="int32"), attrs=sweep_number_dict)
    sweep_mode = xr.DataArray(np.array(["spectra"] * 1), attrs=sweep_mode_dict)
    fixed_angle = xr.DataArray(
        np.array([0.75] * 1, dtype="float32"), attrs=fixed_angle_dict
    )
    sweep_start_ray_index = xr.DataArray(
        np.array([0], dtype="int32"), attrs=sweep_start_ray_index_dict
    )
    sweep_end_ray_index = xr.DataArray(
        np.array([len(time.values) - 1], dtype="int32"), attrs=sweep_end_ray_index_dict
    )

    azimuth = xr.DataArray(
        np.arange(nrays, dtype="float32"), attrs=azimuth_dict, dims="time"
    )
    elevation = xr.DataArray(
        np.array([0.75] * nrays, dtype="float32"), attrs=elevation_dict, dims="time"
    )
    npulses_max = np.arange(npulses_max, dtype="float32")
    velocity_bins = xr.DataArray(
        np.linspace(-10.0, 10.0, len(npulses_max)), attrs={}, dims="npulses_max"
    )

    fields = np.zeros((len(time.values), len(_range.values), len(npulses_max)))

    return RadarSpectra(
        time,
        _range,
        fields,
        metadata,
        scan_type,
        latitude,
        longitude,
        altitude,
        sweep_number,
        sweep_mode,
        fixed_angle,
        sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth,
        elevation,
        npulses_max,
        velocity_bins,
        instrument_parameters=None,
    )


def make_target_spectra_radar():
    """
    Return a spectra radar with a target like spectra field.
    """
    if not _XARRAY_AVAILABLE:
        raise MissingOptionalDependency(
            "Xarray is required to use make_target_spectra_radar "
            "but is not installed!"
        )

    radar = make_empty_spectra_radar(10, 20, 50)
    fdata = np.zeros((10, 20, 50), dtype="float32")
    max_value = 10 ** (-10 / 10)
    fdata[:, :, :] = 10 * np.log10(scipy.signal.gaussian(50, std=7) * max_value)
    radar.ds["spectra"].values = fdata
    return radar
