"""
Functions for gridding and processing radar sweeps in a Cartesian grid

"""

import numpy as np
import xarray as xr
from ..map import grid_from_radars


def grid_ppi_sweeps(radar,
                    target_sweeps=None,
                    grid_size=801,
                    grid_limits='auto',
                    max_z=12000.,
                    el_rounding_frac=0.5,
                    add_grid_altitude=True,
                    **kwargs,
                   ):
    """
    Separately grid PPI sweeps to an X-Y plane considering only horizontal distances
    in grid RoI and weighting function.
    Gridding is performed using the `grid_from_radars` method, which receives any
    additional input parameters.
    Note that `h_factor` and `dist_factor` should not be included in kwargs (required
    for valid gridding results)

    Parameters
    ----------
    radar : Radar
        Radar volume containing PPI sweeps.
    target_sweeps : int or list
        sweeps to grid. Using all sweeps in `radar` if None.
    grid_size: int or 2-tuple
        grid dimension size. Using sizes for the X-Y plane if tuple.
        This input parameter is ignored if `grid_shape` is given
        explicitly via kwargs.
    grid_limits: 3-tuple with 2-tuple elements or 'auto'
        if 'auto' using the maximum horizontal range rounded up to the nearest kilometer
        and limiting vertically up to `max_z`.
    max_z: float
        maximum height to consider in gridding (only used if `grid_size` is 'auto')
    el_rounding_frac: float
        A fraction for rounding the elevation elements. This variables is also used to
        represent the sweep for altitude estimation.
    add_grid_altitude: bool
        adding a sweep-dependent altitude estimate corresponding to the X-Y plane if True.
        This output field is useful considering the slanted PPI scans.

    Returns
    -------
    radar_ds : xarray.Dataset
        Radar data gridded to the X-Y plane with a third dimension
        representing the different sweep elevations.

    """
    if target_sweeps is None:
        target_sweeps = radar.sweep_number['data']
    elif isinstance(target_sweeps, int):
        target_sweeps = [target_sweeps]

    # Set grid shape
    if not 'grid_shape' in kwargs.keys():
        if isinstance(grid_size, int):
            grid_shape = (1, grid_size, grid_size)
        elif isinstance(grid_size, tuple):
            grid_shape = (1, *grid_size)
        else:
            raise TypeError("'grid_shape' must be of type int or tuple")

    # Set grid limits in 'auto' option
    if isinstance (grid_limits, str):
        if grid_limits == 'auto':
            max_xy = np.max([np.max(radar.get_gate_x_y_z(sweep=sw)[0]) for sw in target_sweeps])
            max_xy = np.ceil(max_xy / 1e3) * 1e3
            grid_limits = ((0., max_z), (-max_xy, max_xy), (-max_xy, max_xy))
        else:
            raise ValueError(f"Unknown 'grid_limits' processing string {grid_limits}")

    # Calling the gridding method
    radar_ds = None
    for sweep in target_sweeps:
        radar_sw = radar.extract_sweeps([sweep])
        sweep_grid = grid_from_radars(
            (radar_sw,),
            grid_shape=grid_shape,
            grid_limits=grid_limits,
            h_factor=(0.0, 1.0, 1.0),
            dist_factor=(0.0, 1.0, 1.0),
            **kwargs,
        )

        # Convert to xarray.Dataset and finalize
        el_round = (np.mean(radar_sw.elevation['data']) /
                    el_rounding_frac).round() * el_rounding_frac
        radar_ds_tmp = sweep_grid.to_xarray().squeeze()
        if add_grid_altitude:
            alt_est = ((radar_ds_tmp["x"] ** 2 + radar_ds_tmp["y"] ** 2) ** 0.5 *
                       np.tan(el_round * np.pi / 180.))
        radar_ds_tmp["altitude_est"] = xr.DataArray(
            alt_est,
            coords=radar_ds_tmp.coords,
            dims=radar_ds_tmp.dims,
            attrs={"long_name": "Estimated altitude in PPI scan",
                   "units": "m"}
        )
        radar_ds_tmp = radar_ds_tmp.expand_dims(elevation=[el_round])
        radar_ds["elevation"].attrs = {"long_name": "Elevation angle",
                                       "units": "deg"}
        if radar_ds is None:
            radar_ds = radar_ds_tmp
        else:
            radar_ds = xr.concat((radar_ds, radar_ds_tmp), dim="elevation")

    return radar_ds


def grid_rhi_sweeps(radar,
                    target_sweeps=None,
                    grid_size=801,
                    grid_limits='auto',
                    max_z=12000.,
                    az_rounding_frac=0.5,
                    **kwargs,
                   ):
    """
    Separately grid RHI sweeps to a Y-Z plane considering only cross-sectional distances
    in grid RoI and weighting function.
    Gridding is performed using the `grid_from_radars` method, which receives any
    additional input parameters.
    Note that `h_factor` and `dist_factor` should not be included in kwargs (required
    for valid gridding results)

    Parameters
    ----------
    radar : Radar
        Radar volume containing PPI sweeps.
    target_sweeps : int or list
        sweeps to grid. Using all sweeps in `radar` if None.
    grid_size: int or 2-tuple
        grid dimension size. Using sizes for the Y-Z plane if tuple.
        This input parameter is ignored if `grid_shape` is given
        explicitly via kwargs.
    max_z: float
        maximum height in grid (only used if `grid_size` is 'auto').
    grid_limits: 3-tuple with 2-tuple elements or 'auto'
        if 'auto' using the maximum horizontal range and limiting vertically up to 12 km.
    az_rounding_frac: float
        A fraction for rounding the azimuth elements.

    Returns
    -------
    radar_ds : xarray.Dataset
        Radar data gridded to the Y-Z plane with a third dimension
        representing the different sweep azimuths.

    """
    if target_sweeps is None:
        target_sweeps = radar.sweep_number['data']
    elif isinstance(target_sweeps, int):
        target_sweeps = [target_sweeps]

    # Set grid shape
    if not 'grid_shape' in kwargs.keys():
        if isinstance(grid_size, int):
            grid_shape = (grid_shape, grid_size, 1)
        elif isinstance(grid_size, tuple):
            grid_shape = (*grid_size, 1)
        else:
            raise TypeError("'grid_shape' must be of type int or tuple")

    # Set grid limits in 'auto' option
    if isinstance (grid_limits, str):
        if grid_limits == 'auto':
            max_xy = np.max([np.max(radar.get_gate_x_y_z(sweep=sw)[0]) for sw in target_sweeps])
            max_xy = np.ceil(max_xy / 1e3) * 1e3
            grid_limits = ((0., max_z), (-max_xy, max_xy), (-max_xy, max_xy))
        else:
            raise ValueError(f"Unknown 'grid_limits' processing string {grid_limits}")

    # Calling the gridding method
    radar_ds = None
    for sweep in target_sweeps:
        az_round = (np.mean(radar_sw.azimuth['data']) /
                    az_rounding_frac).round() * az_rounding_frac
        radar_sw = radar.extract_sweeps([sweep])
        radar_sw.azimuth['data'] -= az_round  # centering azimuth values to 0 to maximize grid range
        sweep_grid = grid_from_radars(
            (radar_sw,),
            grid_shape=grid_shape,
            grid_limits=grid_limits,
            h_factor=(1.0, 1.0, 0.0),
            dist_factor=(1.0, 1.0, 0.0),
            **kwargs,
        )

        # Convert to xarray.Dataset and finalize
        radar_ds_tmp = sweep_grid.to_xarray().squeeze()
        radar_ds_tmp = radar_ds_tmp.expand_dims(azimuth=[az_round])
        radar_ds["azimuth"].attrs = {"long_name": "Azimuth angle",
                                     "units": "deg"}
        if radar_ds is None:
            radar_ds = radar_ds_tmp
        else:
            radar_ds = xr.concat((radar_ds, radar_ds_tmp), dim="azimuth")

    return radar_ds
