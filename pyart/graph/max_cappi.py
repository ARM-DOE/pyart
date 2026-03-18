"""
Plot Max-CAPPI

This module provides a function to plot a Maximum Constant Altitude Plan Position Indicator (Max-CAPPI)
from radar data using an xarray dataset. The function includes options for adding map features, range rings,
color bars, and customized visual settings.

Author: Syed Hamid Ali (@syedhamidali)
"""

__all__ = ["plot_maxcappi"]

import os
import warnings

import cartopy.crs as ccrs
import cartopy.feature as feat
import matplotlib.pyplot as plt
import numpy as np
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from matplotlib.ticker import NullFormatter

warnings.filterwarnings("ignore")


def plot_maxcappi(
    grid,
    field,
    cmap=None,
    vmin=None,
    vmax=None,
    title=None,
    lat_lines=None,
    lon_lines=None,
    add_map=True,
    projection=None,
    colorbar=True,
    range_rings=False,
    dpi=100,
    savedir=None,
    show_figure=True,
    add_slogan=False,
    **kwargs,
):
    """
    Plot a Maximum Constant Altitude Plan Position Indicator (Max-CAPPI) using an xarray Dataset.

    Parameters
    ----------
    grid : pyart.core.Grid
        The grid object containing the radar data to be plotted.
    field : str
        The radar field to be plotted (e.g., "REF", "VEL", "WIDTH").
    cmap : str or matplotlib colormap, optional
        Colormap to use for the plot. Default is "HomeyerRainbow".
    vmin : float, optional
        Minimum value for the color scaling. Default is set to the minimum value of the data if not provided.
    vmax : float, optional
        Maximum value for the color scaling. Default is set to the maximum value of the data if not provided.
    title : str, optional
        Title of the plot. If None, the title is set to "Max-{field}".
    lat_lines : array-like, optional
        Latitude lines to be included in the plot. Default is calculated based on dataset coordinates.
    lon_lines : array-like, optional
        Longitude lines to be included in the plot. Default is calculated based on dataset coordinates.
    add_map : bool, optional
        Whether to include a map background in the plot. Default is True.
    projection : cartopy.crs.Projection, optional
        The map projection for the plot. Default is automatically determined based on dataset coordinates.
    colorbar : bool, optional
        Whether to include a colorbar in the plot. Default is True.
    range_rings : bool, optional
        Whether to include range rings at 50 km intervals. Default is False.
    dpi : int, optional
        DPI (dots per inch) for the plot. Default is 100.
    savedir : str, optional
        Directory where the plot will be saved. If None, the plot is not saved.
    show_figure : bool, optional
        Whether to display the plot. Default is True.
    add_slogan : bool, optional
        Whether to add a slogan like "Powered by Py-ART" to the plot. Default is False.
    **kwargs : dict, optional
        Additional keyword arguments to pass to matplotlib's `pcolormesh` function.

    Returns
    -------
    None
        This function does not return any value. It generates and optionally displays or saves a plot.

    Notes
    -----
    - The function extracts the maximum value across the altitude (z) dimension to create the Max-CAPPI.
    - It supports customizations such as map projections, color scales, and range rings.
    - If the radar_name attribute in the dataset is a byte string, it will be decoded and limited to 4 characters.
    - If add_map is True, map features and latitude/longitude lines are included.
    - The plot can be saved to a specified directory in PNG format.

    Author: Syed Hamid Ali (@syedhamidali)
    """

    ds = grid.to_xarray().squeeze()

    if lon_lines is None:
        lon_lines = np.arange(int(ds.lon.min().values), int(ds.lon.max().values) + 1)
    if lat_lines is None:
        lat_lines = np.arange(int(ds.lat.min().values), int(ds.lat.max().values) + 1)

    plt.rcParams.copy()
    plt.rcParams.update(
        {
            "font.weight": "bold",
            "axes.labelweight": "bold",
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 10,
            "ytick.major.size": 10,
            "xtick.minor.size": 7,
            "ytick.minor.size": 7,
            "font.size": 14,
            "axes.linewidth": 2,
            "ytick.labelsize": 12,
            "xtick.labelsize": 12,
        }
    )

    max_c = ds[field].max(dim="z")
    max_x = ds[field].max(dim="y")
    max_y = ds[field].max(dim="x").T

    trgx = ds["x"].values
    trgy = ds["y"].values
    trgz = ds["z"].values

    max_height = int(np.floor(trgz.max()) / 1e3)
    sideticks = np.arange(max_height / 4, max_height + 1, max_height / 4).astype(int)

    if cmap is None:
        cmap = "HomeyerRainbow"
    if vmin is None:
        vmin = grid.fields[field]["data"].min()
    if vmax is None:
        vmax = grid.fields[field]["data"].max()
    if title is None:
        title = f"Max-{field.upper()[:3]}"

    def plot_range_rings(ax_xy, max_range):
        """
        Plots range rings at 50 km intervals.

        Parameters
        ----------
        ax_xy : matplotlib.axes.Axes
            The axis on which to plot the range rings.
        max_range : float
            The maximum range for the range rings.

        Returns
        -------
        None
        """
        background_color = ax_xy.get_facecolor()
        color = "k" if sum(background_color[:3]) / 3 > 0.5 else "w"

        for i, r in enumerate(np.arange(5e4, np.floor(max_range) + 1, 5e4)):
            label = f"Ring Dist. {int(r/1e3)} km" if i == 0 else None
            ax_xy.plot(
                r * np.cos(np.arange(0, 360) * np.pi / 180),
                r * np.sin(np.arange(0, 360) * np.pi / 180),
                color=color,
                ls="--",
                linewidth=0.4,
                alpha=0.3,
                label=label,
            )

        ax_xy.legend(loc="upper right", prop={"weight": "normal", "size": 8})

    def _get_projection(ds):
        """
        Determine the central latitude and longitude from a dataset
        and return the corresponding projection.

        Parameters
        ----------
        ds : xarray.Dataset
            The dataset from which to extract latitude and longitude
            information.

        Returns
        -------
        projection : cartopy.crs.Projection
            A Cartopy projection object centered on the extracted or
            calculated latitude and longitude.
        """

        def get_coord_or_attr(ds, coord_name, attr_name):
            """Helper function to get a coordinate or attribute, or
            calculate median if available.
            """
            if coord_name in ds:
                return (
                    ds[coord_name].values.item()
                    if ds[coord_name].values.ndim == 0
                    else ds[coord_name].values[0]
                )
            if f"origin_{coord_name}" in ds.coords:
                return ds.coords[f"origin_{coord_name}"].median().item()
            if f"radar_{coord_name}" in ds.coords:
                return ds.coords[f"radar_{coord_name}"].median().item()
            return ds.attrs.get(attr_name, None)

        lat_0 = get_coord_or_attr(
            ds, "latitude", "origin_latitude"
        ) or get_coord_or_attr(ds, "radar_latitude", "origin_latitude")
        lon_0 = get_coord_or_attr(
            ds, "longitude", "origin_longitude"
        ) or get_coord_or_attr(ds, "radar_longitude", "origin_longitude")

        if lat_0 is None or lon_0 is None:
            lat_0 = ds.lat.mean().item()
            lon_0 = ds.lon.mean().item()

        projection = ccrs.LambertAzimuthalEqualArea(lon_0, lat_0)
        return projection

    projection = _get_projection(ds)

    # FIG
    fig = plt.figure(figsize=[10.3, 10])
    left, bottom, width, height = 0.1, 0.1, 0.6, 0.2
    ax_xy = plt.axes((left, bottom, width, width), projection=projection)
    ax_x = plt.axes((left, bottom + width, width, height))
    ax_y = plt.axes((left + width, bottom, height, width))
    ax_cnr = plt.axes((left + width, bottom + width, left + left, height))
    if colorbar:
        ax_cb = plt.axes((left - 0.015 + width + height + 0.02, bottom, 0.02, width))

    # Set axis label formatters
    ax_x.xaxis.set_major_formatter(NullFormatter())
    ax_y.yaxis.set_major_formatter(NullFormatter())
    ax_cnr.yaxis.set_major_formatter(NullFormatter())
    ax_cnr.xaxis.set_major_formatter(NullFormatter())
    ax_x.set_ylabel("Height (km)", size=13)
    ax_y.set_xlabel("Height (km)", size=13)

    # Draw CAPPI
    plt.sca(ax_xy)
    xy = ax_xy.pcolormesh(trgx, trgy, max_c, cmap=cmap, vmin=vmin, vmax=vmax, **kwargs)

    # Add map features
    if add_map:
        map_features(ax_xy, lat_lines, lon_lines)

    ax_xy.minorticks_on()

    if range_rings:
        plot_range_rings(ax_xy, trgx.max())

    ax_xy.set_xlim(trgx.min(), trgx.max())
    ax_xy.set_ylim(trgx.min(), trgx.max())

    # Draw colorbar
    if colorbar:
        cb = plt.colorbar(xy, cax=ax_cb)
        cb.set_label(ds[field].attrs["units"], size=15)

    background_color = ax_xy.get_facecolor()
    color = "k" if sum(background_color[:3]) / 3 > 0.5 else "w"

    plt.sca(ax_x)
    plt.pcolormesh(trgx / 1e3, trgz / 1e3, max_x, cmap=cmap, vmin=vmin, vmax=vmax)
    plt.yticks(sideticks)
    ax_x.set_xlim(trgx.min() / 1e3, trgx.max() / 1e3)
    ax_x.grid(axis="y", lw=0.5, color=color, alpha=0.5, ls=":")
    ax_x.minorticks_on()

    plt.sca(ax_y)
    plt.pcolormesh(trgz / 1e3, trgy / 1e3, max_y, cmap=cmap, vmin=vmin, vmax=vmax)
    ax_y.set_xticks(sideticks)
    ax_y.set_ylim(trgx.min() / 1e3, trgx.max() / 1e3)
    ax_y.grid(axis="x", lw=0.5, color=color, alpha=0.5, ls=":")
    ax_y.minorticks_on()

    plt.sca(ax_cnr)
    plt.tick_params(
        axis="both",  # changes apply to both axes
        which="both",  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        left=False,
        right=False,
        labelbottom=False,
    )

    # Initialize an empty list to store the processed radar names
    full_title = []

    # Check if radar_name is a list (or list-like) or a simple string
    if isinstance(ds.attrs["radar_name"], list):
        # Iterate over each radar name in the list
        for name in ds.attrs["radar_name"]:
            # Decode if it's a byte string and take the first 4 characters
            if isinstance(name, bytes):
                site_title = name.decode("utf-8")[:4]
            else:
                site_title = name[:4]
            full_title.append(site_title)
    else:
        # Handle the case where radar_name is a single string
        site_title = ds.attrs["radar_name"][:4]
        full_title.append(site_title)

    # Join the processed radar names into a single string with commas separating them
    formatted_title = ", ".join(full_title)

    # Center-align text in the corner box
    plt.text(
        0.5,
        0.90,
        f"{formatted_title}",
        size=13,
        weight="bold",
        ha="center",
        va="center",
    )
    plt.text(0.5, 0.76, title, size=13, weight="bold", ha="center", va="center")
    plt.text(
        0.5,
        0.63,
        f"Max Range: {np.floor(trgx.max() / 1e3)} km",
        size=11,
        ha="center",
        va="center",
    )
    plt.text(
        0.5,
        0.47,
        f"Max Height: {np.floor(trgz.max() / 1e3)} km",
        size=11,
        ha="center",
        va="center",
    )
    plt.text(
        0.5,
        0.28,
        ds["time"].dt.strftime("%H:%M:%S Z").values.item(),
        weight="bold",
        size=16,
        ha="center",
        va="center",
    )
    plt.text(
        0.5,
        0.13,
        ds["time"].dt.strftime("%d %b, %Y UTC").values.item(),
        size=13.5,
        ha="center",
        va="center",
    )
    ax_xy.set_aspect("auto")

    if add_slogan:
        fig.text(
            0.1,
            0.06,
            "Powered by Py-ART",  # Coordinates close to (0, 0) for lower-left corner
            fontsize=9,
            fontname="Courier New",
            # bbox=dict(facecolor='none', boxstyle='round,pad=0.5')
        )

    if savedir is not None:
        radar_name = ds.attrs.get("instrument_name", "Radar")
        time_str = ds["time"].dt.strftime("%Y%m%d%H%M%S").values.item()
        figname = f"{savedir}{os.sep}{title}_{radar_name}_{time_str}.png"
        plt.savefig(fname=figname, dpi=dpi, bbox_inches="tight")
        print(f"Figure(s) saved as {figname}")

    # plt.rcParams.update(original_rc_params)
    plt.rcdefaults()

    if show_figure:
        plt.show()
    else:
        plt.close()


def map_features(ax, lat_lines, lon_lines):
    """
    Adds map features and gridlines to the plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis on which to add map features and gridlines.
    lat_lines : array-like
        Latitude lines for the gridlines.
    lon_lines : array-like
        Longitude lines for the gridlines.

    Returns
    -------
    None
    """
    background_color = ax.get_facecolor()
    color = "k" if sum(background_color[:3]) / 3 > 0.5 else "w"

    # Labeling gridlines depending on the projection
    if isinstance(ax.projection, (ccrs.PlateCarree, ccrs.Mercator)):
        gl = ax.gridlines(
            xlocs=lon_lines,
            ylocs=lat_lines,
            linewidth=1,
            alpha=0.5,
            linestyle="--",
            draw_labels=True,
        )
        gl.top_labels = False
        gl.right_labels = False
        gl.xlines = False
        gl.ylines = False
        gl.xlabel_style = {"color": color}
        gl.ylabel_style = {"color": color}
        ax.add_feature(feat.COASTLINE, alpha=0.8, lw=1, ec=color)
        ax.add_feature(feat.BORDERS, alpha=0.7, lw=0.7, ls="--", ec=color)
        ax.add_feature(
            feat.STATES.with_scale("10m"), alpha=0.6, lw=0.5, ls=":", ec=color
        )

    elif isinstance(
        ax.projection, (ccrs.LambertConformal, ccrs.LambertAzimuthalEqualArea)
    ):
        ax.figure.canvas.draw()
        gl = ax.gridlines(
            crs=ccrs.PlateCarree(),
            xlocs=lon_lines,
            ylocs=lat_lines,
            linewidth=1,
            alpha=0.5,
            linestyle="--",
            draw_labels=False,
        )
        gl.xlines = False
        gl.ylines = False
        gl.xlabel_style = {"color": color}
        gl.ylabel_style = {"color": color}
        ax.add_feature(feat.COASTLINE, alpha=0.8, lw=1, ec=color)
        ax.add_feature(feat.BORDERS, alpha=0.7, lw=0.7, ls="--", ec=color)
        ax.add_feature(
            feat.STATES.with_scale("10m"), alpha=0.6, lw=0.5, ls=":", ec=color
        )
        # Label the end-points of the gridlines using custom tick makers
        ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
        ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
        from pyart.graph.gridmapdisplay import lambert_xticks, lambert_yticks

        lambert_xticks(ax, lon_lines)
        lambert_yticks(ax, lat_lines)
    else:
        ax.gridlines(xlocs=lon_lines, ylocs=lat_lines)
