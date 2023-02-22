"""
A class for plotting grid objects using xarray plotting
and cartopy.

"""

import warnings

import matplotlib.pyplot as plt
import numpy as np

try:
    import cartopy  # noqa

    _CARTOPY_AVAILABLE = True
except ImportError:
    _CARTOPY_AVAILABLE = False

try:
    import metpy  # noqa

    _METPY_AVAILABLE = True
except ImportError:
    _METPY_AVAILABLE = False

from pyart.core.transforms import _interpolate_axes_edges
from pyart.exceptions import MissingOptionalDependency
from pyart.graph import common

try:
    import xarray  # noqa

    _XARRAY_AVAILABLE = True
except ImportError:
    _XARRAY_AVAILABLE = False

try:
    import netCDF4  # noqa

    _NETCDF4_AVAILABLE = True
except ImportError:
    _NETCDF4_AVAILABLE = False

try:
    from copy import copy

    import shapely.geometry as sgeom

    _LAMBERT_GRIDLINES = True
except ImportError:
    _LAMBERT_GRIDLINES = False


class GridMapDisplay:
    """
    A class for creating plots from a grid object using xarray
    with a cartopy projection.

    Parameters
    ----------
    grid : Grid
        Grid with data which will be used to create plots.
    debug : bool
        True to print debugging messages, False to supress them.

    Attributes
    ----------
    grid : Grid
        Grid object.
    debug : bool
        True to print debugging messages, False to supress them.

    """

    def __init__(self, grid, debug=False):
        """initalize the object."""
        # check that cartopy and xarray are available
        if not _CARTOPY_AVAILABLE:
            raise MissingOptionalDependency(
                "Cartopy is required to use GridMapDisplay but is not installed!"
            )
        if not _XARRAY_AVAILABLE:
            raise MissingOptionalDependency(
                "Xarray is required to use GridMapDisplay but is not installed!"
            )
        if not _NETCDF4_AVAILABLE:
            raise MissingOptionalDependency(
                "netCDF4 is required to use GridMapDisplay but is not installed!"
            )

        # set attributes
        self.grid = grid
        self.debug = debug
        self.mappables = []
        self.fields = []
        self.origin = "origin"

    def plot_grid(
        self,
        field,
        level=0,
        vmin=None,
        vmax=None,
        norm=None,
        cmap=None,
        mask_outside=False,
        title=None,
        title_flag=True,
        axislabels=(None, None),
        axislabels_flag=False,
        colorbar_flag=True,
        colorbar_label=None,
        colorbar_orient="vertical",
        ax=None,
        fig=None,
        lat_lines=None,
        lon_lines=None,
        projection=None,
        embellish=True,
        add_grid_lines=True,
        ticks=None,
        ticklabs=None,
        **kwargs
    ):
        """
        Plot the grid using xarray and cartopy.

        Additional arguments are passed to Xarray's pcolormesh function.

        Parameters
        ----------
        field : str
            Field to be plotted.
        level : int
            Index corresponding to the height level to be plotted.

        Other Parameters
        ----------------
        vmin, vmax : float
            Lower and upper range for the colormesh. If either parameter is
            None, a value will be determined from the field attributes (if
            available) or the default values of -8, 64 will be used.
            Parameters are used for luminance scaling.
        norm : Normalize or None, optional
            matplotlib Normalize instance used to scale luminance data. If not
            None the vmax and vmin parameters are ignored. If None, vmin and
            vmax are used for luminance scaling.
        cmap : str or None
            Matplotlib colormap name. None will use default colormap for
            the field being plotted as specified by the Py-ART configuration.
        mask_outside : bool
            True to mask data outside of vmin, vmax. False performs no
            masking.
        title : str
            Title to label plot with, None will use the default generated from
            the field and level parameters. Parameter is ignored if the
            title_flag is False.
        title_flag : bool
            True to add title to plot, False does not add a title.
        axislabels : (str, str)
            2-tuple of x-axis, y-axis labels. None for either label will use
            the default axis label. Parameter is ignored if axislabels_flag is
            False.
        axislabels_flag : bool
            True to add label the axes, False does not label the axes.
        colorbar_flag : bool
            True to add a colorbar with label to the axis. False leaves off
            the colorbar.
        colorbar_label : str
            Colorbar label, None will use a default label generated from the
            field information.
        colorbar_orient : 'vertical' or 'horizontal'
            Colorbar orientation.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.
        lat_lines, lon_lines : array or None
            Location at which to draw latitude and longitude lines.
            None will use default values which are reasonable for maps of
            North America.
        projection : cartopy.crs class
            Map projection supported by cartopy. Used for all subsequent calls
            to the GeoAxes object generated. Defaults to PlateCarree.
        embellish : bool
            True by default. Set to False to supress drawing of coastlines
            etc... Use for speedup when specifying shapefiles.
        add_grid_lines : bool
            True by default. Set to False to supress drawing of lat/lon lines
            Note that lat lon labels only work with certain projections.
        ticks : array
            Colorbar custom tick label locations.
        ticklabs : array
            Colorbar custom tick labels.

        """
        ds = self.grid.to_xarray()

        # parse parameters
        vmin, vmax = common.parse_vmin_vmax(self.grid, field, vmin, vmax)
        cmap = common.parse_cmap(cmap, field)

        # mask the data where outside the limits
        if mask_outside:
            data = ds[field].data
            masked_data = np.ma.masked_invalid(data)
            masked_data = np.ma.masked_outside(masked_data, vmin, vmax)
            ds[field].data = masked_data

        # Define a figure if None is provided.
        if fig is None:
            fig = plt.gcf()

        # initialize instance of GeoAxes if not provided
        if ax is not None:
            if hasattr(ax, "projection"):
                projection = ax.projection
            else:
                if projection is None:
                    # set map projection to Mercator if none is
                    # specified.
                    projection = cartopy.crs.Mercator()
                    warnings.warn(
                        "No projection was defined for the axes."
                        + " Overridding defined axes and using default "
                        + "axes with projection Mercator.",
                        UserWarning,
                    )
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore")
                    ax = plt.axes(projection=projection)

        # Define GeoAxes if None is provided.
        else:
            if projection is None:
                # set map projection to LambertConformal if none is
                # specified.
                projection = cartopy.crs.Mercator()
                warnings.warn(
                    "No projection was defined for the axes."
                    + " Overridding defined axes and using default "
                    + "axes with projection Mercator.",
                    UserWarning,
                )
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                ax = plt.axes(projection=projection)

        # plot the grid using xarray
        if norm is not None:  # if norm is set do not override with vmin/vmax
            vmin = vmax = None

        pm = ds[field][0, level].plot.pcolormesh(
            x="lon",
            y="lat",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            add_colorbar=False,
            **kwargs
        )

        self.mappables.append(pm)
        self.fields.append(field)

        if embellish:
            # Create a feature for States/Admin 1 regions at 1:50m
            # from Natural Earth
            states = self.cartopy_states()
            coastlines = self.cartopy_coastlines()
            ax.add_feature(states, linestyle="-", edgecolor="k", linewidth=2)
            ax.add_feature(coastlines, linestyle="-", edgecolor="k", linewidth=2)

        if add_grid_lines:
            if lon_lines is None:
                lon_lines = np.linspace(
                    np.around(ds.lon.min() - 0.1, decimals=2),
                    np.around(ds.lon.max() + 0.1, decimals=2),
                    5,
                )
            if lat_lines is None:
                lat_lines = np.linspace(
                    np.around(ds.lat.min() - 0.1, decimals=2),
                    np.around(ds.lat.max() + 0.1, decimals=2),
                    5,
                )

            # labeling gridlines poses some difficulties depending on the
            # projection, so we need some projection-specific methods
            if ax.projection in [cartopy.crs.PlateCarree(), cartopy.crs.Mercator()]:
                ax.gridlines(
                    draw_labels=False,
                    linewidth=2,
                    color="gray",
                    alpha=0.5,
                    linestyle="--",
                    xlocs=lon_lines,
                    ylocs=lat_lines,
                )
                ax.set_extent(
                    [
                        lon_lines.min(),
                        lon_lines.max(),
                        lat_lines.min(),
                        lat_lines.max(),
                    ],
                    crs=projection,
                )
                ax.set_xticks(lon_lines, crs=projection)
                ax.set_yticks(lat_lines, crs=projection)

            elif isinstance(ax.projection, cartopy.crs.LambertConformal):
                ax.figure.canvas.draw()
                ax.gridlines(xlocs=lon_lines, ylocs=lat_lines)

                # Label the end-points of the gridlines using the custom
                # tick makers:
                ax.xaxis.set_major_formatter(cartopy.mpl.gridliner.LONGITUDE_FORMATTER)
                ax.yaxis.set_major_formatter(cartopy.mpl.gridliner.LATITUDE_FORMATTER)
                if _LAMBERT_GRIDLINES:
                    lambert_xticks(ax, lon_lines)
                    lambert_yticks(ax, lat_lines)
            else:
                ax.gridlines(xlocs=lon_lines, ylocs=lat_lines)

        if title_flag:
            if title is None:
                ax.set_title(self.generate_grid_title(field, level))
            else:
                ax.set_title(title)

        if axislabels_flag:
            self._label_axes_grid(axislabels, ax)

        if colorbar_flag:
            self.plot_colorbar(
                mappable=pm,
                label=colorbar_label,
                orientation=colorbar_orient,
                field=field,
                ax=ax,
                fig=fig,
                ticks=ticks,
                ticklabs=ticklabs,
            )

    def plot_crosshairs(
        self, lon=None, lat=None, linestyle="--", color="r", linewidth=2, ax=None
    ):
        """
        Plot crosshairs at a given longitude and latitude.

        Parameters
        ----------
        lon, lat : float
            Longitude and latitude (in degrees) where the crosshairs should
            be placed. If None the center of the grid is used.
        linestyle : str
            Matplotlib string describing the line style.
        color : str
            Matplotlib string for color of the line.
        linewidth : float
            Width of markers in points.
        ax : axes or None
            Axis to add the crosshairs to, if None the current axis is used.

        """
        # parse the parameters
        ax = common.parse_ax(ax)
        lon, lat = common.parse_lon_lat(self.grid, lon, lat)

        # add crosshairs
        ax.axhline(lat, color=color, linestyle=linestyle, linewidth=linewidth)
        ax.axvline(lon, color=color, linestyle=linestyle, linewidth=linewidth)

    def plot_latitude_slice(self, field, lon=None, lat=None, **kwargs):
        """
        Plot a slice along a given latitude.

        For documentation of additional arguments see
        :py:func:`plot_latitudinal_level`.

        Parameters
        ----------
        field : str
            Field to be plotted.
        lon, lat : float
            Longitude and latitude (in degrees) specifying the slice. If
            None the center of the grid is used.

        """
        # parse parameters
        _, y_index = self._find_nearest_grid_indices(lon, lat)
        self.plot_latitudinal_level(field=field, y_index=y_index, **kwargs)

    def plot_latitudinal_level(
        self,
        field,
        y_index,
        vmin=None,
        vmax=None,
        norm=None,
        cmap=None,
        mask_outside=False,
        title=None,
        title_flag=True,
        axislabels=(None, None),
        axislabels_flag=True,
        colorbar_flag=True,
        colorbar_label=None,
        colorbar_orient="vertical",
        edges=True,
        ax=None,
        fig=None,
        ticks=None,
        ticklabs=None,
        **kwargs
    ):
        """
        Plot a slice along a given latitude.

        Additional arguments are passed to Basemaps's pcolormesh function.

        Parameters
        ----------
        field : str
            Field to be plotted.
        y_index : float
            Index of the latitudinal level to plot.
        vmin, vmax : float
            Lower and upper range for the colormesh. If either parameter is
            None, a value will be determined from the field attributes (if
            available) or the default values of -8, 64 will be used.
            Parameters are ignored is norm is not None.
        norm : Normalize or None, optional
            matplotlib Normalize instance used to scale luminance data. If not
            None the vmax and vmin parameters are ignored. If None, vmin and
            vmax are used for luminance scaling.
        cmap : str or None
            Matplotlib colormap name. None will use the default colormap for
            the field being plotted as specified by the Py-ART configuration.
        mask_outside : bool
            True to mask data outside of vmin, vmax. False performs no
            masking.
        title : str
            Title to label plot with, None to use default title generated from
            the field and lat,lon parameters. Parameter is ignored if
            title_flag is False.
        title_flag : bool
            True to add a title to the plot, False does not add a title.
        axislabels : (str, str)
            2-tuple of x-axis, y-axis labels. None for either label will use
            the default axis label. Parameter is ignored if axislabels_flag is
            False.
        axislabels_flag : bool
            True to add label the axes, False does not label the axes.
        colorbar_flag : bool
            True to add a colorbar with label to the axis.  False leaves off
            the colorbar.
        colorbar_label : str
            Colorbar label, None will use a default label generated from the
            field information.
        colorbar_orient : 'vertical' or 'horizontal'
            Colorbar orientation.
        edges : bool
            True will interpolate and extrapolate the gate edges from the
            range, azimuth and elevations in the radar, treating these
            as specifying the center of each gate. False treats these
            coordinates themselved as the gate edges, resulting in a plot
            in which the last gate in each ray and the entire last ray are not
            not plotted.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.
        ticks : array
            Colorbar custom tick label locations.
        ticklabs : array
            Colorbar custom tick labels.

        """
        # parse parameters
        ax, fig = common.parse_ax_fig(ax, fig)
        vmin, vmax = common.parse_vmin_vmax(self.grid, field, vmin, vmax)
        cmap = common.parse_cmap(cmap, field)

        data = self.grid.fields[field]["data"][:, y_index, :]

        # mask the data where outside the limits
        if mask_outside:
            data = np.ma.masked_invalid(data)
            data = np.ma.masked_outside(data, vmin, vmax)

        # plot the grid
        x_1d = self.grid.x["data"] / 1000
        z_1d = self.grid.z["data"] / 1000

        if edges:
            if len(x_1d) > 1:
                x_1d = _interpolate_axes_edges(x_1d)
            if len(z_1d) > 1:
                z_1d = _interpolate_axes_edges(z_1d)
        xd, yd = np.meshgrid(x_1d, z_1d)
        if norm is not None:  # if norm is set do not override with vmin, vmax
            vmin = vmax = None

        pm = ax.pcolormesh(
            xd, yd, data, vmin=vmin, vmax=vmax, norm=norm, cmap=cmap, **kwargs
        )

        self.mappables.append(pm)
        self.fields.append(field)

        if title_flag:
            if title is None:
                ax.set_title(
                    common.generate_latitudinal_level_title(self.grid, field, y_index)
                )
            else:
                ax.set_title(title)

        if axislabels_flag:
            self._label_axes_latitude(axislabels, ax)

        if colorbar_flag:
            self.plot_colorbar(
                mappable=pm,
                label=colorbar_label,
                orientation=colorbar_orient,
                field=field,
                ax=ax,
                fig=fig,
                ticks=ticks,
                ticklabs=ticklabs,
            )

    def plot_longitude_slice(self, field, lon=None, lat=None, **kwargs):
        """
        Plot a slice along a given longitude.

        For documentation of additional arguments see
        :py:func:`plot_longitudinal_level`.

        Parameters
        ----------
        field : str
            Field to be plotted.
        lon, lat : float
            Longitude and latitude (in degrees) specifying the slice.  If
            None the center of the grid is used.

        """
        # parse parameters
        x_index, _ = self._find_nearest_grid_indices(lon, lat)
        self.plot_longitudinal_level(field=field, x_index=x_index, **kwargs)

    def plot_longitudinal_level(
        self,
        field,
        x_index,
        vmin=None,
        vmax=None,
        norm=None,
        cmap=None,
        mask_outside=False,
        title=None,
        title_flag=True,
        axislabels=(None, None),
        axislabels_flag=True,
        colorbar_flag=True,
        colorbar_label=None,
        colorbar_orient="vertical",
        edges=True,
        ax=None,
        fig=None,
        ticks=None,
        ticklabs=None,
        **kwargs
    ):
        """
        Plot a slice along a given longitude.

        Additional arguments are passed to Basemaps's pcolormesh function.

        Parameters
        ----------
        field : str
            Field to be plotted.
        x_index : float
            Index of the longitudinal level to plot.
        vmin, vmax : float
            Lower and upper range for the colormesh. If either parameter is
            None, a value will be determined from the field attributes (if
            available) or the default values of -8, 64 will be used.
            Parameters are ignored is norm is not None.
        norm : Normalize or None, optional
            matplotlib Normalize instance used to scale luminance data. If not
            None the vmax and vmin parameters are ignored. If None, vmin and
            vmax are used for luminance scaling.
        cmap : str or None
            Matplotlib colormap name. None will use the default colormap for
            the field being plotted as specified by the Py-ART configuration.
        mask_outside : bool
            True to mask data outside of vmin, vmax. False performs no
            masking.
        title : str
            Title to label plot with, None to use default title generated from
            the field and lat,lon parameters. Parameter is ignored if
            title_flag is False.
        title_flag : bool
            True to add a title to the plot, False does not add a title.
        axislabels : (str, str)
            2-tuple of x-axis, y-axis labels. None for either label will use
            the default axis label. Parameter is ignored if axislabels_flag is
            False.
        axislabels_flag : bool
            True to add label the axes, False does not label the axes.
        colorbar_flag : bool
            True to add a colorbar with label to the axis. False leaves off
            the colorbar.
        colorbar_label : str
            Colorbar label, None will use a default label generated from the
            field information.
        colorbar_orient : 'vertical' or 'horizontal'
            Colorbar orientation.
        edges : bool
            True will interpolate and extrapolate the gate edges from the
            range, azimuth and elevations in the radar, treating these
            as specifying the center of each gate. False treats these
            coordinates themselved as the gate edges, resulting in a plot
            in which the last gate in each ray and the entire last ray are not
            not plotted.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.
        ticks : array
            Colorbar custom tick label locations.
        ticklabs : array
            Colorbar custom tick labels.

        """
        # parse parameters
        ax, fig = common.parse_ax_fig(ax, fig)
        vmin, vmax = common.parse_vmin_vmax(self.grid, field, vmin, vmax)
        cmap = common.parse_cmap(cmap, field)

        data = self.grid.fields[field]["data"][:, :, x_index]

        # mask the data where outside the limits
        if mask_outside:
            data = np.ma.masked_invalid(data)
            data = np.ma.masked_outside(data, vmin, vmax)

        # plot the grid
        y_1d = self.grid.y["data"] / 1000
        z_1d = self.grid.z["data"] / 1000

        if edges:
            if len(y_1d) > 1:
                y_1d = _interpolate_axes_edges(y_1d)
            if len(z_1d) > 1:
                z_1d = _interpolate_axes_edges(z_1d)
        xd, yd = np.meshgrid(y_1d, z_1d)

        if norm is not None:  # if norm is set do not override with vmin, vmax
            vmin = vmax = None

        pm = ax.pcolormesh(
            xd, yd, data, vmin=vmin, vmax=vmax, norm=norm, cmap=cmap, **kwargs
        )

        self.mappables.append(pm)
        self.fields.append(field)

        if title_flag:
            if title is None:
                ax.set_title(
                    common.generate_longitudinal_level_title(self.grid, field, x_index)
                )
            else:
                ax.set_title(title)

        if axislabels_flag:
            self._label_axes_longitude(axislabels, ax)

        if colorbar_flag:
            self.plot_colorbar(
                mappable=pm,
                label=colorbar_label,
                orientation=colorbar_orient,
                field=field,
                ax=ax,
                fig=fig,
                ticks=ticks,
                ticklabs=ticklabs,
            )

    def plot_cross_section(
        self,
        field,
        start,
        end,
        steps=100,
        interp_type="linear",
        x_axis=None,
        vmin=None,
        vmax=None,
        norm=None,
        cmap=None,
        title=None,
        title_flag=True,
        axislabels_flag=True,
        colorbar_flag=True,
        colorbar_label=None,
        colorbar_orient="vertical",
        ax=None,
        fig=None,
        ticks=None,
        ticklabs=None,
        **kwargs
    ):
        """
        Plot a cross section through a set of given points (latitude,
        longitude).

        This uses the MetPy cross section interpolation function.

        Additional arguments are passed to Matplotlib's pcolormesh function.

        Parameters
        ----------
        field : str
            Field to be plotted.
        start : tuple
            A latitude-longitude pair designating the start point of the cross
            section (units are degrees north and degrees east).
        end : tuple
            A latitude-longitude pair designating the end point of the cross
            section (units are degrees north and degrees east).
        steps: int
            The number of points along the geodesic between the start and the
            end point (including the end points) to use in the cross section.
            Defaults to 100.
        interp_type: str
            The interpolation method, either ‘linear’ or ‘nearest’
            (see xarray.DataArray.interp() for details). Defaults to ‘linear’.
        x_axis: str
            Field to use for plotting along the x-axis (ex. Latitude).
            Defaults to number of points from the first point.
        vmin, vmax : float
            Lower and upper range for the colormesh. If either parameter is
            None, a value will be determined from the field attributes (if
            available) or the default values of -8, 64 will be used.
            Parameters are ignored is norm is not None.
        norm : Normalize or None, optional
            matplotlib Normalize instance used to scale luminance data. If not
            None the vmax and vmin parameters are ignored. If None, vmin and
            vmax are used for luminance scaling.
        cmap : str or None
            Matplotlib colormap name. None will use the default colormap for
            the field being plotted as specified by the Py-ART configuration.
        mask_outside : bool
            True to mask data outside of vmin, vmax. False performs no
            masking.
        title : str
            Title to label plot with, None to use default title generated from
            the field and lat,lon parameters. Parameter is ignored if
            title_flag is False.
        title_flag : bool
            True to add a title to the plot, False does not add a title.
        axislabels : (str, str)
            2-tuple of x-axis, y-axis labels. None for either label will use
            the default axis label. Parameter is ignored if axislabels_flag is
            False.
        axislabels_flag : bool
            True to add label the axes, False does not label the axes.
        colorbar_flag : bool
            True to add a colorbar with label to the axis. False leaves off
            the colorbar.
        colorbar_label : str
            Colorbar label, None will use a default label generated from the
            field information.
        colorbar_orient : 'vertical' or 'horizontal'
            Colorbar orientation.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.
        ticks : array
            Colorbar custom tick label locations.
        ticklabs : array
            Colorbar custom tick labels.

        """
        if not _METPY_AVAILABLE:
            raise MissingOptionalDependency(
                "MetPy is required to use plot_cross_section but is not " + "installed!"
            )

        from metpy.interpolate import cross_section

        # parse parameters
        ax, fig = common.parse_ax_fig(ax, fig)
        vmin, vmax = common.parse_vmin_vmax(self.grid, field, vmin, vmax)
        cmap = common.parse_cmap(cmap, field)

        # Convert the grid into an xarray object
        ds = self.grid.to_xarray()

        # Extract the proj parameters
        proj_params = self.grid.get_projparams()

        # Convert the projection information into cartopy
        radar_crs = cartopy.crs.AzimuthalEquidistant(
            central_longitude=proj_params["lon_0"],
            central_latitude=proj_params["lat_0"],
        )

        # Now, convert that to cf-compliant coordinate information and assign
        # it to data
        projection_info = radar_crs.to_cf()
        ds = ds.metpy.assign_crs(projection_info)

        # Calculate the cross section, which returns a dataset
        ds = cross_section(ds, start, end, steps, interp_type).set_coords(
            ("lat", "lon")
        )

        # Convert from meters to km for the different variables
        ds["z"] = ds["z"] / 1000
        ds.z.attrs["units"] = "Distance above radar (km)"

        if x_axis == "y":
            ds["y"] = ds["y"] / 1000
            ds.y.attrs["units"] = "North South distance from radar (km)"

        if x_axis == "x":
            ds["x"] = ds["x"] / 1000
            ds.y.attrs["units"] = "East West distance from radar (km)"

        # Plot the data
        plot = ds[field].plot(
            y="z",
            x=x_axis,
            vmin=vmin,
            vmax=vmax,
            norm=norm,
            add_colorbar=False,
            ax=ax,
            cmap=cmap,
            **kwargs
        )

        self.mappables.append(plot)
        self.fields.append(field)

        if axislabels_flag:
            ax.set_ylabel(ds.z.attrs["units"])

        if title_flag:
            if title is None:
                ax.set_title(
                    common.generate_cross_section_title(self.grid, field, start, end)
                )
            else:
                ax.set_title(title)

        if colorbar_flag:
            self.plot_colorbar(
                mappable=plot,
                label=colorbar_label,
                orientation=colorbar_orient,
                field=field,
                ax=ax,
                fig=fig,
                ticks=ticks,
                ticklabs=ticklabs,
            )

    def plot_colorbar(
        self,
        mappable=None,
        orientation="horizontal",
        label=None,
        cax=None,
        ax=None,
        fig=None,
        field=None,
        ticks=None,
        ticklabs=None,
    ):
        """
        Plot a colorbar.

        Parameters
        ----------
        mappable : Image, ContourSet, etc.
            Image, ContourSet, etc to which the colorbar applied. If None the
            last mappable object will be used.
        field : str
            Field to label colorbar with.
        label : str
            Colorbar label. None will use a default value from the last field
            plotted.
        orient : str
            Colorbar orientation, either 'vertical' [default] or 'horizontal'.
        cax : Axis
            Axis onto which the colorbar will be drawn. None is also valid.
        ax : Axes
            Axis onto which the colorbar will be drawn. None is also valid.
        fig : Figure
            Figure to place colorbar on. None will use the current figure.
        ticks : array
            Colorbar custom tick label locations.
        ticklabs : array
            Colorbar custom tick labels.

        """
        if fig is None:
            fig = plt.gcf()

        if mappable is None:
            if len(self.mappables) == 0:
                raise ValueError("mappable must be specified.")
            mappable = self.mappables[-1]

        if label is None:
            if len(self.fields) == 0:
                raise ValueError("field must be specified.")

            field = self.grid.fields[self.fields[-1]]
            if "long_name" in field and "units" in field:
                label = field["long_name"] + "(" + field["units"] + ")"
            else:
                label = ""

        # plot the colorbar and set the label.
        cb = fig.colorbar(mappable, orientation=orientation, ax=ax, cax=cax)
        if ticks is not None:
            cb.set_ticks(ticks)
        if ticklabs is not None:
            cb.set_ticklabels(ticklabs)
        cb.set_label(label)

    def _find_nearest_grid_indices(self, lon, lat):
        """Find the nearest x, y grid indices for a given latitude and
        longitude."""

        # A similar method would make a good addition to the Grid class itself
        lon, lat = common.parse_lon_lat(self.grid, lon, lat)
        grid_lons, grid_lats = self.grid.get_point_longitude_latitude()
        diff = (grid_lats - lat) ** 2 + (grid_lons - lon) ** 2
        y_index, x_index = np.unravel_index(diff.argmin(), diff.shape)
        return x_index, y_index

    ##########################
    # Plot adjusting methods #
    ##########################

    def _get_label_x(self):
        """Get default label for x units."""

        return "East West distance from " + self.origin + " (km)"

    def _get_label_y(self):
        """Get default label for y units."""
        return "North South distance from " + self.origin + " (km)"

    def _get_label_z(self):
        """Get default label for z units."""
        return "Distance Above " + self.origin + " (km)"

    def _label_axes_grid(self, axis_labels, ax):
        """Set the x and y axis labels for a grid plot."""
        x_label, y_label = axis_labels
        if x_label is None:
            x_label = self._get_label_x()
        if y_label is None:
            y_label = self._get_label_y()
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

    def _label_axes_longitude(self, axis_labels, ax):
        """Set the x and y axis labels for a longitude slice."""
        x_label, y_label = axis_labels
        if x_label is None:
            x_label = self._get_label_y()
        if y_label is None:
            y_label = self._get_label_z()
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

    def _label_axes_latitude(self, axis_labels, ax):
        """Set the x and y axis labels for a latitude slice."""
        x_label, y_label = axis_labels
        if x_label is None:
            x_label = self._get_label_x()
        if y_label is None:
            y_label = self._get_label_z()
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

    ##########################
    # name generator methods #
    ##########################

    def generate_filename(self, field, level, ext="png"):
        """
        Generate a filename for a grid plot.

        Generated filename has form:
            grid_name_field_level_time.ext

        Parameters
        ----------
        field : str
            Field plotted.
        level : int
            Level plotted.
        ext : str
            Filename extension.

        Returns
        -------
        filename : str
            Filename suitable for saving a plot.

        """
        return common.generate_grid_filename(self.grid, field, level, ext)

    def generate_grid_title(self, field, level):
        """
        Generate a title for a plot.

        Parameters
        ----------
        field : str
            Field plotted.
        level : int
            Vertical level plotted.

        Returns
        -------
        title : str
            Plot title.

        """
        return common.generate_grid_title(self.grid, field, level)

    def generate_latitudinal_level_title(self, field, level):
        """
        Generate a title for a plot.

        Parameters
        ----------
        field : str
            Field plotted.
        level : int
            Latitudinal level plotted.

        Returns
        -------
        title : str
            Plot title.

        """
        return common.generate_latitudinal_level_title(self.grid, field, level)

    def generate_longitudinal_level_title(self, field, level):
        """
        Generate a title for a plot.

        Parameters
        ----------
        field : str
            Field plotted.
        level : int
            Longitudinal level plotted.

        Returns
        -------
        title : str
            Plot title.

        """
        return common.generate_longitudinal_level_title(self.grid, field, level)

    def cartopy_states(self):
        """Get state boundaries using cartopy."""
        return cartopy.feature.NaturalEarthFeature(
            category="cultural",
            name="admin_1_states_provinces_lines",
            scale="50m",
            facecolor="none",
        )

    def cartopy_political_boundaries(self):
        """Get political boundaries using cartopy."""
        return cartopy.feature.NaturalEarthFeature(
            category="cultural",
            name="admin_0_boundary_lines_land",
            scale="50m",
            facecolor="none",
        )

    def cartopy_coastlines(self):
        """Get coastlines using cartopy."""
        return cartopy.feature.NaturalEarthFeature(
            category="physical", name="coastline", scale="10m", facecolor="none"
        )


# These methods are a hack to allow gridlines when the projection is lambert
# https://nbviewer.jupyter.org/gist/ajdawson/dd536f786741e987ae4e


def find_side(ls, side):
    """
    Given a shapely LineString which is assumed to be rectangular, return the
    line corresponding to a given side of the rectangle.
    """
    minx, miny, maxx, maxy = ls.bounds
    points = {
        "left": [(minx, miny), (minx, maxy)],
        "right": [(maxx, miny), (maxx, maxy)],
        "bottom": [(minx, miny), (maxx, miny)],
        "top": [(minx, maxy), (maxx, maxy)],
    }
    return sgeom.LineString(points[side])


def lambert_xticks(ax, ticks):
    """Draw ticks on the bottom x-axis of a Lambert Conformal projection."""

    def te(xy):
        return xy[0]

    def lc(t, n, b):
        return np.vstack((np.zeros(n) + t, np.linspace(b[2], b[3], n))).T

    xticks, xticklabels = _lambert_ticks(ax, ticks, "bottom", lc, te)
    ax.xaxis.tick_bottom()
    ax.set_xticks(xticks)
    ax.set_xticklabels([ax.xaxis.get_major_formatter()(xtick) for xtick in xticklabels])


def lambert_yticks(ax, ticks):
    """Draw ticks on the left y-axis of a Lambert Conformal projection."""

    def te(xy):
        return xy[1]

    def lc(t, n, b):
        return np.vstack((np.linspace(b[0], b[1], n), np.zeros(n) + t)).T

    yticks, yticklabels = _lambert_ticks(ax, ticks, "left", lc, te)
    ax.yaxis.tick_left()
    ax.set_yticks(yticks)
    ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for ytick in yticklabels])


def _lambert_ticks(ax, ticks, tick_location, line_constructor, tick_extractor):
    """
    Get the tick locations and labels for a Lambert Conformal projection.
    """
    outline_patch = sgeom.LineString(ax.spines["geo"].get_path().vertices.tolist())
    axis = find_side(outline_patch, tick_location)
    n_steps = 30
    extent = ax.get_extent(cartopy.crs.PlateCarree())
    _ticks = []
    for t in ticks:
        xy = line_constructor(t, n_steps, extent)
        proj_xyz = ax.projection.transform_points(
            cartopy.crs.Geodetic(), xy[:, 0], xy[:, 1]
        )
        xyt = proj_xyz[..., :2]
        ls = sgeom.LineString(xyt.tolist())
        locs = axis.intersection(ls)
        if not locs:
            tick = [None]
        else:
            tick = tick_extractor(locs.xy)
        _ticks.append(tick[0])
    # Remove ticks that aren't visible:
    ticklabels = copy(ticks)
    while True:
        try:
            index = _ticks.index(None)
        except ValueError:
            break
        _ticks.pop(index)
        ticklabels = np.delete(ticklabels, index)
    return _ticks, ticklabels
