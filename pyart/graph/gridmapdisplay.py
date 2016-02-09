"""
pyart.graph.gridmapdisplay
==========================

A class for plotting grid objects with a basemap.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    GridMapDisplay

"""

from __future__ import print_function

import warnings

import numpy as np
import matplotlib.pyplot as plt
try:
    from mpl_toolkits.basemap import Basemap
    from mpl_toolkits.basemap import pyproj
    _BASEMAP_AVAILABLE = True
except ImportError:
    _BASEMAP_AVAILABLE = False

from . import common
from ..exceptions import MissingOptionalDependency, DeprecatedAttribute
from ..core.transforms import _interpolate_axes_edges


class GridMapDisplay(object):
    """
    A class for creating plots from a grid object on top of a Basemap.

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
        True to print debugging messages, False to supressed them.
    basemap : Basemap
        Last plotted basemap, None when no basemap has been plotted.
    mappables : list
        List of ContourSet, etc. which have been plotted, useful
        when adding colorbars.
    fields : list
        List of fields which have been plotted.

    """

    def __init__(self, grid, debug=False):
        """ initalize the object. """
        # check that basemap is available
        if not _BASEMAP_AVAILABLE:
            raise MissingOptionalDependency(
                "Basemap is required to use GridMapDisplay but is not " +
                "installed")

        # set attributes
        self.grid = grid
        self.debug = debug
        self.mappables = []
        self.fields = []
        self.origin = 'origin'
        self.basemap = None

    @property
    def proj(self):
        """ Deprecated proj attribute. """
        warnings.warn(
            "The 'proj' attribute has been deprecated and will be removed "
            "in future versions of Py-ART", category=DeprecatedAttribute)
        lat0 = self.grid.origin_latitude['data'][0]
        lon0 = self.grid.origin_longitude['data'][0]
        return pyproj.Proj(proj='aeqd', datum='NAD83', lat_0=lat0, lon_0=lon0)

    @property
    def grid_lats(self):
        """ Deprecated grid_lats attribute. """
        warnings.warn(
            "The 'grid_lats' attribute has been deprecated and will be "
            "removed in future versions of Py-ART",
            category=DeprecatedAttribute)
        return self.grid.point_latitude['data'][0]

    @property
    def grid_lons(self):
        """ Deprecated grid_lons attribute. """
        warnings.warn(
            "The 'grid_lons' attribute has been deprecated and will be "
            "removed in future versions of Py-ART",
            category=DeprecatedAttribute)
        return self.grid.point_latitude['data'][0]

    def plot_basemap(
            self, lat_lines=None, lon_lines=None, resolution='l',
            area_thresh=10000, auto_range=True, min_lon=-92, max_lon=-86,
            min_lat=40, max_lat=44, ax=None, **kwargs):
        """
        Plot a basemap.

        Parameters
        ----------
        lat_lines, lon_lines : array or None
            Locations at which to draw latitude and longitude lines.
            None will use default values which are resonable for maps of
            North America.
        auto_range : bool
            True to determine map ranges from the latitude and longitude
            limits of the grid. False will use the min_lon, max_lon, min_lat,
            and max_lat parameters for the map range.
        min_lat, max_lat, min_lon, max_lon : float
            Latitude and longitude ranges for the map projection region in
            degrees.  These parameter are not used if auto_range is True.
        resolution : 'c', 'l', 'i', 'h', or 'f'.
            Resolution of boundary database to use. See Basemap documentation
            for details.
        area_thresh : int
            Basemap area_thresh parameter. See Basemap documentation.
        ax : axes or None.
            Axis to add the basemap to, if None the current axis is used.
        kwargs: Basemap options
            Options to be passed to Basemap. If projection is not specified
            here it uses proj='merc' (mercator).

        """
        # make basemap
        self._make_basemap(resolution, area_thresh, auto_range,
                           min_lon, max_lon, min_lat, max_lat, ax, **kwargs)

        # parse the parameters
        if lat_lines is None:
            lat_lines = np.arange(30, 46, 1)
        if lon_lines is None:
            lon_lines = np.arange(-110, -75, 1)

        self.basemap.drawcoastlines(linewidth=1.25)
        self.basemap.drawstates()
        self.basemap.drawparallels(
            lat_lines, labels=[True, False, False, False])
        self.basemap.drawmeridians(
            lon_lines, labels=[False, False, False, True])

    def plot_grid(
            self, field, level=0,
            vmin=None, vmax=None, norm=None, cmap=None,
            mask_outside=False, title=None, title_flag=True,
            axislabels=(None, None), axislabels_flag=False,
            colorbar_flag=True, colorbar_label=None,
            colorbar_orient='vertical', edges=True,
            ax=None, fig=None, **kwargs):
        """
        Plot the grid onto the current basemap.

        Additional arguments are passed to Basemaps's pcolormesh function.

        Parameters
        ----------
        field : str
            Field to be plotted.
        level : int
            Index corresponding to the height level to be plotted.
        vmin, vmax : float
            Lower and upper range for the colormesh.  If either parameter is
            None, a value will be determined from the field attributes (if
            available) or the default values of -8, 64 will be used.
            Parameters are ignored is norm is not None.
        norm : Normalize or None, optional
            matplotlib Normalize instance used to scale luminance data.  If not
            None the vmax and vmin parameters are ignored.  If None, vmin and
            vmax are used for luminance scaling.
        cmap : str or None
            Matplotlib colormap name. None will use the default colormap for
            the field being plotted as specified by the Py-ART configuration.
        mask_outside : bool
            True to mask data outside of vmin, vmax.  False performs no
            masking.
        title : str
            Title to label plot with, None to use default title generated from
            the field and level parameters. Parameter is ignored if title_flag
            is False.
        title_flag : bool
            True to add a title to the plot, False does not add a title.
        axislabels : (str, str)
            2-tuple of x-axis, y-axis labels.  None for either label will use
            the default axis label.  Parameter is ignored if axislabels_flag is
            False.
        axislabel_flag : bool
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
            as specifying the center of each gate.  False treats these
            coordinates themselved as the gate edges, resulting in a plot
            in which the last gate in each ray and the entire last ray are not
            not plotted.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.

        """
        # parse parameters
        ax, fig = common.parse_ax_fig(ax, fig)
        norm, vmin, vmax = common.parse_norm_vmin_vmax(
            norm, self.grid, field, vmin, vmax)
        cmap = common.parse_cmap(cmap, field)

        basemap = self.get_basemap()

        data = self.grid.fields[field]['data'][level]

        # mask the data where outside the limits
        if mask_outside:
            data = np.ma.masked_invalid(data)
            data = np.ma.masked_outside(data, vmin, vmax)

        # plot the grid
        lons, lats = self.grid.get_point_longitude_latitude(edges=edges)
        pm = basemap.pcolormesh(
            lons, lats, data, vmin=vmin, vmax=vmax, cmap=cmap, norm=norm,
            latlon=True, **kwargs)
        self.mappables.append(pm)
        self.fields.append(field)

        if title_flag:
            if title is None:
                ax.set_title(self.generate_grid_title(field, level))
            else:
                ax.set_title(title)

        if axislabels_flag:
            self._label_axes_grid(axislabels, ax)

        if colorbar_flag:
            self.plot_colorbar(
                mappable=pm, label=colorbar_label, orientation=colorbar_orient,
                field=field, ax=ax, fig=fig)

        return

    def plot_crosshairs(
            self, lon=None, lat=None, line_style='r--', linewidth=2, ax=None):
        """
        Plot crosshairs at a given longitude and latitude.

        Parameters
        ----------
        lon, lat : float
            Longitude and latitude (in degrees) where the crosshairs should
            be placed.  If None the center of the grid is used.
        line_style : str
            Matplotlib string describing the line style.
        linewidth : float
            Width of markers in points.
        ax : axes or None.
            Axis to add the crosshairs to, if None the current axis is used.

        """
        # parse the parameters
        ax = common.parse_ax(ax)
        lon, lat = common.parse_lon_lat(self.grid, lon, lat)
        basemap = self.get_basemap()

        # add crosshairs.
        x_lon, y_lon = basemap(
            np.array([lon, lon]),
            np.array([basemap.latmin, basemap.latmax]))
        x_lat, y_lat = basemap(
            np.array([basemap.lonmin, basemap.lonmax]),
            np.array([lat, lat]))
        ax.plot(x_lon, y_lon, line_style, linewidth=linewidth)
        ax.plot(x_lat, y_lat, line_style, linewidth=linewidth)
        return

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
            Longitude and latitude (in degrees) specifying the slice.  If
            None the center of the grid is used.

        """
        # parse parameters
        _, y_index = self._find_nearest_grid_indices(lon, lat)
        self.plot_latitudinal_level(field=field, y_index=y_index, **kwargs)

    def plot_latitudinal_level(
            self, field, y_index,
            vmin=None, vmax=None, norm=None, cmap=None,
            mask_outside=False, title=None, title_flag=True,
            axislabels=(None, None), axislabels_flag=True, colorbar_flag=True,
            colorbar_label=None, colorbar_orient='vertical', edges=True,
            ax=None, fig=None, **kwargs):
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
            Lower and upper range for the colormesh.  If either parameter is
            None, a value will be determined from the field attributes (if
            available) or the default values of -8, 64 will be used.
            Parameters are ignored is norm is not None.
        norm : Normalize or None, optional
            matplotlib Normalize instance used to scale luminance data.  If not
            None the vmax and vmin parameters are ignored.  If None, vmin and
            vmax are used for luminance scaling.
        cmap : str or None
            Matplotlib colormap name. None will use the default colormap for
            the field being plotted as specified by the Py-ART configuration.
        mask_outside : bool
            True to mask data outside of vmin, vmax.  False performs no
            masking.
        title : str
            Title to label plot with, None to use default title generated from
            the field and lat,lon parameters. Parameter is ignored if
            title_flag is False.
        title_flag : bool
            True to add a title to the plot, False does not add a title.
        axislabels : (str, str)
            2-tuple of x-axis, y-axis labels.  None for either label will use
            the default axis label.  Parameter is ignored if axislabels_flag is
            False.
        axislabel_flag : bool
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
            as specifying the center of each gate.  False treats these
            coordinates themselved as the gate edges, resulting in a plot
            in which the last gate in each ray and the entire last ray are not
            not plotted.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.

        """
        # parse parameters
        ax, fig = common.parse_ax_fig(ax, fig)
        norm, vmin, vmax = common.parse_norm_vmin_vmax(
            norm, self.grid, field, vmin, vmax)
        cmap = common.parse_cmap(cmap, field)

        data = self.grid.fields[field]['data'][:, y_index, :]

        # mask the data where outside the limits
        if mask_outside:
            data = np.ma.masked_invalid(data)
            data = np.ma.masked_outside(data, vmin, vmax)

        # plot the grid
        x_1d = self.grid.x['data'] / 1000.
        z_1d = self.grid.z['data'] / 1000.
        if edges:
            if len(x_1d) > 1:
                x_1d = _interpolate_axes_edges(x_1d)
            if len(z_1d) > 1:
                z_1d = _interpolate_axes_edges(z_1d)
        xd, yd = np.meshgrid(x_1d, z_1d)
        pm = ax.pcolormesh(
            xd, yd, data, vmin=vmin, vmax=vmax, norm=norm, cmap=cmap, **kwargs)
        self.mappables.append(pm)
        self.fields.append(field)

        if title_flag:
            if title is None:
                ax.set_title(common.generate_latitudinal_level_title(
                    self.grid, field, y_index))
            else:
                ax.set_title(title)

        if axislabels_flag:
            self._label_axes_latitude(axislabels, ax)

        if colorbar_flag:
            self.plot_colorbar(
                mappable=pm, label=colorbar_label, orientation=colorbar_orient,
                field=field, ax=ax, fig=fig)
        return

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
        x_index, _ = self._find_nearest_grid_indices(lon, lat)
        self.plot_longitudinal_level(field=field, x_index=x_index, **kwargs)

    def plot_longitudinal_level(
            self, field, x_index,
            vmin=None, vmax=None, norm=None, cmap=None,
            mask_outside=False, title=None, title_flag=True,
            axislabels=(None, None), axislabels_flag=True, colorbar_flag=True,
            colorbar_label=None, colorbar_orient='vertical', edges=True,
            ax=None, fig=None, **kwargs):
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
            Lower and upper range for the colormesh.  If either parameter is
            None, a value will be determined from the field attributes (if
            available) or the default values of -8, 64 will be used.
            Parameters are ignored is norm is not None.
        norm : Normalize or None, optional
            matplotlib Normalize instance used to scale luminance data.  If not
            None the vmax and vmin parameters are ignored.  If None, vmin and
            vmax are used for luminance scaling.
        cmap : str or None
            Matplotlib colormap name. None will use the default colormap for
            the field being plotted as specified by the Py-ART configuration.
        mask_outside : bool
            True to mask data outside of vmin, vmax.  False performs no
            masking.
        title : str
            Title to label plot with, None to use default title generated from
            the field and lat,lon parameters. Parameter is ignored if
            title_flag is False.
        title_flag : bool
            True to add a title to the plot, False does not add a title.
        axislabels : (str, str)
            2-tuple of x-axis, y-axis labels.  None for either label will use
            the default axis label.  Parameter is ignored if axislabels_flag is
            False.
        axislabel_flag : bool
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
            as specifying the center of each gate.  False treats these
            coordinates themselved as the gate edges, resulting in a plot
            in which the last gate in each ray and the entire last ray are not
            not plotted.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.

        """
        # parse parameters
        ax, fig = common.parse_ax_fig(ax, fig)
        norm, vmin, vmax = common.parse_norm_vmin_vmax(
            norm, self.grid, field, vmin, vmax)
        cmap = common.parse_cmap(cmap, field)

        data = self.grid.fields[field]['data'][:, :, x_index]

        # mask the data where outside the limits
        if mask_outside:
            data = np.ma.masked_invalid(data)
            data = np.ma.masked_outside(data, vmin, vmax)

        # plot the grid
        y_1d = self.grid.y['data'] / 1000.
        z_1d = self.grid.z['data'] / 1000.
        if edges:
            if len(y_1d) > 1:
                y_1d = _interpolate_axes_edges(y_1d)
            if len(z_1d) > 1:
                z_1d = _interpolate_axes_edges(z_1d)
        xd, yd = np.meshgrid(y_1d, z_1d)
        pm = ax.pcolormesh(
            xd, yd, data, vmin=vmin, vmax=vmax, cmap=cmap, norm=norm, **kwargs)
        self.mappables.append(pm)
        self.fields.append(field)

        if title_flag:
            if title is None:
                ax.set_title(common.generate_longitudinal_level_title(
                    self.grid, field, x_index))
            else:
                ax.set_title(title)

        if axislabels_flag:
            self._label_axes_longitude(axislabels, ax)

        if colorbar_flag:
            self.plot_colorbar(
                mappable=pm, label=colorbar_label, orientation=colorbar_orient,
                field=field, ax=ax, fig=fig)
        return

    def plot_colorbar(
            self, mappable=None, orientation='horizontal', label=None,
            cax=None, ax=None, fig=None, field=None):
        """
        Plot a colorbar.

        Parameters
        ----------
        mappable : Image, ContourSet, etc.
            Image, ContourSet, etc to which the colorbar applied.  If None the
            last mappable object will be used.
        field : str
            Field to label colorbar with.
        label : str
            Colorbar label.  None will use a default value from the last field
            plotted.
        orient : str
            Colorbar orientation, either 'vertical' [default] or 'horizontal'.
        cax : Axis
            Axis onto which the colorbar will be drawn.  None is also valid.
        ax : Axes
            Axis onto which the colorbar will be drawn. None is also valid.
        fig : Figure
            Figure to place colorbar on.  None will use the current figure.

        """
        if fig is None:
            fig = plt.gcf()

        if mappable is None:
            if len(self.mappables) == 0:
                raise ValueError('mappable must be specified.')
            else:
                mappable = self.mappables[-1]

        if label is None:
            if len(self.fields) == 0:
                raise ValueError('field must be specified.')
            field = self.grid.fields[self.fields[-1]]
            if 'long_name' in field and 'units' in field:
                label = field['long_name'] + '(' + field['units'] + ')'
            else:
                label = ''

        # plot the colorbar and set the label.
        cb = fig.colorbar(mappable, orientation=orientation, ax=ax, cax=cax)
        cb.set_label(label)

        return

    def _make_basemap(
            self, resolution='l', area_thresh=10000, auto_range=True,
            min_lon=-92, max_lon=-86, min_lat=40, max_lat=44, ax=None,
            **kwargs):
        """
        Make a basemap.

        Parameters
        ----------
        auto_range : bool
            True to determine map ranges from the latitude and longitude limits
            of the grid. False will use the min_lon, max_lon, min_lat, and
            max_lat parameters for the map range.
        min_lat, max_lat, min_lon, max_lon : float
            Latitude and longitude ranges for the map projection region in
            degrees.  These parameter are not used if auto_range is True.
        resolution : 'c', 'l', 'i', 'h', or 'f'.
            Resolution of boundary database to use. See Basemap documentation
            for details.
        area_thresh : int
            Basemap area_thresh parameter. See Basemap documentation.
        ax : axes or None.
            Axis to add the basemap to, if None the current axis is used.
        kwargs: Basemap options
            Options to be passed to Basemap. If projection is not specified
            here it uses proj='merc' (mercator).
        """
        # parse the parameters
        ax = common.parse_ax(ax)

        # determine map region
        if auto_range:
            max_lat = self.grid.point_latitude['data'][0].max()
            max_lon = self.grid.point_longitude['data'][0].max()
            min_lat = self.grid.point_latitude['data'][0].min()
            min_lon = self.grid.point_longitude['data'][0].min()

        if self.debug:
            print("Maximum latitude: ", max_lat)
            print("Maximum longitude: ", max_lon)
            print("Minimum latitude: ", min_lat)
            print("Minimum longitute: ", min_lon)

        # determine plot center
        lat_0 = self.grid.origin_latitude['data'][0]
        lon_0 = self.grid.origin_longitude['data'][0]

        default_args = {
            'lat_0': lat_0, 'lon_0': lon_0, 'lat_ts': lat_0,
            'projection': 'merc', 'area_thresh': area_thresh,
            'resolution': resolution, 'ax': ax}

        using_corners = (None not in [min_lon, min_lat, max_lon, max_lat])
        if using_corners:
            default_args['llcrnrlon'] = min_lon
            default_args['llcrnrlat'] = min_lat
            default_args['urcrnrlon'] = max_lon
            default_args['urcrnrlat'] = max_lat
        else:
            # determine width and height of the plot
            x = self.grid.x['data'][0]
            y = self.grid.y['data'][0]
            default_args['width'] = (x.max() - x.min())
            default_args['height'] = (y.max() - y.min())

        for key in default_args.keys():
            if key not in kwargs:
                kwargs[key] = default_args[key]

        # plot the basemap
        self.basemap = Basemap(**kwargs)

        return self.basemap

    def _find_nearest_grid_indices(self, lon, lat):
        """
        Find the nearest x, y grid indices for a given latitude and longitude.
        """
        # A similar method would make a good addition to the Grid class itself
        lon, lat = common.parse_lon_lat(self.grid, lon, lat)
        grid_lons, grid_lats = self.grid.get_point_longitude_latitude()
        diff = (grid_lats - lat)**2 + (grid_lons - lon)**2
        y_index, x_index = np.unravel_index(diff.argmin(), diff.shape)
        return x_index, y_index

    ##########################
    # Plot adjusting methods #
    ##########################

    def _get_label_x(self):
        """ Get default label for x units. """
        return 'East West distance from ' + self.origin + ' (km)'

    def _get_label_y(self):
        """ Get default label for y units. """
        return 'North South distance from ' + self.origin + ' (km)'

    def _get_label_z(self):
        """ Get default label for z units. """
        return 'Distance Above ' + self.origin + '  (km)'

    def _label_axes_grid(self, axis_labels, ax):
        """ Set the x and y axis labels for a grid plot. """
        x_label, y_label = axis_labels
        if x_label is None:
            x_label = self._get_label_x()
        if y_label is None:
            y_label = self._get_label_y()
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

    def _label_axes_longitude(self, axis_labels, ax):
        """ Set the x and y axis labels for a longitude slice. """
        x_label, y_label = axis_labels
        if x_label is None:
            x_label = self._get_label_y()
        if y_label is None:
            y_label = self._get_label_z()
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

    def _label_axes_latitude(self, axis_labels, ax):
        """ Set the x and y axis labels for a latitude slice. """
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

    def generate_filename(self, field, level, ext='png'):
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
            Verical level plotted.

        Returns
        -------
        title : str
            Plot title.

        """
        return common.generate_grid_title(self.grid, field, level)

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
        return common.generate_longitudinal_level_title(
            self.grid, field, level)

    def generate_latitudinal_level_title(self, field, level):
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
        return common.generate_latitudinal_level_title(
            self.grid, field, level)

    ##########################
    #      get methods       #
    ##########################

    def get_basemap(self):
        """ get basemap of the plot """
        if self.basemap is None:
            self._make_basemap()

        return self.basemap
