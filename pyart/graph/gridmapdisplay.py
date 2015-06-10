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

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import pyproj


class GridMapDisplay():
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
    proj : Proj
        Object for performing cartographic transformations specific to the
        grid.
    grid_lons : array
        Grid longitudes in degrees.
    grid_lats : array
        Grid latitudes in degress.
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
        self.grid = grid
        self.debug = debug

        # set up the projection
        lat0 = grid.axes['lat']['data'][0]
        lon0 = grid.axes['lon']['data'][0]
        self.proj = pyproj.Proj(proj='lcc', datum='NAD83',
                                lat_0=lat0, lon_0=lon0,
                                x_0=0.0, y_0=0.0)

       # determine grid latitudes and longitudes.
        x_1d = grid.axes['x_disp']['data']
        y_1d = grid.axes['y_disp']['data']
        x_2d, y_2d = np.meshgrid(x_1d, y_1d)
        self.grid_lons, self.grid_lats = self.proj(x_2d, y_2d, inverse=True)

        # set attributes
        self.mappables = []
        self.fields = []
        self.basemap = None

    def plot_basemap(self, lat_lines=None, lon_lines=None,
                     resolution='l', area_thresh=10000,
                     auto_range=True,
                     min_lon=-92, max_lon=-86, min_lat=40, max_lat=44,
                     ax=None):
        """
        Plot a basemap.

        Parameters
        ----------
        lat_lines, lon_lines : array or None
            Locations at which to draw latitude and longitude lines.
            None will use default values which are resonable for maps of
            North America.
        auto_range : bool
            True to determine map ranges from the grid_lats and grid_lons
            attribute.  False will use the min_lon, max_lon, min_lat, and
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

        """
        # parse the parameters
        ax = self._parse_ax(ax)
        if lat_lines is None:
            lat_lines = np.arange(30, 46, 1)
        if lon_lines is None:
            lon_lines = np.arange(-110, -75, 1)

        # determine map region
        if auto_range:
            max_lat = self.grid_lats.max()
            max_lon = self.grid_lons.max()
            min_lat = self.grid_lats.min()
            min_lon = self.grid_lons.min()

        if self.debug:
            print("Maximum latitude: ", max_lat)
            print("Maximum longitude: ", max_lon)
            print("Minimum latitude: ", min_lat)
            print("Minimum longitute: ", min_lon)

        # plot the basemap
        basemap = Basemap(llcrnrlon=min_lon, llcrnrlat=min_lat,
                          urcrnrlon=max_lon, urcrnrlat=max_lat,
                          projection='mill', area_thresh=area_thresh,
                          resolution=resolution, ax=ax)
        basemap.drawcoastlines(linewidth=1.25)
        basemap.drawstates()
        basemap.drawparallels(lat_lines, labels=[True, False, False, False])
        basemap.drawmeridians(lon_lines, labels=[False, False, False, True])
        self.basemap = basemap

    def plot_grid(self, field, level=0, vmin=None, vmax=None, cmap='jet'):
        """
        Plot the grid onto the current basemap.

        Parameters
        ----------
        field : str
            Field to be plotted.
        level : int
            Index of the z level to be plotted.
        vmin, vmax : float
            Lower and upper range for the colormesh.  If either parameter is
            None, a value will be determined from the field attributes (if
            available) or the default values of -8, 64 will be used.
        cmap : str
            Matplotlib colormap name or colormap object.

        """
        # parse parameters
        vmin, vmax = self._parse_vmin_vmax(field, vmin, vmax)
        self._check_basemap()

        # plot the grid
        xd, yd = self.basemap(self.grid_lons, self.grid_lats)
        self.mappables.append(self.basemap.pcolormesh(
            xd, yd, self.grid.fields[field]['data'][level],
            vmin=vmin, vmax=vmax, cmap=cmap))
        self.fields.append(field)
        return

    def plot_crosshairs(self, lon=None, lat=None,
                        line_style='r--', linewidth=2, ax=None):
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
        ax = self._parse_ax(ax)
        lon, lat = self._parse_lon_lat(lon, lat)
        self._check_basemap()

        # add crosshairs.
        x_lon, y_lon = self.basemap(
            np.array([lon, lon]),
            np.array([self.basemap.latmin, self.basemap.latmax]))
        x_lat, y_lat = self.basemap(
            np.array([self.basemap.lonmin, self.basemap.lonmax]),
            np.array([lat, lat]))
        ax.plot(x_lon, y_lon, line_style, linewidth=linewidth)
        ax.plot(x_lat, y_lat, line_style, linewidth=linewidth)
        return

    def plot_latitude_slice(self, field, lon=None, lat=None,
                            vmin=None, vmax=None, cmap='jet', ax=None):
        """
        Plot a slice along a given latitude.

        Parameters
        ----------
        field : str
            Field to be plotted.
        lon, lat : float
            Longitude and latitude (in degrees) specifying the slice.  If
            None the center of the grid is used.
        vmin, vmax : float
            Lower and upper range for the colormesh.  If either parameter is
            None, a value will be determined from the field attributes (if
            available) or the default values of -8, 64 will be used.
        cmap : str
            Matplotlib colormap name or colormap object.
        ax : axes or None.
            Where to create the plot, if None the current axis is used.


        """
        # parse the parameters
        ax = self._parse_ax(ax)
        lon, lat = self._parse_lon_lat(lon, lat)
        vmin, vmax = self._parse_vmin_vmax(field, vmin, vmax)

        # plot the slice.
        x_index, y_index = self._find_nearest_grid_indices(lon, lat)
        self.mappables.append(ax.pcolormesh(
            self.grid.axes['x_disp']['data'] / 1000.0,
            self.grid.axes['z_disp']['data'] / 1000.0,
            self.grid.fields[field]['data'][:, y_index, :],
            vmin=vmin, vmax=vmax, cmap=cmap))
        self.fields.append(field)
        ax.set_ylabel('Height (km)')
        ax.set_title('Slice at ' + str(lat) + ' Latitude')
        return

    def plot_longitude_slice(self, field, lon=None, lat=None,
                             vmin=None, vmax=None, cmap='jet', ax=None):
        """
        Plot a slice along a given longitude.

        Parameters
        ----------
        field : str
            Field to be plotted.
        lon, lat : float
            Longitude and latitude (in degrees) specifying the slice.  If
            None the center of the grid is used.
        vmin, vmax : float
            Lower and upper range for the colormesh.  If either parameter is
            None, a value will be determined from the field attributes (if
            available) or the default values of -8, 64 will be used.
        cmap : str
            Matplotlib colormap name or colormap object.
        ax : axes or None.
            Where to create the plot, if None the current axis is used.

        """
        # parse the parameters
        ax = self._parse_ax(ax)
        lon, lat = self._parse_lon_lat(lon, lat)
        vmin, vmax = self._parse_vmin_vmax(field, vmin, vmax)

        # plot the slice
        x_index, y_index = self._find_nearest_grid_indices(lon, lat)
        self.mappables.append(ax.pcolormesh(
            self.grid.axes['y_disp']['data'] / 1000.0,
            self.grid.axes['z_disp']['data'] / 1000.0,
            self.grid.fields[field]['data'][:, :, x_index],
            vmin=vmin, vmax=vmax, cmap=cmap))
        self.fields.append(field)
        ax.set_ylabel('Height (km)')
        ax.set_title('Slice at ' + str(lon) + ' Longitude')
        return

    def plot_colorbar(self, mappable=None, orientation='horizontal',
                      label=None, cax=None):
        """
        Plot a colorbar to an axis.

        Parameters
        ----------
        mappable : ContourSet, etc.
            Matplotlib object to create colorbar for.  None will use
            the last mappable plotted.
        orientation : 'vertical' or 'horizontal'
            Colorbar orientation.
        label : str
            Colorbar label.  None will use: 'long_name' ('units') for
            the last field plotted or '' if the field does not have these
            keys.
        cax : axis
            Axes object into which the colorbar will be drawn.
            None is also allowed.

        """
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
        cb = plt.colorbar(cax=cax, mappable=mappable, orientation=orientation)
        cb.set_label(label)

        return

    def _check_basemap(self):
        """ Check that basemap is not None, raise ValueError if it is. """
        if self.basemap is None:
            raise ValueError('no basemap plotted')

    def _find_nearest_grid_indices(self, lon, lat):
        """
        Find the nearest x, y grid indices for a given latitude and longitude.
        """
        lon, lat = self._parse_lon_lat(lon, lat)
        x_cut, y_cut = self.proj(lon, lat)

        if self.debug:
            print("x_cut: ", x_cut)
            print("y_cut: ", y_cut)

        x_index = np.abs(self.grid.axes['x_disp']['data'] - x_cut).argmin()
        y_index = np.abs(self.grid.axes['y_disp']['data'] - y_cut).argmin()

        if self.debug:
            print("x_index", x_index)
            print("y_index", y_index)

        return x_index, y_index

    def _parse_vmin_vmax(self, field, vmin, vmax):
        """ Parse vmin and vmax parameters. """
        if vmin is None:
            if 'valid_min' in self.grid.fields[field]:
                vmin = self.grid.fields[field]['valid_min']
            else:
                vmin = -8

        if vmax is None:
            if 'valid_max' in self.grid.fields[field]:
                vmax = self.grid.fields[field]['valid_max']
            else:
                vmax = 64
        return vmin, vmax

    def _parse_lon_lat(self, lon, lat):
        """ Parse lat and lon parameters """
        if lat is None:
            lat = self.grid.axes['lat']['data'][0]
        if lon is None:
            lon = self.grid.axes['lon']['data'][0]
        return lon, lat

    def _parse_ax(self, ax):
        """ Parse the ax parameter. """
        if ax is None:
            return plt.gca()
        return ax
