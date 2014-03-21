"""
pyart.graph.gridmapdisplay
==========================

A class for plotting grid objects with a basemap.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    GridMapDisplay

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import pyproj

from ..util import datetime_utils


class GridMapDisplay:
    """
    A class for creating plots from a grid object on top of a Basemap.

    Parameters
    ----------
    grid : Grid
        Grid with data which will be used to create plots.
    proj : str
        The projection to be used.
    datum : str
    debug : bool
        True to print debugging information.

    Attributes
    ----------
    fields : dict
        Field data and metadata dictionaries.
    proj : Proj
        Object for performing cartographic transformations specific to the
        grid.
    grid_lons : array
        Grid longitudes in degrees.
    grid_lats : array
        Grid latitudes in degrees.
    start_time : datetime
        Start time of gridded data.
    x, y, z : np.ndarray
        Eastward, Northward, and vertical Cartesian axes in meters.
    lat, lon : float
        Latitude and longitude location of the grid origin.
    source : str
        The source of the gridded data. For gridded radar data, this should
        be the instrument name of the radar.
    basemap : Basemap
        Last plotted Basemap. None when no Basemap has been plotted.
    plots : list
        List of plots created.
    plot_vars : list
        List of fields plotted, order matches plot list.
    debug : bool
        True to print debugging information.

    """

    def __init__(self, grid, proj='lcc', datum='NAD83', debug=False):
        """ Initialize graph object. """
        
        self.debug = debug
        
        # Populate the fields attribute
        self.fields = grid.fields
        
        # Set source or instrument attribute
        if 'radar_0_instrument_name' in grid.metadata:
            self.source = str(grid.metadata['radar_0_instrument_name'])
        else:
            self.source = 'Model'
        
        # Populate axes attributes
        self.start_time = datetime_utils.datetime_from_grid(grid)
        self.x = grid.axes['x_disp']['data']
        self.y = grid.axes['y_disp']['data']
        self.z = grid.axes['z_disp']['data']
        self.lat = grid.axes['lat']['data'][0]
        self.lon = grid.axes['lon']['data'][0]

        # Set up the projection
        self.proj = pyproj.Proj(proj=proj, datum=datum, lat_0=self.lat,
                                lon_0=self.lon, x_0=0.0, y_0=0.0)
        self.projection = proj
        self.datum = datum

       # Determine grid latitudes and longitudes
        X, Y = np.meshgrid(self.x, self.y)
        self.grid_lons, self.grid_lats = self.proj(X, Y, inverse=True)

        # Populate graph attributes
        self.plots = []
        self.plot_vars = []
        self.basemap = None

    def create_basemap(self, lat_lines=None, lon_lines=None, resolution='l',
                       area_thresh=10000, auto_range=True, min_lon=-92,
                       max_lon=-86, min_lat=40, max_lat=44, ax=None):
        """
        Create a Basemap instance.

        Parameters
        ----------
        lat_lines, lon_lines : array_like
            Locations at which to draw latitude and longitude lines.
            None will use default values which are reasonable for maps of
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
            Basemap area_thresh parameter. See Basemap documentation for more
            details.
        ax : Axes
            Axis to add the Basemap to. If None the current axis is used.

        """
        
        # Parse the parameters
        ax = self._parse_ax(ax)
        if lat_lines is None:
            lat_lines = np.arange(30, 46, 1)
        if lon_lines is None:
            lon_lines = np.arange(-110, -75, 1)

        # Determine the map region
        if auto_range:
            max_lat = self.grid_lats.max()
            max_lon = self.grid_lons.max()
            min_lat = self.grid_lats.min()
            min_lon = self.grid_lons.min()

        if self.debug:
            print "Maximum latitude is %.2f " %max_lat
            print "Maximum longitude is %.2f " %max_lon
            print "Minimum latitude is %.2f " %min_lat
            print "Minimum longitude is %.2f " %min_lon

        # Create the Basemap instance
        m = Basemap(llcrnrlon=min_lon, llcrnrlat=min_lat, urcrnrlon=max_lon,
                    urcrnrlat=max_lat, lat_0=self.lat, lon_0=self.lon,
                    area_thresh=area_thresh, resolution=resolution,
                    projection=self.projection, ax=ax)
        m.drawcoastlines(linewidth=1.25)
        m.drawstates()
        m.drawparallels(lat_lines, labels=[True, False, False, False])
        m.drawmeridians(lon_lines, labels=[False, False, False, True])
        
        self.basemap = m

    def plot_basemap(self, field, level=0, vmin=None, vmax=None, cmap='jet',
                     title=None, title_flag=True, colorbar_label=None,
                     colorbar_flag=True, orientation='vertical', fig=None,
                     ax=None):
        """
        Plot the grid onto the current Basemap.

        Parameters
        ----------
        field : str
            Field to be plotted.
        level : int
            Index of the height level to be plotted.
        vmin, vmax : float
            Lower and upper range for the color bar.  If either parameter is
            None, a value will be determined from the field attributes (if
            available) or the default reflectivity values will be used.
        cmap : str or Colormap
            Color map name or Colormap instance.

        """
        
        self._check_basemap()
        
        # Parse parameters
        ax = self._parse_ax(ax)
        vmin, vmax = self._parse_vmin_vmax(field, vmin, vmax)
        x, y = self.basemap(self.grid_lons, self.grid_lats)

        # Plot the grid
        qm = self.basemap.pcolormesh(x, y, self.fields[field]['data'][level],
                                     vmin=vmin, vmax=vmax, cmap=cmap)
        
        self.plots.append(qm)
        self.plot_vars.append(field)
        
        # Plot color bar
        if colorbar_flag:
            self.plot_colorbar(field, mappable=qm, orientation=orientation,
                               label=colorbar_label, fig=fig, ax=ax)
        
        # Set title
        if title_flag:
            self._set_xy_title(field, level, title, ax)
        
        return
    
    def plot_xy_slice(self, field, level=0, vmin=None, vmax=None, cmap='jet',
                      title=None, title_flag=True, colorbar_label=None,
                      colorbar_flag=True, orientation='vertical',
                      axislabels=(None,None), axislabels_flag=True,
                      fig=None, ax=None):
        """
        Plot the grid on a horizontal cross section plot.
        
        Parameters
        ----------
        
        """
        
        # Parse the parameters
        ax = self._parse_ax(ax)
        vmin, vmax = self._parse_vmin_vmax(field, vmin, vmax)
        
        qm = ax.pcolormesh(self.x/1000, self.y/1000,
                           self.fields[field]['data'][level],
                           vmin=vmin, vmax=vmax, cmap=cmap)
        
        self.plots.append(qm)
        self.plot_vars.append(field)
        
        # Set axis labels
        if axislabels_flag:
            self._label_axes_xy(axislabels, ax)
            
        # Plot color bar
        if colorbar_flag:
            self.plot_colorbar(field, mappable=qm, orientation=orientation,
                               label=colorbar_label, fig=fig, ax=ax)
            
        # Set title
        if title_flag:
            self._set_xy_title(field, level, title, ax)
            
        return

    def plot_latitude_slice(self, field, lon=None, lat=None, vmin=None,
                            vmax=None, cmap='jet', title=None,
                            title_flag=True, colorbar_label=None,
                            colorbar_flag=True, orientation='vertical',
                            fig=None, ax=None):
        """
        Plot a slice along a given latitude.

        Parameters
        ----------
        field : str
            Field to be plotted.
        lon, lat : float
            Longitude and latitude in degrees specifying the slice. If None
            then the center of the grid is used.
        vmin, vmax : float
            Lower and upper range for the color bar.  If either parameter is
            None, a value will be determined from the field attributes (if
            available) or the default reflectivity values will be used.
        cmap : str or Colormap
            Color map name or Colormap instance.
        ax : Axes
            Axis to add the plot to. If None the current axis is used.
        """
        
        # Parse the parameters
        ax = self._parse_ax(ax)
        lon, lat = self._parse_lon_lat(lon, lat)
        vmin, vmax = self._parse_vmin_vmax(field, vmin, vmax)
        i, j = self._find_nearest_grid_indices(lon, lat)
        
        qm = ax.pcolormesh(self.x/1000, self.z/1000,
                           self.fields[field]['data'][:,j,:],
                           vmin=vmin, vmax=vmax, cmap=cmap)
            
        self.plots.append(qm)
        self.plot_vars.append(field)
        
        # Plot color bar
        if colorbar_flag:
            self.plot_colorbar(field, mappable=qm, orientation=orientation,
                               label=colorbar_label, fig=fig, ax=ax)
        
        # Set title
        if title_flag:
            self._set_lat_slice_title(field, lat, title, ax)
        
        return

    def plot_longitude_slice(self, field, lon=None, lat=None, vmin=None,
                             vmax=None, cmap='jet', title=None,
                             title_flag=True, colorbar_label=None,
                             colorbar_flag=True, orientation='vertical',
                             fig=None, ax=None):
        """
        Plot a slice along a given longitude.

        Parameters
        ----------
        field : str
            Field to be plotted.
        lon, lat : float
            Longitude and latitude in degrees specifying the slice. If None
            then the center of the grid is used.
        vmin, vmax : float
            Lower and upper range for the color bar.  If either parameter is
            None, a value will be determined from the field attributes (if
            available) or the default reflectivity values will be used.
        cmap : str or Colormap
            Color map name or Colormap instance.
        ax : Axes
            Axis to add the plot to. If None the current axis is used.

        """
        
        # Parse the parameters
        ax = self._parse_ax(ax)
        lon, lat = self._parse_lon_lat(lon, lat)
        vmin, vmax = self._parse_vmin_vmax(field, vmin, vmax)
        i, j = self._find_nearest_grid_indices(lon, lat)

        qm = ax.pcolormesh(self.y/1000, self.z/1000,
                           self.fields[field]['data'][:,:,i],
                           vmin=vmin, vmax=vmax, cmap=cmap)
            
        self.plots.append(qm)
        self.plot_vars.append(field)
        
        # Plot color bar
        if colorbar_flag:
            self.plot_colorbar(field, mappable=qm, orientation=orientation,
                               label=colorbar_label, fig=fig, ax=ax)
        
        # Set title
        if title_flag:
            self._set_lon_slice_title(field, lat, title, ax)
        
        return

    def plot_colorbar(self, field, mappable=None, orientation='vertical',
                      label=None, fig=None, ax=None, cax=None):
        """
        Plot a color bar to an axis.

        Parameters
        ----------
        mappable : ScalarMappable
            None will use the last mappable plotted.
        orientation : 'vertical' or 'horizontal'
            Color bar orientation on plot.
        label : str
            Color bar label.  None will use 'long_name' ['units'] for the
            last field plotted or an empty string if the field does not have
            these keys.
        fig : Figure
        ax : Axes
        """
        
        # Parse the parameters
        mappable = self._parse_mappable(mappable)
        
        # Create color bar instance
        cb = plt.colorbar(mappable=mappable, orientation=orientation,
                          ax=ax, cax=cax)

        # Set label
        self._set_colorbar_label(field, label, cb)

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
    
    def generate_xy_title(self, field, level):
        """
        Generate a title for a horizontal cross section plot.
        
        Parameters
        ----------
        field : str
            The field to be plotted.
        level : int
            The level to be plotted.
        
        Returns
        -------
        title : str
            The title for the plot.
        """
        
        date = self.start_time.isoformat() + 'Z'
        height = '%.2f km ' %(self.z[level] / 1000)
        field_name = self.fields[field]['long_name']
        
        return self.source + ' ' + height + date + '\n' + field_name
    
    def generate_lat_slice_title(self, field, lat):
        """
        Generate a title for a latitude slice plot.
        """
        
        # Parse paremeters
        lon, lat = self._parse_lon_lat(None, lat)
        
        date = self.start_time.isoformat() + 'Z'
        field_name = self.fields[field]['long_name']
        lat = '%.2f deg ' %lat
        
        return self.source + ' ' + lat + date + '\n' + field_name
    
    def generate_lon_slice_title(self, field, lon):
        """
        Generate a title for a latitude slice plot.
        """
        
        # Parse paremeters
        lon, lat = self._parse_lon_lat(lon, None)
        
        date = self.start_time.isoformat() + 'Z'
        field_name = self.fields[field]['long_name']
        lon = '%.2f deg ' %lon
        
        return self.source + ' ' + lon + date + '\n' + field_name
    
    def generate_colorbar_label(self, field):
        """
        Generate the color bar label for a plot.
        """
        
        field_name = self.fields[field]['long_name']
        units = self.fields[field]['units']
        
        return field_name + ' ' + '[' + units + ']'
    
    def _label_axes_xy(self, labels, ax):
        """
        Label the axes for horizontal cross section plots.
        """
        
        ax = self._parse_ax(ax)
        
        xlabel, ylabel = labels
        
        if xlabel is None:
            ax.set_xlabel('Eastward Distance from Origin [km]')
        else:
            ax.set_xlabel(xlabel)
        if ylabel is None:
            ax.set_ylabel('Northward Distance from Origin [km]')
        else:
            ax.set_ylabel(ylabel)
    
    def _set_xy_title(self, field, level, title, ax):
        """
        Set figure title for horizontal cross section plots.
        """
        
        ax = self._parse_ax(ax)
        
        if title is None:
            ax.set_title(self.generate_xy_title(field, level))
        else:
            ax.set_title(title)
            
    def _set_lat_slice_title(self, field, lat, title, ax):
        """
        Set the figure title for latitude slice plots.
        """
        
        ax = self._parse_ax(ax)
        
        if title is None:
            ax.set_title(self.generate_lat_slice_title(field, lat))
        else:
            ax.set_title(title)
            
    def _set_lon_slice_title(self, field, lon, title, ax):
        """
        Set the figure title for longitude slice plots.
        """
        
        ax = self._parse_ax(ax)
        
        if title is None:
            ax.set_title(self.generate_lon_slice_title(field, lon))
        else:
            ax.set_title(title)
            
    def _set_colorbar_label(self, field, label, cb):
        """
        """
        
        if label is None:
            cb.set_label(self.generate_colorbar_label(field))
        else:
            cb.set_label(label)

    def _check_basemap(self):
        """
        If Basemap is None raise a ValueError.
        """
        
        if self.basemap is None:
            raise ValueError('No basemap has been plotted yet')

    def _find_nearest_grid_indices(self, lon, lat):
        """
        Find the nearest (x, y) grid indices for a given latitude and
        longitude.
        """
        
        lon, lat = self._parse_lon_lat(lon, lat)
        x_cut, y_cut = self.proj(lon, lat)

        if self.debug:
            print "x = %.1f m at longitude = %.2f deg" %(x_cut, lon)
            print "y = %.1f m at latitude = %.2f deg" %(y_cut, lat)

        x_index = np.abs(self.x - x_cut).argmin()
        y_index = np.abs(self.y - y_cut).argmin()

        if self.debug:
            print "x index is %i", x_index
            print "y index is %i", y_index

        return x_index, y_index

    def _parse_vmin_vmax(self, field, vmin, vmax):
        """
        Parse vmin and vmax parameters. If these parameters are not
        supplied and the field does not have the valid min and max
        metadata, then vmin and vmax default to reflectivity values.
        """
        
        if vmin is None:
            
            if 'valid_min' in self.fields[field]:
                vmin = self.fields[field]['valid_min']
            else:
                vmin = -8

        if vmax is None:
            
            if 'valid_max' in self.fields[field]:
                vmax = self.fields[field]['valid_max']
            else:
                vmax = 64
                
        return vmin, vmax

    def _parse_lon_lat(self, lon, lat):
        """
        Parse latitude and longitude parameters.
        """
        
        if lat is None:
            lat = self.lat
        if lon is None:
            lon = self.lon
            
        return lon, lat

    def _parse_ax(self, ax):
        """
        Parse the axis parameter.
        """
        
        if ax is None:
            ax = plt.gca()
            
        return ax
    
    def _parse_mappable(self, mappable):
        """
        Parse the ScalarMappable.
        """
        
        if mappable is None:
            if len(self.plots) == 0:
                raise ValueError('No mappables (plots) found')
            else:
                mappable = self.plots[-1]
                
        return mappable