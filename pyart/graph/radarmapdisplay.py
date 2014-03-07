"""
pyart.graph.radarmapdisplay
=========================

Class for creating plots on a geographic map using a Radar objects.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    RadarMapDisplay

"""

import numpy as np
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import pyproj

from .radardisplay import RadarDisplay


class RadarMapDisplay(RadarDisplay):
    """
    A display object for creating plots on a geographic map from data in a
    Radar object.

    This class is still an in progress addition to Py-ART.  Some
    functionality may not work correctly. Please report any problems to
    the Py-ART GitHub Issue Tracker.

    Parameters
    ----------
    radar : Radar
        Radar object to use for creating plots.
    shift : (float, float)
        Shifts in km to offset the calculated x and y locations.

    Attributes
    ----------
    plots : list
        List of plots created.
    plot_vars : list
        List of fields plotted, order matches plot list.
    cbs : list
        List of colorbars created.
    radar_name : str
        Name of radar.
    origin : str
        'Origin' or 'Radar'.
    shift : (float, float)
        Shift in meters.
    x, y, z : array
        Cartesian location of a sweep in meters.
    loc : (float, float)
        Latitude and Longitude of radar in degrees.
    time_begin : datetime
        Beginning time of first radar scan.
    starts : array
        Starting ray index for each sweep.
    ends : array
        Ending ray index for each sweep.
    fields : dict
        Radar fields.
    scan_type : str
        Scan type.
    ranges : array
        Gate ranges in meters.
    azimuths : array
        Azimuth angle in degrees.
    elevations : array
        Elevations in degrees.
    fixed_angle : array
        Scan angle in degrees.
    proj : Proj
        Object for performing cartographic transformations specific to the
        geographic map plotted.
    basemap : Basemap
        Last plotted basemap, None when no basemap has been plotted.


    """

    def __init__(self, radar, shift=(0.0, 0.0)):
        """ Initialize the object. """

        # initalize the base class
        RadarDisplay.__init__(self, radar, shift=shift)

        # additional Map attributes.
        self.basemap = None
        self.proj = pyproj.Proj(proj='lcc', datum='NAD83',
                                lat_0=self.loc[0], lon_0=self.loc[1],
                                x_0=0.0, y_0=0.0)
        return

    def plot_ppi_map(self, field, sweep=0, mask_tuple=None,
                     vmin=None, vmax=None, cmap='jet', mask_outside=True,
                     title=None, title_flag=True,
                     colorbar_flag=True, colorbar_label=None,
                     ax=None, fig=None,
                     lat_lines=None, lon_lines=None,
                     auto_range=True,
                     min_lon=-92, max_lon=-86, min_lat=40, max_lat=44,
                     resolution='h', shapefile=None):
        """
        Plot a PPI volume sweep onto a geographic map.

        Parameters
        ----------
        field : str
            Field to plot.
        sweep : int, optional
            Sweep number to plot.

        Other Parameters
        ----------------
        mask_tuple : (str, float)
            Tuple containing the field name and value below which to mask
            field prior to plotting, for example to mask all data where
            NCP < 0.5 set mask_tuple to ['NCP', 0.5]. None performs no masking.
        vmin : float
            Luminance minimum value, None for default value.
        vmax : float
            Luminance maximum value, None for default value.
        cmap : str
            Matplotlib colormap name.
        mask_outside : bool
            True to mask data outside of vmin, vmax.  False performs no
            masking.
        title : str
            Title to label plot with, None to use default title generated from
            the field and tilt parameters. Parameter is ignored if title_flag
            is False.
        title_flag : bool
            True to add a title to the plot, False does not add a title.
        colorbar_flag : bool
            True to add a colorbar with label to the axis.  False leaves off
            the colorbar.
        colorbar_label : str
            Colorbar label, None will use a default label generated from the
            field information.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.
        lat_lines, lon_lines : array or None
            Locations at which to draw latitude and longitude lines.
            None will use default values which are resonable for maps of
            North America.
        auto_range : bool
            True to determine map ranges from the extend of the radar data.
            False will use the min_lon, max_lon, min_lat, and max_lat
            parameters for the map range.
        min_lat, max_lat, min_lon, max_lon : float
            Latitude and longitude ranges for the map projection region in
            degrees.  These parameter are not used if auto_range is True.
        shapefile : str
            Filename for a ESRI shapefile as background (untested).
        resolution : 'c', 'l', 'i', 'h', or 'f'.
            Resolution of boundary database to use. See Basemap documentation
            for details.

        """
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
        vmin, vmax = self._parse_vmin_vmax(field, vmin, vmax)
        if lat_lines is None:
            lat_lines = np.arange(30, 46, 1)
        if lon_lines is None:
            lon_lines = np.arange(-110, -75, 1)

        # get data for the plot
        data = self._get_data(field, sweep, mask_tuple)
        x, y = self._get_x_y(field, sweep)

        # mask the data where outside the limits
        if mask_outside:
            data = np.ma.masked_outside(data, vmin, vmax)

        # determine map region
        lon, lat = self.proj(x * 1000.0, y * 1000.0, inverse=True)

        if auto_range:
            max_lat = lat.max()
            max_lon = lon.max()
            min_lat = lat.min()
            min_lon = lon.min()

        # plot the basemap
        basemap = Basemap(llcrnrlon=min_lon, llcrnrlat=min_lat,
                          urcrnrlon=max_lon, urcrnrlat=max_lat,
                          projection='mill', area_thresh=1000,
                          resolution=resolution, ax=ax)
        basemap.drawcoastlines(linewidth=1.25)
        basemap.drawstates()
        basemap.drawparallels(lat_lines, labels=[True, False, False, False])
        basemap.drawmeridians(lon_lines, labels=[False, False, False, True])
        self.basemap = basemap

        # plot the data and optionally the shape file
        xd, yd = basemap(lon, lat)
        pm = basemap.pcolormesh(xd, yd, data, vmin=vmin, vmax=vmax, cmap=cmap)

        if shapefile is not None:
            basemap.readshapefile(shapefile, 'shapefile', ax=ax)

        if title_flag:
            self._set_title(field, sweep, title, ax)

        # add plot and field to lists
        self.plots.append(pm)
        self.plot_vars.append(field)

        if colorbar_flag:
            self.plot_colorbar(mappable=pm, label=colorbar_label,
                               field=field, fig=fig)
        return
