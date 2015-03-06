"""
pyart.graph.radarmapdisplay
===========================

Class for creating plots on a geographic map using a Radar object.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    RadarMapDisplay

"""

import numpy as np
from mpl_toolkits.basemap import Basemap

from .radardisplay import RadarDisplay


class RadarMapDisplay(RadarDisplay):
    """
    A display object for creating plots on a geographic map from data in a
    Radar object.

    This class is still a work in progress.  Some functionality may not work
    correctly. Please report any problems to the Py-ART GitHub Issue Tracker.

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

        # additional attributes needed for plotting on a basemap.
        self.basemap = None
        self._x0 = None     # x axis radar location in map coords (meters)
        self._y0 = None     # y axis radar location in map coords (meters)
        return

    def _check_basemap(self):
        """ Check that basemap is not None, raise ValueError if it is. """
        if self.basemap is None:
            raise ValueError('no basemap plotted')

    def plot_ppi_map(self, field, sweep=0, mask_tuple=None,
                     vmin=None, vmax=None, cmap='jet', mask_outside=True,
                     title=None, title_flag=True,
                     colorbar_flag=True, colorbar_label=None,
                     ax=None, fig=None,
                     lat_lines=None, lon_lines=None,
                     projection='lcc', area_thresh=10000,
                     min_lon=None, max_lon=None, min_lat=None, max_lat=None,
                     width=None, height=None, lon_0=None, lat_0=None,
                     resolution='h', shapefile=None, edges=True,
                     filter_transitions=True, embelish=True, **kwargs):
        """
        Plot a PPI volume sweep onto a geographic map.

        Additional arguments are passed to Basemap.

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
        projection : str
            Map projection supported by basemap.  The use of cylindrical
            projections (mill, merc, cyl, etc) is not recommended as they
            exhibit large distortions at high latitudes.  Equal area
            (aea, laea), conformal (lcc, tmerc, stere) or equidistant
            projection (aeqd, cass) work well even at high latitudes.
        area_thresh : float
            Coastline or lake with an area smaller than area_thresh in
            km^2 will not be plotted.
        min_lat, max_lat, min_lon, max_lon : float
            Latitude and longitude ranges for the map projection region in
            degrees.
        width, height : float
            Width and height of map domain in meters.
            Only this set of parameters or the previous set of parameters
            (min_lat, max_lat, min_lon, max_lon) should be specified.
            If neither set is specified then the map domain will be determined
            from the extend of the radar gate locations.
        lon_0, lat_0 : float
            Center of the map domain in degrees.  If the default, None is used
            the latitude and longitude of the radar will be used.
        shapefile : str
            Filename for a ESRI shapefile as background (untested).
        resolution : 'c', 'l', 'i', 'h', or 'f'.
            Resolution of boundary database to use. See Basemap documentation
            for details.
        filter_transitions : bool
            True to remove rays where the antenna was in transition between
            sweeps from the plot.  False will include these rays in the plot.
            No rays are filtered when the antenna_transition attribute of the
            underlying radar is not present.
        edges : bool
            True will interpolate and extrapolate the gate edges from the
            range, azimuth and elevations in the radar, treating these
            as specifying the center of each gate.  False treats these
            coordinates themselved as the gate edges, resulting in a plot
            in which the last gate in each ray and the entire last ray are not
            not plotted.
        embelish: bool
            True by default. Set to false to supress drawing of coastlines
            etc.. Use for speedup when specifying shapefiles.

        """
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
        vmin, vmax = self._parse_vmin_vmax(field, vmin, vmax)
        if lat_lines is None:
            lat_lines = np.arange(30, 46, 1)
        if lon_lines is None:
            lon_lines = np.arange(-110, -75, 1)
        if lat_0 is None:
            lat_0 = self.loc[0]
        if lon_0 is None:
            lon_0 = self.loc[1]

        # get data for the plot
        data = self._get_data(field, sweep, mask_tuple, filter_transitions)
        x, y = self._get_x_y(field, sweep, edges, filter_transitions)

        # mask the data where outside the limits
        if mask_outside:
            data = np.ma.masked_outside(data, vmin, vmax)

        # plot the basemap
        using_corners = (None not in [min_lon, min_lat, max_lon, max_lat])
        if using_corners:
            basemap = Basemap(llcrnrlon=min_lon, llcrnrlat=min_lat,
                              urcrnrlon=max_lon, urcrnrlat=max_lat,
                              lat_0=lat_0, lon_0=lon_0,
                              projection=projection, area_thresh=area_thresh,
                              resolution=resolution, ax=ax, **kwargs)
        else:   # using width and height
            # map domain determined from location of radar gates
            if width is None:
                width = (x.max() - y.min()) * 1000.
            if height is None:
                height = (y.max() - y.min()) * 1000.
            basemap = Basemap(width=width, height=height,
                              lon_0=lon_0, lat_0=lat_0,
                              projection=projection, area_thresh=area_thresh,
                              resolution=resolution, ax=ax, **kwargs)

        # add embelishments
        if embelish == True:
            basemap.drawcoastlines(linewidth=1.25)
            basemap.drawstates()
            basemap.drawparallels(lat_lines, labels=[True, False, False, False])
            basemap.drawmeridians(lon_lines, labels=[False, False, False, True])
        self.basemap = basemap
        self._x0, self._y0 = basemap(self.loc[1], self.loc[0])

        # plot the data and optionally the shape file
        # we need to convert the radar gate locations (x and y) which are in
        # km to meters as well as add the map coordiate radar location
        # which is given by self._x0, self._y0.
        pm = basemap.pcolormesh(self._x0 + x * 1000., self._y0 + y * 1000.,
                                data, vmin=vmin, vmax=vmax, cmap=cmap)

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

    def plot_point(self, lon, lat, symbol='ro', label_text=None,
                   label_offset=(None, None), **kwargs):
        """
        Plot a point on the current map.

        Additional arguments are passed to basemap.plot.

        Parameters
        ----------
        lon : float
            Longitude of point to plot.
        lat : float
            Latitude of point to plot.
        symbol : str
            Matplotlib compatible string which specified the symbol of the
            point.
        label_text : str, optional.
            Text to label symbol with.  If None no label will be added.
        label_offset : [float, float]
            Offset in lon, lat degrees for the bottom left corner of the label
            text relative to the point. A value of None will use 0.01 de

        """
        self._check_basemap()
        lon_offset, lat_offset = label_offset
        if lon_offset is None:
            lon_offset = 0.01
        if lat_offset is None:
            lat_offset = 0.01
        self.basemap.plot(lon, lat, symbol, latlon=True, **kwargs)
        if label_text is not None:
            # basemap does not have a text method so we must determine
            # the x and y points and plot them on the basemap's axis.
            x_text, y_text = self.basemap(lon + lon_offset,
                                          lat + lat_offset)
            self.basemap.ax.text(x_text, y_text, label_text)

    def plot_line_geo(self, line_lons, line_lats, line_style='r-', **kwargs):
        """
        Plot a line segments on the current map given values in lat and lon.

        Additional arguments are passed to basemap.plot.

        Parameters
        ----------
        line_lons : array
            Longitude of line segment to plot.
        line_lats : array
            Latitude of line segment to plot.
        line_style : str
            Matplotlib compatible string which specifies the line style.

        """
        self._check_basemap()
        self.basemap.plot(line_lons, line_lats, line_style, latlon=True,
                          **kwargs)

    def plot_line_xy(self, line_x, line_y, line_style='r-', **kwargs):
        """
        Plot a line segments on the current map given radar x, y values.

        Additional arguments are passed to basemap.plot.

        Parameters
        ----------
        line_x : array
            X location of points to plot in meters from the radar.
        line_y : array
            Y location of points to plot in meters from the radar.
        line_style : str, optional
            Matplotlib compatible string which specifies the line style.

        """
        self._check_basemap()
        self.basemap.plot(self._x0 + line_x, self._y0 + line_y, line_style,
                          latlon=False, **kwargs)

    def plot_range_ring(self, range_ring_location_km, npts=360,
                        line_style='k-', **kwargs):
        """
        Plot a single range ring on the map.

        Additional arguments are passed to basemap.plot.

        Parameters
        ----------
        range_ring_location_km : float
            Location of range ring in km.
        npts: int
            Number of points in the ring, higher for better resolution.
        line_style : str
            Matplotlib compatible string which specified the line
            style of the ring.

        """
        self._check_basemap()
        angle = np.linspace(0., 2.0 * np.pi, npts)
        xpts = range_ring_location_km * 1000. * np.sin(angle)
        ypts = range_ring_location_km * 1000. * np.cos(angle)
        self.plot_line_xy(xpts, ypts, line_style=line_style, **kwargs)
