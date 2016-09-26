"""
pyart.graph.radarmapdisplay_cartopy
===================================

Class for creating plots on a geographic map using a Radar object using Cartopy
for drawing maps.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    RadarMapDisplayCartopy

"""

import numpy as np
try:
    import cartopy
    _CARTOPY_AVAILABLE = True
except ImportError:
    _CARTOPY_AVAILABLE = False

from .radardisplay import RadarDisplay
from .common import parse_ax_fig, parse_norm_vmin_vmax, parse_cmap
from ..exceptions import MissingOptionalDependency


class RadarMapDisplayCartopy(RadarDisplay):
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
    origin : str
        'Origin' or 'Radar'.
    shift : (float, float)
        Shift in meters.
    loc : (float, float)
        Latitude and Longitude of radar in degrees.
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
        # check that basemap is available
        if not _CARTOPY_AVAILABLE:
            raise MissingOptionalDependency(
                "Cartopy is required to use RadarMapDisplayCartopy but is " +
                "not installed")

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

    def plot_ppi_map(
            self, field, sweep=0, mask_tuple=None,
            vmin=None, vmax=None, cmap=None, norm=None, mask_outside=False,
            title=None, title_flag=True,
            colorbar_flag=True, colorbar_label=None, ax=None, fig=None,
            lat_lines=None, lon_lines=None,
            projection=None, area_thresh=10000,
            min_lon=None, max_lon=None, min_lat=None, max_lat=None,
            width=None, height=None, lon_0=None, lat_0=None,
            resolution='110m', shapefile=None, edges=True, gatefilter=None,
            basemap=None, filter_transitions=True, embelish=True,
            ticks=None, ticklabs=None, **kwargs):
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
            Parameter is ignored is norm is not None.
        vmax : float
            Luminance maximum value, None for default value.
            Parameter is ignored is norm is not None.
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
            the field and tilt parameters. Parameter is ignored if title_flag
            is False.
        title_flag : bool
            True to add a title to the plot, False does not add a title.
        colorbar_flag : bool
            True to add a colorbar with label to the axis.  False leaves off
            the colorbar.
        ticks : array
            Colorbar custom tick label locations.
        ticklabs : array
            Colorbar custom tick labels.
        colorbar_label : str
            Colorbar label, None will use a default label generated from the
            field information.
        ax : Axis
            Axis to plot on. If provided, must have a Cartopy crs projection.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.
        lat_lines, lon_lines : array or None
            Locations at which to draw latitude and longitude lines.
            None will use default values which are resonable for maps of
            North America.
        projection : cartopy.crs class
            Map projection supported by cartopy.
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
        resolution : '10m', '50m', '110m'.
            Resolution of NaturalEarthFeatures to use. See Cartopy
            documentation for details.
        gatefilter : GateFilter
            GateFilter instance. None will result in no gatefilter mask being
            applied to data.
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
        basemap: Cartopy GeoAxes instance
            If None, create GeoAxes instance using other keyword info.
            If not None, use the user-specifed basemap instance.

        """
        # parse parameters
        ax, fig = parse_ax_fig(ax, fig)
        norm, vmin, vmax = parse_norm_vmin_vmax(
            norm, self._radar, field, vmin, vmax)
        cmap = parse_cmap(cmap, field)
        if lat_lines is None:
            lat_lines = np.arange(30, 46, 1)
        if lon_lines is None:
            lon_lines = np.arange(-110, -75, 1)
        if lat_0 is None:
            lat_0 = self.loc[0]
        if lon_0 is None:
            lon_0 = self.loc[1]

        # get data for the plot
        data = self._get_data(
            field, sweep, mask_tuple, filter_transitions, gatefilter)
        x, y = self._get_x_y(sweep, edges, filter_transitions)

        # mask the data where outside the limits
        if mask_outside:
            data = np.ma.masked_outside(data, vmin, vmax)

        proj = cartopy.crs.LambertConformal(
            central_longitude=lon_0, central_latitude=lat_0,
            standard_parallels=(lat_0, ))
        self.proj = proj

        # create the basemap if not provided
        import matplotlib.pyplot as plt

        if hasattr(ax, 'projection'):
            basemap = ax
        elif hasattr(ax, 'projection'):
            basemap = basemap
        else:
            if projection is None:
                basemap = plt.axes(projection=proj)
            else:
                basemap = plt.axes(projection=projection)
        self.basemap = basemap
        if min_lon:
            basemap.set_extent([min_lon, max_lon, min_lat, max_lat],
                               crs=cartopy.crs.PlateCarree())

        # add embelishments
        if embelish is True:
            # shpfile = cartopy.io.shapereader.gshhs(resolution)
            # shp = cartopy.io.shapereader.Reader(shpfile)
            # basemap.add_geometries(
            #    shp.geometries(), cartopy.crs.PlateCarree(), edgecolor='black',
            #    facecolor='none')
            #basemap.drawstates()
            if basemap.projection in [cartopy.crs.PlateCarree(),
                                      cartopy.crs.Mercator()]:
                gl = basemap.gridlines(xlocs=lon_lines, ylocs=lat_lines,
                                       draw_labels=True)
                gl.xlabels_top = False
                gl.ylabels_right = False

            elif basemap.projection == proj:
                fig.canvas.draw()
                basemap.gridlines(xlocs=lon_lines, ylocs=lat_lines)

                # Label the end-points of the gridlines using the custom
                # tick makers:
                basemap.xaxis.set_major_formatter(
                    cartopy.mpl.gridliner.LONGITUDE_FORMATTER)
                basemap.yaxis.set_major_formatter(
                    cartopy.mpl.gridliner.LATITUDE_FORMATTER)
                lambert_xticks(basemap, lon_lines)
                lambert_yticks(basemap, lat_lines)

            else:
                basemap.gridlines(xlocs=lon_lines, ylocs=lat_lines)

            # Create a feature for States/Admin 1 regions at 1:resolution
            # from Natural Earth
            states_provinces = cartopy.feature.NaturalEarthFeature(
                category='cultural',
                name='admin_1_states_provinces_lines',
                scale=resolution,
                facecolor='none')
            basemap.coastlines(resolution=resolution)
            basemap.add_feature(states_provinces, edgecolor='gray')

        # plot the data and optionally the shape file
        # we need to convert the radar gate locations (x and y) which are in
        # km to meters as well as add the map coordiate radar location
        # which is given by self._x0, self._y0.
        pm = basemap.pcolormesh(
            x * 1000., y * 1000., data,
            vmin=vmin, vmax=vmax, cmap=cmap, norm=norm, transform=proj)

        if shapefile is not None:
            from cartopy.io.shapereader import Reader
            basemap.add_geometries(Reader(shapefile).geometries(),
                                   cartopy.crs.PlateCarree())

        if title_flag:
            self._set_title(field, sweep, title, basemap)

        # add plot and field to lists
        self.plots.append(pm)
        self.plot_vars.append(field)

        if colorbar_flag:
            self.plot_colorbar(
                mappable=pm, label=colorbar_label, field=field, fig=fig,
                ax=basemap, ticks=ticks, ticklabs=ticklabs)
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
        self.basemap.plot(lon, lat, symbol,
                          transform=cartopy.crs.PlateCarree(), **kwargs)
        if label_text is not None:
            # the "plt.annotate call" does not have a "transform=" keyword,
            # so for this one we transform the coordinates with a Cartopy call.
            x_text, y_text = self.basemap.projection.transform_point(
                                            lon + lon_offset,
                                            lat + lat_offset,
                                            src_crs=cartopy.crs.PlateCarree())
            self.basemap.annotate(label_text, xy=(x_text, y_text))

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
        self.basemap.plot(
            line_lons, line_lats, line_style,
            transform=cartopy.crs.PlateCarree(), **kwargs)

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
        self.basemap.plot(line_x, line_y, line_style, transform=self.proj,
                          **kwargs)

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
        # The RadarDisplay.plot_range_rings uses a col parameter to specify
        # the line color, deal with this here.
        if 'col' in kwargs:
            color = kwargs.pop('col')
            kwargs['c'] = color
        if 'ax' in kwargs:
            kwargs.pop('ax')
        self._check_basemap()
        angle = np.linspace(0., 2.0 * np.pi, npts)
        xpts = range_ring_location_km * 1000. * np.sin(angle)
        ypts = range_ring_location_km * 1000. * np.cos(angle)
        self.plot_line_xy(xpts, ypts, line_style=line_style, **kwargs)

# Thses methods are a hack to allow gridlines when the projection is lambert
# http://nbviewer.jupyter.org/gist/ajdawson/dd536f786741e987ae4e
def find_side(ls, side):
    """
    Given a shapely LineString which is assumed to be rectangular, return the
    line corresponding to a given side of the rectangle.

    """
    import shapely.geometry as sgeom
    minx, miny, maxx, maxy = ls.bounds
    points = {'left': [(minx, miny), (minx, maxy)],
              'right': [(maxx, miny), (maxx, maxy)],
              'bottom': [(minx, miny), (maxx, miny)],
              'top': [(minx, maxy), (maxx, maxy)],}
    return sgeom.LineString(points[side])


def lambert_xticks(ax, ticks):
    """Draw ticks on the bottom x-axis of a Lambert Conformal projection."""
    te = lambda xy: xy[0]
    lc = lambda t, n, b: np.vstack((np.zeros(n) + t,
                                    np.linspace(b[2], b[3], n))).T
    xticks, xticklabels = _lambert_ticks(ax, ticks, 'bottom', lc, te)
    ax.xaxis.tick_bottom()
    ax.set_xticks(xticks)
    ax.set_xticklabels([ax.xaxis.get_major_formatter()(xtick) for
                        xtick in xticklabels])


def lambert_yticks(ax, ticks):
    """Draw ricks on the left y-axis of a Lamber Conformal projection."""
    te = lambda xy: xy[1]
    lc = lambda t, n, b: np.vstack((np.linspace(b[0], b[1], n),
                                    np.zeros(n) + t)).T
    yticks, yticklabels = _lambert_ticks(ax, ticks, 'left', lc, te)
    ax.yaxis.tick_left()
    ax.set_yticks(yticks)
    ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for
                        ytick in yticklabels])

def _lambert_ticks(ax, ticks, tick_location, line_constructor, tick_extractor):
    """Get the tick locations and labels for a Lambert Conformal projection."""
    import shapely.geometry as sgeom
    from copy import copy
    outline_patch = sgeom.LineString(
                                ax.outline_patch.get_path().vertices.tolist())
    axis = find_side(outline_patch, tick_location)
    n_steps = 30
    extent = ax.get_extent(cartopy.crs.PlateCarree())
    _ticks = []
    for t in ticks:
        xy = line_constructor(t, n_steps, extent)
        proj_xyz = ax.projection.transform_points(cartopy.crs.Geodetic(),
                                                  xy[:, 0], xy[:, 1])
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
