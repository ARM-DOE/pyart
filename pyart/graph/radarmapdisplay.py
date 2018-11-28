"""
pyart.graph.radarmapdisplay
===========================

Class for creating plots on a geographic map using a Radar object using Cartopy
for drawing maps.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    RadarMapDisplay

"""

import warnings

import numpy as np
import matplotlib.pyplot as plt
try:
    import cartopy
    _CARTOPY_AVAILABLE = True
except ImportError:
    _CARTOPY_AVAILABLE = False
try:
    import shapely.geometry as sgeom
    from copy import copy
    _LAMBERT_GRIDLINES = True
except ImportError:
    _LAMBERT_GRIDLINES = False

from .radardisplay import RadarDisplay
from .common import parse_ax_fig, parse_vmin_vmax, parse_cmap
from ..exceptions import MissingOptionalDependency


class RadarMapDisplay(RadarDisplay):
    """
    A display object for creating plots on a geographic map from data in a
    Radar object.

    This class is still a work in progress. Some functionality may not work
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
    grid_projection : cartopy.crs
        AzimuthalEquidistant cartopy projection centered on radar.
        Used to transform points into map projection

    """

    def __init__(self, radar, shift=(0.0, 0.0), grid_projection=None):
        """ Initialize the object. """
        # check that cartopy is available
        if not _CARTOPY_AVAILABLE:
            raise MissingOptionalDependency(
                "Cartopy is required to use RadarMapDisplay but is " +
                "not installed")

        # initalize the base class
        RadarDisplay.__init__(self, radar, shift=shift)

        # additional attributes needed for plotting on a cartopy map.
        if grid_projection is None:
            lat_0 = self.loc[0]
            lon_0 = self.loc[1]
            grid_projection = cartopy.crs.AzimuthalEquidistant(
                central_longitude=lon_0, central_latitude=lat_0)

        elif not isinstance(grid_projection, cartopy.crs.Projection):
            raise TypeError("grid_projection keyword must " +
                            "be a cartopy.crs object")
        self.grid_projection = grid_projection
        self.ax = None
        self._x0 = None     # x axis radar location in map coords (meters)
        self._y0 = None     # y axis radar location in map coords (meters)
        return

    def _check_ax(self):
        """ Check that a GeoAxes object exists, raise ValueError if not """
        if self.ax is None:
            raise ValueError('no GeoAxes plotted')

    def plot_ppi_map(
            self, field, sweep=0, mask_tuple=None,
            vmin=None, vmax=None, cmap=None, norm=None, mask_outside=False,
            title=None, title_flag=True,
            colorbar_flag=True, colorbar_label=None, ax=None, fig=None,
            lat_lines=None, lon_lines=None, projection=None,
            min_lon=None, max_lon=None, min_lat=None, max_lat=None,
            width=None, height=None, lon_0=None, lat_0=None,
            resolution='110m', shapefile=None, shapefile_kwargs=None,
            edges=True, gatefilter=None,
            filter_transitions=True, embelish=True, raster=False,
            ticks=None, ticklabs=None, alpha=None):
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
        ax : Cartopy GeoAxes instance
            If None, create GeoAxes instance using other keyword info.
            If provided, ax must have a Cartopy crs projection and projection
            kwarg below is ignored.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.
        lat_lines, lon_lines : array or None
            Locations at which to draw latitude and longitude lines.
            None will use default values which are resonable for maps of
            North America.
        projection : cartopy.crs class
            Map projection supported by cartopy. Used for all subsequent calls
            to the GeoAxes object generated. Defaults to LambertConformal
            centered on radar.
        min_lat, max_lat, min_lon, max_lon : float
            Latitude and longitude ranges for the map projection region in
            degrees.
        width, height : float
            Width and height of map domain in meters.
            Only this set of parameters or the previous set of parameters
            (min_lat, max_lat, min_lon, max_lon) should be specified.
            If neither set is specified then the map domain will be determined
            from the extend of the radar gate locations.
        shapefile : str
            Filename for a shapefile to add to map.
        shapefile_kwargs : dict
            Key word arguments used to format shapefile. Projection defaults
            to lat lon (cartopy.crs.PlateCarree())
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
            True by default. Set to False to supress drawing of coastlines
            etc.. Use for speedup when specifying shapefiles.
            Note that lat lon labels only work with certain projections.
        raster : bool
            False by default.  Set to true to render the display as a raster
            rather than a vector in call to pcolormesh.  Saves time in plotting
            high resolution data over large areas.  Be sure to set the dpi
            of the plot for your application if you save it as a vector format
            (i.e., pdf, eps, svg).
        alpha : float or None
            Set the alpha tranparency of the radar plot. Useful for
            overplotting radar over other datasets.

        """
        # parse parameters
        ax, fig = parse_ax_fig(ax, fig)
        vmin, vmax = parse_vmin_vmax(self._radar, field, vmin, vmax)
        cmap = parse_cmap(cmap, field)
        if lat_lines is None:
            lat_lines = np.arange(30, 46, 1)
        if lon_lines is None:
            lon_lines = np.arange(-110, -75, 1)
        lat_0 = self.loc[0]
        lon_0 = self.loc[1]

        # get data for the plot
        data = self._get_data(
            field, sweep, mask_tuple, filter_transitions, gatefilter)
        x, y = self._get_x_y(sweep, edges, filter_transitions)

        # mask the data where outside the limits
        if mask_outside:
            data = np.ma.masked_outside(data, vmin, vmax)

        # initialize instance of GeoAxes if not provided
        if hasattr(ax, 'projection'):
            projection = ax.projection
        else:
            if projection is None:
                # set map projection to LambertConformal if none is specified
                projection = cartopy.crs.LambertConformal(
                    central_longitude=lon_0, central_latitude=lat_0)
                warnings.warn("No projection was defined for the axes."
                              + " Overridding defined axes and using default "
                              + "axes.",
                              UserWarning)
            ax = plt.axes(projection=projection)

        if min_lon:
            ax.set_extent([min_lon, max_lon, min_lat, max_lat],
                          crs=cartopy.crs.PlateCarree())
        elif width:
            ax.set_extent([-width/2., width/2., -height/2., height/2.],
                          crs=self.grid_projection)

        # plot the data
        if norm is not None:  # if norm is set do not override with vmin/vmax
            vmin = vmax = None
        pm = ax.pcolormesh(x * 1000., y * 1000., data, alpha=alpha,
                           vmin=vmin, vmax=vmax, cmap=cmap,
                           norm=norm, transform=self.grid_projection)

        # plot as raster in vector graphics files
        if raster:
            pm.set_rasterized(True)

        # add embelishments
        if embelish is True:
            # Create a feature for States/Admin 1 regions at 1:resolution
            # from Natural Earth
            states_provinces = cartopy.feature.NaturalEarthFeature(
                category='cultural',
                name='admin_1_states_provinces_lines',
                scale=resolution,
                facecolor='none')
            ax.coastlines(resolution=resolution)
            ax.add_feature(states_provinces, edgecolor='gray')

            # labeling gridlines poses some difficulties depending on the
            # projection, so we need some projection-spectific methods
            if ax.projection in [cartopy.crs.PlateCarree(),
                                 cartopy.crs.Mercator()]:
                gl = ax.gridlines(xlocs=lon_lines, ylocs=lat_lines,
                                  draw_labels=True)
                gl.xlabels_top = False
                gl.ylabels_right = False

            elif isinstance(ax.projection, cartopy.crs.LambertConformal):
                fig.canvas.draw()
                ax.gridlines(xlocs=lon_lines, ylocs=lat_lines)

                # Label the end-points of the gridlines using the custom
                # tick makers:
                ax.xaxis.set_major_formatter(
                    cartopy.mpl.gridliner.LONGITUDE_FORMATTER)
                ax.yaxis.set_major_formatter(
                    cartopy.mpl.gridliner.LATITUDE_FORMATTER)
                if _LAMBERT_GRIDLINES:
                    lambert_xticks(ax, lon_lines)
                    lambert_yticks(ax, lat_lines)
            else:
                ax.gridlines(xlocs=lon_lines, ylocs=lat_lines)

        # plot the data and optionally the shape file
        # we need to convert the radar gate locations (x and y) which are in
        # km to meters we also need to give the original projection of the
        # data which is stored in self.grid_projection

        if shapefile is not None:
            from cartopy.io.shapereader import Reader
            if shapefile_kwargs is None:
                shapefile_kwargs = {}
            if 'crs' not in shapefile_kwargs:
                shapefile_kwargs['crs'] = cartopy.crs.PlateCarree()
            ax.add_geometries(Reader(shapefile).geometries(),
                              **shapefile_kwargs)

        if title_flag:
            self._set_title(field, sweep, title, ax)

        # add plot and field to lists
        self.plots.append(pm)
        self.plot_vars.append(field)

        if colorbar_flag:
            self.plot_colorbar(
                mappable=pm, label=colorbar_label, field=field, fig=fig,
                ax=ax, ticks=ticks, ticklabs=ticklabs)
        # keep track of this GeoAxes object for later
        self.ax = ax
        return

    def plot_point(self, lon, lat, symbol='ro', label_text=None,
                   label_offset=(None, None), **kwargs):
        """
        Plot a point on the current map.

        Additional arguments are passed to ax.plot.

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
            text relative to the point. A value of None will use 0.01.

        """
        self._check_ax()
        lon_offset, lat_offset = label_offset
        if lon_offset is None:
            lon_offset = 0.01
        if lat_offset is None:
            lat_offset = 0.01
        if 'transform' not in kwargs.keys():
            kwargs['transform'] = cartopy.crs.PlateCarree()
        self.ax.plot(lon, lat, symbol, **kwargs)
        if label_text is not None:
            # the "plt.annotate call" does not have a "transform=" keyword,
            # so for this one we transform the coordinates with a Cartopy call.
            x_text, y_text = self.ax.projection.transform_point(
                lon + lon_offset, lat + lat_offset,
                src_crs=kwargs['transform'])
            self.ax.annotate(label_text, xy=(x_text, y_text))

    def plot_line_geo(self, line_lons, line_lats, line_style='r-', **kwargs):
        """
        Plot a line segments on the current map given values in lat and lon.

        Additional arguments are passed to ax.plot.

        Parameters
        ----------
        line_lons : array
            Longitude of line segment to plot.
        line_lats : array
            Latitude of line segment to plot.
        line_style : str
            Matplotlib compatible string which specifies the line style.

        """
        self._check_ax()
        if 'transform' not in kwargs.keys():
            kwargs['transform'] = cartopy.crs.PlateCarree()
        self.ax.plot(line_lons, line_lats, line_style, **kwargs)

    def plot_line_xy(self, line_x, line_y, line_style='r-', **kwargs):
        """
        Plot a line segments on the current map given radar x, y values.

        Additional arguments are passed to ax.plot.

        Parameters
        ----------
        line_x : array
            X location of points to plot in meters from the radar.
        line_y : array
            Y location of points to plot in meters from the radar.
        line_style : str, optional
            Matplotlib compatible string which specifies the line style.

        """
        self._check_ax()
        if 'transform' not in kwargs.keys():
            kwargs['transform'] = self.grid_projection
        self.ax.plot(line_x, line_y, line_style, **kwargs)

    def plot_range_ring(self, range_ring_location_km, npts=360,
                        line_style='k-', **kwargs):
        """
        Plot a single range ring on the map.

        Additional arguments are passed to ax.plot.

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
        self._check_ax()
        angle = np.linspace(0., 2.0 * np.pi, npts)
        xpts = range_ring_location_km * 1000. * np.sin(angle)
        ypts = range_ring_location_km * 1000. * np.cos(angle)
        self.plot_line_xy(xpts, ypts, line_style=line_style, **kwargs)


# These methods are a hack to allow gridlines when the projection is lambert
# http://nbviewer.jupyter.org/gist/ajdawson/dd536f786741e987ae4e

def find_side(ls, side):
    """
    Given a shapely LineString which is assumed to be rectangular, return the
    line corresponding to a given side of the rectangle.
    """
    minx, miny, maxx, maxy = ls.bounds
    points = {'left': [(minx, miny), (minx, maxy)],
              'right': [(maxx, miny), (maxx, maxy)],
              'bottom': [(minx, miny), (maxx, miny)],
              'top': [(minx, maxy), (maxx, maxy)]}
    return sgeom.LineString(points[side])


def lambert_xticks(ax, ticks):
    """Draw ticks on the bottom x-axis of a Lambert Conformal projection."""
    def te(xy):
        return xy[0]

    def lc(t, n, b):
        return np.vstack((np.zeros(n) + t, np.linspace(b[2], b[3], n))).T

    xticks, xticklabels = _lambert_ticks(ax, ticks, 'bottom', lc, te)
    ax.xaxis.tick_bottom()
    ax.set_xticks(xticks)
    ax.set_xticklabels([ax.xaxis.get_major_formatter()(xtick) for
                        xtick in xticklabels])


def lambert_yticks(ax, ticks):
    """Draw ticks on the left y-axis of a Lambert Conformal projection."""
    def te(xy):
        return xy[1]

    def lc(t, n, b):
        return np.vstack((np.linspace(b[0], b[1], n), np.zeros(n) + t)).T

    yticks, yticklabels = _lambert_ticks(ax, ticks, 'left', lc, te)
    ax.yaxis.tick_left()
    ax.set_yticks(yticks)
    ax.set_yticklabels([ax.yaxis.get_major_formatter()(ytick) for
                        ytick in yticklabels])


def _lambert_ticks(ax, ticks, tick_location, line_constructor, tick_extractor):
    """Get the tick locations and labels for a Lambert Conformal projection."""
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
