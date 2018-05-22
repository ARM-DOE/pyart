"""
pyart.graph.crfmapdisplay
===========================

Function for creating plots on a geographic map using a Radar object or Radar objects.

"""

import numpy as np
try:
    from mpl_toolkits.basemap import Basemap
    _BASEMAP_AVAILABLE = True
except ImportError:
    _BASEMAP_AVAILABLE = False

from ..map import grid_from_radars
from .radarmapdisplay import RadarMapDisplay
from .common import parse_ax_fig, parse_vmin_vmax, parse_cmap
            
def plot_crf(
        radars, field, grid_shape, grid_limits,
        vmin=None, vmax=None,
        mask_outside=False, title=None, title_flag=True,
        cmap=None, norm=None,
        colorbar_flag=True, colorbar_label=None, ax=None, fig=None,
        lat_lines=None, lon_lines=None,
        projection='lcc', area_thresh=10000,
        min_lon=None, max_lon=None, min_lat=None, max_lat=None,
        width=None, height=None, lon_0=None, lat_0=None,
        resolution='h', shapefile=None, gatefilter=None,
        basemap=None, embelish=True,
        ticks=None, ticklabs=None, alpha=None, **kwargs):
    """
    Plot CRF onto a geographic map.

    Additional arguments are passed to Basemap.

    Parameters
    ----------
    radars : Radar or tuple of Radar objects.
        Radar objects which will be mapped to the Cartesian grid.
    field : str
        Field to plot.
    grid_shape : 3-tuple of floats
        Number of points in the grid (z, y, x).
    grid_limits : 3-tuple of 2-tuples
        Minimum and maximum grid location (inclusive) in meters for the
        z, y, x coordinates.
    gridding_algo : 'map_to_grid' or 'map_gates_to_grid'
        Algorithm to use for gridding.  'map_to_grid' finds all gates within
        a radius of influence for each grid point, 'map_gates_to_grid' maps
        each radar gate onto the grid using a radius of influence and is
        typically significantly faster.
    Other Parameters
    ----------------
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
        Axis to plot on. None will use the current axis.
    fig : Figure
        Figure to add the colorbar to. None will use the current figure.
    lat_lines, lon_lines : array or None
        Locations at which to draw latitude and longitude lines.
        None will use default values which are resonable for maps of
        North America.
    projection : str
        Map projection supported by basemap.  The use of cylindrical
        projections (mill, merc, etc) is not recommended as they
        exhibit large distortions at high latitudes.  Equal area
        (aea, laea), conformal (lcc, tmerc, stere) or equidistant
        projection (aeqd, cass) work well even at high latitudes.
        The cylindrical equidistant projection (cyl) is not supported as
        coordinate transformations cannot be performed.
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
    gatefilter : GateFilter
        GateFilter instance. None will result in no gatefilter mask being
        applied to data.
    embelish: bool
        True by default. Set to false to supress drawing of coastlines
        etc.. Use for speedup when specifying shapefiles.
    basemap: Basemap instance
        If None, create basemap instance using other keyword info.
        If not None, use the user-specifed basemap instance.
    alpha : float or None
        Set the alpha tranparency of the radar plot. Useful for
        overplotting radar over other datasets.

    """
    # parse parameters
    first_radar = radars[0]
    ax, fig = parse_ax_fig(ax, fig)
    vmin, vmax = parse_vmin_vmax(first_radar, field, vmin, vmax)
    if lat_lines is None:
        lat_lines = np.arange(30, 46, 1)
    if lon_lines is None:
        lon_lines = np.arange(-110, -75, 1)
    if lat_0 is None:
        lat_0 = RadarMapDisplay(first_radar).loc[0]
    if lon_0 is None:
        lon_0 = RadarMapDisplay(first_radar).loc[1]

    if gatefilter == None:
        gatefilter = (None,)*len(radars)

    # get data for the plot
    data = grid_from_radars(
                        radars,
                        gatefilters=gatefilter,
                        grid_shape=grid_shape,
                        grid_limits=grid_limits,
                        fields=[field])

    # mask the data where outside the limits
    if mask_outside:
        data = np.ma.masked_outside(data, vmin, vmax)

    # create the basemap if not provided
    if type(basemap) != Basemap:
        using_corners = (None not in [min_lon, min_lat, max_lon, max_lat])
        if using_corners:
            basemap = Basemap(
                llcrnrlon=min_lon, llcrnrlat=min_lat,
                urcrnrlon=max_lon, urcrnrlat=max_lat,
                lat_0=lat_0, lon_0=lon_0, projection=projection,
                area_thresh=area_thresh, resolution=resolution, ax=ax,
                **kwargs)
        else:   # using width and height
            # map domain determined from location of radar gates
            if width is None:
                width = data.x['data'].max() - data.x['data'].min()
            if height is None:
                height = data.y['data'].max() - data.y['data'].min()
            basemap = Basemap(
                width=width, height=height, lon_0=lon_0, lat_0=lat_0,
                projection=projection, area_thresh=area_thresh,
                resolution=resolution, ax=ax, **kwargs)

    # The cylindrical equidistant projection does not support conversions
    # from geographic (lon/lat) to map projection (x/y) coordinates and
    # therefore cannot be used.
    if basemap.projection == 'cyl':
        raise ValueError(
            'The cylindrical equidistant projection is not supported')

    # add embelishments
    if embelish is True:
        # if you don't have shapefile, you can uncomment two lines below.
        # basemap.drawcoastlines(linewidth=1.25)
        # basemap.drawstates()
        basemap.drawparallels(
            lat_lines, labels=[True, False, False, False])
        basemap.drawmeridians(
            lon_lines, labels=[False, False, False, True])
    RadarMapDisplay(first_radar).basemap = basemap
    x0, y0 = basemap(RadarMapDisplay(first_radar).loc[1], RadarMapDisplay(first_radar).loc[0])

    # plot the data and optionally the shape file
    # we need to read the radar gate locations (x and y) which are in
    # meters as well as add the map coordiate radar location
    # which is given by self._x0, self._y0.
    if norm is not None:  # if level is set do not override with vmin/vmax
        vmin = vmax = None
    pm = basemap.pcolormesh(
        x0 + data.x['data'], y0 + data.y['data'], data.fields[field]['data'].max(axis=0),
        vmin=vmin, vmax=vmax, norm=norm, cmap=cmap, alpha=alpha)        

    if shapefile is not None:
        basemap.readshapefile(shapefile, 'shapefile', ax=ax)

    if title_flag:
        RadarMapDisplay(first_radar)._set_crf_title('composite reflectivity', title, ax)

    # add plot and field to lists
    RadarMapDisplay(first_radar).plots.append(pm)
    RadarMapDisplay(first_radar).plot_vars.append(field)

    if colorbar_flag:
        RadarMapDisplay(first_radar).plot_colorbar(
            mappable=pm, label=colorbar_label, field='composite reflectivity', fig=fig,
            ax=ax, ticks=ticks, ticklabs=ticklabs)
    return
