#Utilities to plot grid objects

import pyart
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pyproj
from netCDF4 import num2date


def mapplotgrid_3p(grid, fig, **kwargs):
    """
    Function to plot a three panel display from a grid object

    Parameters
    ----------
    grid : A pyart grid instance
    fig : a matplotlib figure object

    Keywords
    ----------
    time_step :
        the index of the time step to be plotted defaults to 0
    level :
        Index of the level to be potted in the map window defaults to 0
    field :
        Field to be plotted defaults to REF
    valid_min :
        Lower range of the colormesh. Tries to get from attributes,
        else defaults to -8
    valid_max :
        Lower range of the colormesh Tries to get from attributes,
        else defaults to 64
    lat_cut, lon_cut :
        Slicing locations defaults to grid origin
    min_lat, max_lat, min_lon, max_lon:
        ranges for the map window defaults to range determined from the data
    debug:
        verbose state, defaults to false
    overplot:
        set this to run more commands after the map is generated
        called as overplot(mapobject, gridobject)
        useful for using opendap et al to overplot analysis data

    """

    # keywords
    debug = kwargs.get('debug', False)

    time_step = kwargs.get('time_step', 0)
    level = kwargs.get('level', 0)
    lat_lines = kwargs.get('lat_lines', np.arange(30, 46, 1))
    lon_lines = kwargs.get('lon_lines', np.arange(-110, -75, 1))

    field = kwargs.get('field', 'REF')

    try:
        valid_min = kwargs.get('vmin', grid.fields[field]['valid_min'])
        valid_max = kwargs.get('vmax', grid.fields[field]['valid_max'])
    except KeyError:
        valid_min = -8
        valid_max = 64

    # lat and lon of the cut

    lat_cut = kwargs.get('lat_cut', grid.axes['lat']['data'][0])
    lon_cut = kwargs.get('lon_cut', grid.axes['lon']['data'][0])

    # set up projection

    lat0 = grid.axes['lat']['data'][0]
    lon0 = grid.axes['lon']['data'][0]
    pnyc = pyproj.Proj(proj='lcc', datum='NAD83', lat_0=lat0, lon_0=lon0,
                       x_0=0.0, y_0=0.0)

    # get x,y of cut
    x_cut, y_cut = pnyc(lon_cut, lat_cut)

    if debug:
        print "X cut: ", x_cut, " Y cut: ", y_cut

    # get element of cut
    element_x_cut = np.abs(grid.axes['x_disp']['data'] - x_cut).argsort()[0]
    element_y_cut = np.abs(grid.axes['y_disp']['data'] - y_cut).argsort()[0]

    if debug:
        print("X cut element: ", element_x_cut,
              "Y cut element: ", element_y_cut)

    # determine lat, lons of  grid points
    x_1d = grid.axes['x_disp']['data']
    y_1d = grid.axes['y_disp']['data']
    x_2d, y_2d = np.meshgrid(x_1d, y_1d)
    lon, lat = pnyc(x_2d, y_2d, inverse=True)

    # map ranges setup
    if kwargs.get('map_setup', 'auto') == 'auto':
        # setup the map based on the ranges of lat and lons
        max_lat = lat.max()
        max_lon = lon.max()
        min_lat = lat.min()
        min_lon = lon.min()
    else:
        min_lon = kwargs.get('min_lon', -92)
        max_lon = kwargs.get('max_lon', -86)
        min_lat = kwargs.get('min_lat', 40)
        max_lat = kwargs.get('max_lat', 44)

    if debug:
        print("Maximum latitude: ", max_lat,
              "Maximum longitude: ", max_lon,
              "Minimum latitude: ", min_lat,
              "Minimum longitute: ", min_lon)

    # set up the figure axes
    map_panel_axes = [0.05, 0.05, .4, .80]
    x_cut_panel_axes = [0.55, 0.10, .4, .30]
    y_cut_panel_axes = [0.55, 0.50, .4, .30]
    colorbar_panel_axes = [0.05, 0.90, .4, .03]

    ax1 = fig.add_axes(map_panel_axes)

    m = Basemap(llcrnrlon=min_lon, llcrnrlat=min_lat,
                urcrnrlon=max_lon, urcrnrlat=max_lat,
                projection='mill', area_thresh=10000,
                resolution='l')

    # Map projection for data
    xd, yd = m(lon, lat)

    #map projection for cut lines
    x_lon_cut, y_lon_cut = m(np.array([lon_cut, lon_cut]),
                             np.array([min_lat, max_lat]))
    x_lat_cut, y_lat_cut = m(np.array([min_lon, max_lon]),
                             np.array([lat_cut, lat_cut]))

    pcmesh_map = m.pcolormesh(xd, yd, grid.fields[field]['data'][level, :, :],
                              vmin=valid_min, vmax=valid_max)

    # if there is an overplotting command call this..
    # this can be used to, for example, overplot
    # reanalysis data
    if 'overplot' in kwargs.keys():
        if debug:
            print('overplotting')
        kwargs['overplot'](m, grid, **kwargs)

    plt.plot(x_lon_cut, y_lon_cut, 'r--', linewidth=2)
    plt.plot(x_lat_cut, y_lat_cut, 'r--', linewidth=2)

    # set up map navigation
    m.drawcoastlines(linewidth=1.25)
    m.drawstates()
    m.drawparallels(lat_lines, labels=[1, 1, 0, 0])
    m.drawmeridians(lon_lines, labels=[0, 0, 0, 1])

    # cut along constant lon, aka xcut
    ax2 = fig.add_axes(x_cut_panel_axes)
    pcolormesh_lon = plt.pcolormesh(
        grid.axes['y_disp']['data']/1000.0,
        grid.axes['z_disp']['data']/1000.0,
        grid.fields[field]['data'][:, :, element_x_cut],
        vmin=valid_min, vmax=valid_max)

    plt.xlabel('Distance from SGP CF (km)')
    plt.ylabel('Height (km)')
    plt.title('Slice at ' + str(lon_cut) + ' Longitude')

    # cut along constant lat, aka ycut
    ax3 = fig.add_axes(y_cut_panel_axes)
    pcolormesh_lat = plt.pcolormesh(
        grid.axes['x_disp']['data']/1000.0,
        grid.axes['z_disp']['data']/1000.0,
        grid.fields[field]['data'][:, element_y_cut, :],
        vmin=valid_min, vmax=valid_max)

    plt.ylabel('Height (km)')
    plt.title('Slice at ' + str(lat_cut) + ' Latitude')

    # Colorbar time
    cbax = fig.add_axes(colorbar_panel_axes)
    meshes_colorbar = plt.colorbar(cax=cbax, mappable=pcmesh_map,
                                   orientation='horizontal')
    meshes_colorbar.set_label(kwargs.get(
        'cb_label',
        grid.fields[field]['long_name'] + '(' + grid.fields[field]['units'] + ')'))

    # finally the time label
    slc_height = grid.axes['z_disp']['data'][level]

    dateobj = num2date(grid.axes['time']['data'], grid.axes['time']['units'])

    datestr = dateobj[0].strftime('%H:%M Z on %Y-%m-%d')

    fig.text(0.5, 0.9,
             'Sliced at ' + str(slc_height) + ' meters at ' + datestr)
