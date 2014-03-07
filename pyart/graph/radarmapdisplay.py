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

from .radardisplay import RadarDisplay

#Here we have an added method that inherits the radar
#display and allows the plotting on a basemap object

class RadarMapDisplay(RadarDisplay):
    try:
        from mpl_toolkits.basemap import Basemap, pyproj
    except ImportError:
        print("Requires mpl_toolkits.basemap please install")
        raise
    def plot_ppi_map(self, field, tilt, mask_tuple=None, vmin=None, vmax=None,
                 cmap='jet', mask_outside=True, title=None, title_flag=True,
                 axislabels=(None, None), axislabels_flag=True,
                 colorbar_flag=True, colorbar_label=None, ax=None, fig=None,
                 lat_lines =  np.arange(30,46,1), lon_lines =  np.arange(-110,-75,1),
                 box = None, resolution = 'h', shapefile = None):
        """
        Plot a PPI.

        Parameters
        ----------
        field : str
            Field to plot.
        tilt : int,
            Tilt number to plot.

        Other Parameters
        ----------------
        lat_lines : (array)
            Array contianing the latitude lines to plot
        lon_lines : (array)
            Array contianing the longitude lines to plot
        box : (list)
            [min_lon, max_lon, min_lat, max_lat] The bounds of plot, leave to None
            for the method to auto determine.
        shapefile : (string)
            filename for a esri shapefile as background (untested)
        resolution : (string)
            resolution parameter passed to basemap
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
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.

        """
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
        vmin, vmax = self._parse_vmin_vmax(field, vmin, vmax)

        # get data for the plot
        data = self._get_data(field, tilt, mask_tuple)
        x, y = self._get_x_y(field, tilt)

        # mask the data where outside the limits
        if mask_outside:
            data = np.ma.masked_outside(data, vmin, vmax)

        #build the proj
        pnyc = self.pyproj.Proj(proj='lcc',datum='NAD83',lat_0=self.loc[0],lon_0=self.loc[1],x_0=0.0,y_0=0.0)
        lon, lat = pnyc(x*1000.0, y*1000.0, inverse=True)

        #set up the map

        if box == None :
            #setup the map based on the ranges of lat and lons
            max_lat = lat.max()
            max_lon = lon.max()
            min_lat = lat.min()
            min_lon = lon.min()
        else:
            min_lon = box[0]
            max_lon = box[1]
            min_lat = box[2]
            max_lat = box[3]

        m = self.Basemap(llcrnrlon=min_lon,llcrnrlat=min_lat,
                    urcrnrlon=max_lon,urcrnrlat=max_lat,
                    projection='mill',area_thresh=1000,
                    resolution = resolution)

        xd,yd = m(lon, lat)
        print min_lat, max_lat

        # plot the data

        pm = m.pcolormesh(xd, yd, data, vmin=vmin, vmax=vmax, cmap=cmap)
        m.drawcoastlines(linewidth=1.25)
        m.drawstates()
        m.drawparallels(lat_lines,labels=[1,0,0,0])
        m.drawmeridians(lon_lines,labels=[0,0,0,1])
        if shapefile != None :
            shp_info = m.readshapefile(shapefile,'whatever',drawbounds=True)

        if title_flag:
            self._set_title(field, tilt, title, ax)

        # add plot and field to lists
        self.plots.append(pm)
        self.plot_vars.append(field)

        if colorbar_flag:
            self.plot_colorbar(mappable=pm, label=colorbar_label,
                               field=field, fig=fig)



