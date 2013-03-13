"""
pyart.graph.radar_display
=========================

Class for creating plots from Radar objects.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    RadarDisplay

"""

import matplotlib.pyplot as plt
import numpy as np
import netCDF4

from .common import corner_to_point, radar_coords_to_cart


class RadarDisplay:
    """
    A display object for creating plots from data in a radar object.

    Parameters
    ----------
    radar : Radar
        Radar object to use for creating plots.
    shift : (float, float)
        Shifts in km to offset the calulcated x and y locations.

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
        Shift in km.
    x, y, z : array
        Cartesian location of a sweep in km.
    loc : (float, float)
        Latitude and Longitude of radar in degrees.
    time_begin : datetime
        Beginning time of first radar scan.
    starts : array
        Starting ray index for each tilt.
    ends : array
        Ending ray index for each tilt.
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

    """

    def __init__(self, radar, shift=(0.0, 0.0)):
        """ Initialize the object. """
        # populate attributes from radar object
        self.fields = radar.fields
        self.scan_type = radar.scan_type
        self.ranges = radar.range['data']
        self.azimuths = radar.azimuth['data']
        self.elevations = radar.elevation['data']
        self.fixed_angle = radar.sweep_info['fixed_angle']['data']
        self.radar_name = radar.metadata['instrument_name']

        # origin
        if shift != (0.0, 0.0):
            self.origin = 'origin'
        else:
            self.origin = 'radar'

        # x, y, z attributes: cartesian location for a sweep in km.
        rg, azg = np.meshgrid(self.ranges, self.azimuths)
        rg, eleg = np.meshgrid(self.ranges, self.elevations)
        self.shift = shift
        self.x, self.y, self.z = radar_coords_to_cart(rg / 1000.0, azg, eleg)
        self.x = self.x + self.shift[0]
        self.y = self.y + self.shift[1]

        # radar location in latitude and longitude
        lat = radar.location['latitude']['data']
        lon = radar.location['longitude']['data']
        self.loc = (lat, lon)

        # datetime object describing first sweep time
        times = radar.time['data'][0]
        units = radar.time['units']
        calendar = radar.time['calendar']
        self.time_begin = netCDF4.num2date(times, units, calendar)

        # sweep start and end indices
        self.starts = radar.sweep_info['sweep_start_ray_index']['data']
        self.ends = radar.sweep_info['sweep_end_ray_index']['data']

        # list to hold plots, plotted fields and plotted colorbars
        self.plots = []
        self.plot_vars = []
        self.cbs = []

    ####################
    # Plotting methods #
    ####################

    def plot_ppi(self, field, tilt, mask_tuple=None, vmin=None, vmax=None,
                 cmap='jet', mask_outside=True, title=None, title_flag=True,
                 axislabels=(None, None), axislabels_flag=True,
                 colorbar_flag=True, colorbar_label=None, ax=None, fig=None):
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

        # plot the data
        pm = ax.pcolormesh(x, y, data, vmin=vmin, vmax=vmax, cmap=cmap)

        if title_flag:
            self._set_title(field, tilt, title, ax)

        if axislabels_flag:
            self._label_axes_ppi(axislabels, ax)

        # add plot and field to lists
        self.plots.append(pm)
        self.plot_vars.append(field)

        if colorbar_flag:
            self.plot_colorbar(mappable=pm, label=colorbar_label,
                               field=field, fig=fig)

    def plot_rhi(self, field, tilt, mask_tuple=None, vmin=None, vmax=None,
                 cmap='jet', mask_outside=True, title=None, title_flag=True,
                 axislabels=(None, None), axislabels_flag=True,
                 colorbar_flag=True, colorbar_label= None, ax=None, fig=None):
        """
        Plot a RHI.

        Parameters
        ----------
        field : str
            Field to plot.
        tilt : int,
            Tilt number to plot.

        Other Parameters
        ----------------
        mask_tuple : (str, float)
            2-Tuple containing the field name and value below which to mask
            field prior to plotting, for example to mask all data where
            NCP < 0.5 set mask to ['NCP', 0.5]. None performs no masking.
        vmin : float
            Luminance minimum value, None for default value.
        vmax : float
            Luminance maximum value, None for default value.
        cmap : str
            Matplotlib colormap name.
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
        x, y, z = self._get_x_y_z(field, tilt)

        # mask the data where outside the limits
        if mask_outside:
            data = np.ma.masked_outside(data, vmin, vmax)

        # plot the data
        R = np.sqrt(x ** 2 + y ** 2) * np.sign(y)
        pm = ax.pcolormesh(R, z, data, vmin=vmin, vmax=vmax, cmap=cmap)

        if title_flag:
            self._set_title(field, tilt, title, ax)

        if axislabels_flag:
            self._label_axes_rhi(axislabels, ax)

        # add plot and field to lists
        self.plots.append(pm)
        self.plot_vars.append(field)

        if colorbar_flag:
            self.plot_colorbar(mappable=pm, label=colorbar_label,
                               field=field, fig=fig)

    def plot_range_rings(self, range_rings, ax=None):
        """
        Plot a series of range rings.

        Parameters
        ----------
        range_rings : list
            List of locations in km to draw range rings.
        ax : Axis
            Axis to plot on.  None will use the current axis.

        """
        ax = self._parse_ax(ax)
        for range_ring_location_km in range_rings:
            self.plot_range_ring(range_ring_location_km, ax=ax)

    def plot_range_ring(self, range_ring_location_km, npts=100, ax=None):
        """
        Plot a single range ring.

        Parameters
        ----------
        range_ring_location_km : float
            Location of range ring in km.
        npts: int
            Number of points in the ring, higher for better resolution.
        ax : Axis
            Axis to plot on.  None will use the current axis.

        """
        ax = self._parse_ax(ax)
        theta = np.linspace(0, 2 * np.pi, npts)
        r = np.ones([npts], dtype=np.float32) * range_ring_location_km
        x = r * np.sin(theta)
        y = r * np.cos(theta)
        ax.plot(x, y, 'k-')

    def plot_labels(self, labels, locations, symbols='r+', text_color='k',
                    ax=None):
        """
        Plot symbols and labels at given locations.

        Parameters
        ----------
        labels : list of str
            List of labels to place just above symbols.
        locations : list of 2-tuples
            List of latitude, longitude (in degrees) tuples at which symbols
            will be place.  Labels are placed just above the symbols.
        symbols : list of str or str
            List of matplotlib color+marker strings defining symbols to place
            at given locations.  If a single string is provided, that symbol
            will be placed at all locations.
        text_color : str
            Matplotlib color defining the color of the label text.
        ax : Axis
            Axis to plot on.  None will use the current axis.

        """
        ax = self._parse_ax(ax)

        if type(symbols) is str:
            symbols = [symbols] * len(labels)
        if len(labels) != len(locations):
            raise ValueError('length of labels and locations must match')
        if len(labels) != len(symbols):
            raise ValueError('length of labels and symbols must match')

        for loc, label, sym in zip(locations, labels, symbols):
            self.plot_label(label, loc, sym, text_color, ax)

    def plot_label(self, label, location, symbol='r+', text_color='k',
                   ax=None):
        """
        Plot a single symbol and label at a given location.

        Parameters
        ----------
        label : str
            Label text to place just above symbol.
        location : 2-tuples
            Tuple of latitude, longitude (in degrees) at which the symbol
            will be place.  The label is placed just above the symbol.
        symbol : str
            Matplotlib color+marker strings defining the symbol to place
            at the given location.
        text_color : str
            Matplotlib color defining the color of the label text.
        ax : Axis
            Axis to plot on.  None will use the current axis.

        """
        ax = self._parse_ax(ax)
        loc_x, loc_y = corner_to_point(self.loc, location)
        loc_x /= 1000.0
        loc_y /= 1000.0
        ax.plot([loc_x], [loc_y], symbol)
        ax.text(loc_x - 5.0, loc_y, label, color=text_color)

    def plot_cross_hair(self, size, npts=100, ax=None):
        """
        Plot a cross-hair on a ppi plot.

        Parameters
        ----------
        size : float
            Size of cross-hair in km.
        npts: int
            Number of points in the cross-hair, higher for better resolution.
        ax : Axis
            Axis to plot on.  None will use the current axis.

        """
        ax = self._parse_ax(ax)
        x = np.zeros(npts, dtype=np.float32)
        y = np.linspace(-size, size, npts)
        ax.plot(x, y, 'k-')  # verticle
        ax.plot(y, x, 'k-')  # horizontal

    def plot_colorbar(self, mappable=None, field=None, label=None, cax=None,
                      fig=None):
        """
        Plot a colorbar.

        Parameters
        ----------
        mappable : Image, ContourSet, etc.
            Image, ContourSet, etc to which the colorbar applied.  If None the
            last mappable object will be used.
        field : str
            Field to label colorbar with.
        label :
            Colorbar label.  None will use a default value from the last field
            plotted.
        cax : Axis
            Axis onto which the colorbar will be drawn.  None is also valid.
        fig : Figure
            Figure to place colorbar on.  None will use the current figure.

        """
        if fig is None:
            fig = plt.gcf()
        if mappable is None:
            mappable = self.plots[-1]
        if label is None:
            if field is None:
                field = self.plot_vars[-1]
            label = self._get_colorbar_label(field)
        cb = fig.colorbar(mappable, cax=cax)
        cb.set_label(label)
        self.cbs.append(cb)

    ##########################
    # Plot adjusting methods #
    ##########################

    def set_limits(self, xlim=None, ylim=None, ax=None):
        """
        Set the display limits.

        Parameters
        ----------
        xlim : tuple, optional
            2-Tuple containing y-axis limits in km. None uses default limits.
        ylim : tuple, optional
            2-Tuple containing x-axis limits in km. None uses default limits.
        ax : Axis
            Axis to adjust.  None will adjust the current axis.

        """
        if ax is None:
            ax = plt.gca()
        if ylim is not None:
            ax.set_ylim(ylim)
        if xlim is not None:
            ax.set_xlim(xlim)

    def label_xaxis_x(self, ax=None):
        """ Label the xaxis with the default label for x units. """
        ax = self._parse_ax(ax)
        ax.set_xlabel('East West distance from ' + self.origin + ' (km)')

    def label_yaxis_y(self, ax=None):
        """ Label the yaxis with the default label for y units. """
        ax = self._parse_ax(ax)
        ax.set_ylabel('North South distance from ' + self.origin + ' (km)')

    def label_xaxis_r(self, ax=None):
        """ Label the xaxis with the default label for r units. """
        ax = self._parse_ax(ax)
        ax.set_xlabel('Distance from ' + self.origin + ' (km)')

    def label_yaxis_z(self, ax=None):
        """ Label the yaxis with the default label for z units. """
        ax = self._parse_ax(ax)
        ax.set_ylabel('Distance Above ' + self.origin + '  (km)')

    def _set_title(self, field, tilt, title, ax):
        """ Set the figure title using a default title. """
        if title is None:
            ax.set_title(self.generate_title(field, tilt))
        else:
            ax.set_title(title)

    def _label_axes_ppi(self, axis_labels, ax):
        """ Set the x and y axis labels for a PPI plot. """
        x_label, y_label = axis_labels
        if x_label is None:
            self.label_xaxis_x(ax)
        else:
            ax.set_xlabel(x_label)
        if y_label is None:
            self.label_yaxis_y(ax)
        else:
            ax.set_ylabel(y_label)

    def _label_axes_rhi(self, axis_labels, ax):
        """ Set the x and y axis labels for a RHI plot. """
        x_label, y_label = axis_labels
        if x_label is None:
            self.label_xaxis_r(ax)
        else:
            ax.set_xlabel(x_label)
        if y_label is None:
            self.label_yaxis_z(ax)
        else:
            ax.set_ylabel(y_label)

    ##########################
    # name generator methods #
    ##########################

    def generate_filename(self, field, tilt, ext='png'):
        """
        Generate a filename for a plot.

        Generated filename has form:
            radar_name_field_tilt_time.ext

        Parameters
        ----------
        field : str
            Field plotted.
        titl : int
            Tilt plotted.
        ext : str
            Filename extension.

        Returns
        -------
        filename : str
            Filename suitable for saving a plot.

        """
        name_s = self.radar_name.replace(' ', '_')
        field_s = field.replace(' ', '_')
        time_s = self.time_begin.strftime('%Y%m%d%H%M%S')
        tilt_s = str(tilt).zfill(2)
        return '%s_%s_%s_%s.%s' % (name_s, field_s, tilt_s, time_s, ext)

    def generate_title(self, field, tilt):
        """
        Generate a title for a plot.

        Parameters
        ----------
        field : str
            Field plotted.
        tilt : int
            Tilt plotted.

        Returns
        -------
        title : str
            Plot title.

        """
        time_str = self.time_begin.isoformat() + 'Z'
        fixed_angle = self.fixed_angle[tilt]
        l1 = "%s %.1f Deg. %s " % (self.radar_name, fixed_angle, time_str)

        field_name = self.fields[field]['standard_name'].replace('_', ' ')
        field_name = field_name[0].upper() + field_name[1:]
        return l1 + '\n' + field_name

    def _generate_colorbar_label(self, standard_name, units):
        """ Generate and return a label for a colorbar. """
        return standard_name.replace('_', ' ') + ' (' + units + ')'

    ####################
    # Parseing methods #
    ####################

    def _parse_ax(self, ax):
        """ Parse and return ax parameter. """
        if ax is None:
            ax = plt.gca()
        return ax

    def _parse_ax_fig(self, ax, fig):
        """ Parse and return ax and fig parameters. """
        if ax is None:
            ax = plt.gca()
        if fig is None:
            fig = plt.gcf()
        return ax, fig

    def _parse_vmin_vmax(self, field, vmin, vmax):
        """ Parse and return vmin and vmax parameters. """
        field_dict = self.fields[field]
        if vmin is None:
            if 'valid_min' in field_dict:
                vmin = field_dict['valid_min']
            else:
                vmin = -6   # default value
        if vmax is None:
            if 'valid_max' in field_dict:
                vmax = field_dict['valid_max']
            else:
                vmax = 100
        return vmin, vmax

    ###############
    # Get methods #
    ###############

    def _get_data(self, field, tilt, mask_tuple):
        """ Retrieve and return data from a plot function. """
        start = self.starts[tilt]
        end = self.ends[tilt] + 1
        data = self.fields[field]['data'][start:end]

        if mask_tuple is not None:  # mask data if mask_tuple provided
            mask_field, mask_value = mask_tuple
            mdata = self.fields[mask_field]['data'][start:end]
            data = np.ma.masked_where(mdata < mask_value, data)
        return data

    def _get_x_y(self, field, tilt):
        """ Retrieve and return x and y coordinate in km. """
        start = self.starts[tilt]
        end = self.ends[tilt] + 1
        return self.x[start:end] / 1000.0, self.y[start:end] / 1000.0

    def _get_x_y_z(self, field, tilt):
        """ Retrieve and return x, y, and z coordinate in km. """
        start = self.starts[tilt]
        end = self.ends[tilt] + 1
        x = self.x[start:end] / 1000.0
        y = self.y[start:end] / 1000.0
        z = self.z[start:end] / 1000.0
        return x, y, z

    def _get_colorbar_label(self, field):
        """ Return a colorbar label for a given field. """
        last_field_dict = self.fields[field]
        if 'standard_name' in last_field_dict:
            standard_name = last_field_dict['standard_name']
        else:
            standard_name = last_field_dict['long_name']
        units = last_field_dict['units']
        return self._generate_colorbar_label(standard_name, units)
