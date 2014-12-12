"""
pyart.graph.radardisplay_airborne
=================================

Class for creating plots from Airborne Radar objects.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    RadarDisplay_Airborne

"""

import matplotlib.pyplot as plt
import numpy as np
import netCDF4
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .common import corner_to_point, radar_coords_to_cart
from .coord_transform import radar_coords_to_cart_track_relative
from .coord_transform import radar_coords_to_cart_earth_relative


class RadarDisplay_Airborne:
    """
    A display object for creating plots from data in a airborne radar object.

    Parameters
    ----------
    radar : Radar
        Radar object to use for creating plots, should be an airborne radar.
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
    rotation : array
        Rotation angle in degrees.
    roll : array
        Roll angle in degrees.
    drift : array
        Drift angle in degrees.
    tilt : array
        Tilt angle in degrees.
    heading : array
        Heading angle in degrees.
    pitch : array
        Pitch angle in degrees.
    altitude : array
        Altitude angle in meters.


    """
    # TODO this class could likely be made a subclass of RadarDisplay.

    def __init__(self, radar, shift=(0.0, 0.0)):
        """ Initialize the object. """
        # populate attributes from radar object
        self.fields = radar.fields
        self.scan_type = radar.scan_type
        self.ranges = radar.range['data']
        self.azimuths = radar.azimuth['data']
        self.elevations = radar.elevation['data']
        self.fixed_angle = radar.fixed_angle['data'][0]
        self.rotation = radar.rotation['data']
        self.roll = radar.roll['data']
        self.drift = radar.drift['data']
        self.tilt = radar.tilt['data']
        self.heading = radar.heading['data']
        self.pitch = radar.pitch['data']
        self.altitude = radar.altitude['data']
        if 'instrument_name' in radar.metadata:
            self.radar_name = radar.metadata['instrument_name']
        else:
            self.radar_name = ''

        # origin
        if shift != (0.0, 0.0):
            self.origin = 'origin'
        else:
            self.origin = 'radar'

        # x, y, z attributes: cartesian location for a sweep in km.
        self.shift = shift

        if radar.metadata['platform_type'] == 'aircraft_belly':
            rg, azg = np.meshgrid(self.ranges, self.azimuths)
            rg, eleg = np.meshgrid(self.ranges, self.elevations)
            self.x, self.y, self.z = radar_coords_to_cart(
                rg / 1000.0, azg, eleg)
            self.x = self.x + self.shift[0]
            self.y = self.y + self.shift[1]
        else:
            rg, rotg = np.meshgrid(self.ranges, self.rotation)
            rg, rollg = np.meshgrid(self.ranges, self.roll)
            rg, driftg = np.meshgrid(self.ranges, self.drift)
            rg, tiltg = np.meshgrid(self.ranges, self.tilt)
            rg, pitchg = np.meshgrid(self.ranges, self.pitch)
            self.x, self.y, self.z = radar_coords_to_cart_track_relative(
                rg / 1000.0, rotg, rollg, driftg, tiltg, pitchg)
            self.x = self.x + self.shift[0]
            self.y = self.y + self.shift[1]

        # radar location in latitude and longitude
        middle_lat = int(radar.latitude['data'].shape[0] / 2)
        middle_lon = int(radar.longitude['data'].shape[0] / 2)
        lat = float(radar.latitude['data'][middle_lat])
        lon = float(radar.longitude['data'][middle_lon])
        self.loc = (lat, lon)

        # datetime object describing first sweep time
        times = radar.time['data'][0]
        units = radar.time['units']
        calendar = radar.time['calendar']
        self.time_begin = netCDF4.num2date(times, units, calendar)

        # sweep start and end indices
        self.starts = radar.sweep_start_ray_index['data']
        self.ends = radar.sweep_end_ray_index['data']

        # list to hold plots, plotted fields and plotted colorbars
        self.plots = []
        self.plot_vars = []
        self.cbs = []

    ####################
    # Plotting methods #
    ####################

    def plot(self, field, **kwargs):
        """
        Create a plot appropiate for the radar.

        This function calls the plotting function corresponding to
        the scan_type of the radar.  Additional keywords can be passed to
        customize the plot, see the appropiate plot function for the
        allowed keywords.

        Parameters
        ----------
        field : str
            Field to plot.

        See Also
        --------
        plot_ppi : Plot a PPI scan
        plot_sweep_grid : Plot a RHI or VPT scan

        """
        if self.scan_type == 'ppi':
            self.plot_sweep_ppi(field, **kwargs)
        elif self.scan_type == 'rhi':
            self.plot_sweep_grid(field, **kwargs)
        elif self.scan_type == 'vpt':
            self.plot_sweep_grid(field, **kwargs)
        else:
            raise ValueError('unknown scan_type % s' % (self.scan_type))
        return

    def plot_ray(self, field, ray, format_str='k-', mask_tuple=None,
                 ray_min=None, ray_max=None, mask_outside=False, title=None,
                 title_flag=True, axislabels=(None, None),
                 axislabels_flag=True, ax=None, fig=None):
        """
        Plot a single ray.

        Parameters
        ----------
        field : str
            Field to plot.
        ray : int
            Ray number to plot.

        Other Parameters
        ----------------
        format_str : str
            Format string defining the line style and marker.
        mask_tuple : (str, float)
            Tuple containing the field name and value below which to mask
            field prior to plotting, for example to mask all data where
            NCP < 0.5 set mask_tuple to ['NCP', 0.5]. None performs no masking.
        ray_min : float
            Minimum ray value, None for default value, ignored if mask_outside
            is False.
        ray_max : float
            Maximum ray value, None for default value, ignored if mask_outside
            is False.
        mask_outside : bool
            True to mask data outside of vmin, vmax.  False performs no
            masking.
        title : str
            Title to label plot with, None to use default title generated from
            the field and ray parameters. Parameter is ignored if title_flag
            is False.
        title_flag : bool
            True to add a title to the plot, False does not add a title.
        axislabels : (str, str)
            2-tuple of x-axis, y-axis labels.  None for either label will use
            the default axis label.  Parameter is ignored if axislabels_flag is
            False.
        axislabel_flag : bool
            True to add label the axes, False does not label the axes.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.

        """
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)

        # get the data and mask
        data = self._get_ray_data(field, ray, mask_tuple)

         # mask the data where outside the limits
        if mask_outside:
            data = np.ma.masked_outside(data, ray_min, ray_max)

        # plot the data
        line, = ax.plot(self.ranges / 1000., data, format_str)

        if title_flag:
            self._set_ray_title(field, ray, title, ax)

        if axislabels_flag:
            self._label_axes_ray(axislabels, field, ax)

        # add plot and field to attribute lists
        self.plots.append(line)
        self.plot_vars.append(field)

    def plot_sweep_grid(self, field, mask_tuple=None, vmin=None, vmax=None,
                        cmap='jet', mask_outside=True, title=None,
                        title_flag=True, axislabels=(None, None),
                        axislabels_flag=True, colorbar_flag=True,
                        colorbar_label=None, colorbar_orient=None,
                        ax=None, fig=None):
        """
        Plot a sweep as a grid.

        Parameters
        ----------
        field : str
            Field to plot.

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
            the field and sweep parameters. Parameter is ignored if title_flag
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
        colorbar_orient : str
            Colorbar orientation, None will use default orientation of
            vertical.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.

        """
        # parse parameters
        ax, fig = self._parse_ax_fig(ax, fig)
        vmin, vmax = self._parse_vmin_vmax(field, vmin, vmax)

        # get data for the plot
        data = self._get_data(field, mask_tuple=mask_tuple)
        x, z = self._get_x_z(field)

        # mask the data where outside the limits
        if mask_outside:
            data = np.ma.masked_outside(data, vmin, vmax)

        # plot the data
        pm = ax.pcolormesh(x, z, data, vmin=vmin, vmax=vmax, cmap=cmap)

        # Set the aspcet ratio
        # ax.axis('scaled')

        if title_flag:
            self._set_title(field, title, ax)

        if axislabels_flag:
            self._label_axes_sweep(axislabels, ax)

        # add plot and field to lists
        self.plots.append(pm)
        self.plot_vars.append(field)

        # colorbar options
        if colorbar_flag:
            self.plot_colorbar(mappable=pm, label=colorbar_label,
                               orient=colorbar_orient,
                               field=field, ax=ax, fig=fig)

    def plot_ppi(self, field, sweep=0, mask_tuple=None, vmin=None, vmax=None,
                 cmap='jet', mask_outside=True, title=None, title_flag=True,
                 axislabels=(None, None), axislabels_flag=True,
                 colorbar_flag=True, colorbar_label=None, ax=None, fig=None):
        """
        Plot a PPI.

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
            the field and sweep parameters. Parameter is ignored if title_flag
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
        data = self._get_data(field, mask_tuple)
        x, y = self._get_x_y(field)

        # mask the data where outside the limits
        if mask_outside:
            data = np.ma.masked_outside(data, vmin, vmax)

        # plot the data
        pm = ax.pcolormesh(x, y, data, vmin=vmin, vmax=vmax, cmap=cmap)

        if title_flag:
            self._set_title(field, title, ax)

        if axislabels_flag:
            self._label_axes_ppi(axislabels, ax)

        # add plot and field to lists
        self.plots.append(pm)
        self.plot_vars.append(field)

        if colorbar_flag:
            self.plot_colorbar(mappable=pm, label=colorbar_label,
                               field=field, ax=ax, fig=fig)

    def plot_range_rings(self, range_rings, ax=None, col=None, ls=None,
                         lw=None):
        """
        Plot a series of range rings.

        Parameters
        ----------
        range_rings : list
            List of locations in km to draw range rings.
        ax : Axis
            Axis to plot on.  None will use the current axis.
        col : str or value
            Color to use for range rings.
        ls : str
            Linestyle to use for range rings.

        """
        ax = self._parse_ax(ax)
        for range_ring_location_km in range_rings:
            self.plot_range_ring(range_ring_location_km, ax=ax, col=col,
                                 ls=ls, lw=lw)

    def plot_range_ring(self, range_ring_location_km, npts=100, ax=None,
                        col=None, ls=None, lw=None):
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
        col : str or value
            Color to use for range rings.
        ls : str
            Linestyle to use for range rings.

        """
        ax = self._parse_ax(ax)
        theta = np.linspace(0, 2 * np.pi, npts)
        r = np.ones([npts], dtype=np.float32) * range_ring_location_km
        x = r * np.sin(theta)
        y = r * np.cos(theta)
        if lw is None:
            lw = 2
        if ls is None:
            ls = '-'
        if col is None:
            col = 'k'
        ax.plot(x, y, c=col, ls=ls, lw=lw)

    def plot_grid_lines(self, ax=None, col=None, ls=None, lw=None):
        """
        Plot grid lines.

        Parameters
        ----------
        ax : Axis
            Axis to plot on.  None will use the current axis.
        col : str or value
            Color to use for grid lines.
        ls : str
            Linestyle to use for grid lines.

        """
        ax = self._parse_ax(ax)
        if lw is None:
            lw = 2
        if ls is None:
            ls = ':'
        if col is None:
            col = 'k'
        ax.grid(c=col, ls=ls)

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

    def plot_colorbar(self, mappable=None, field=None, label=None, orient=None,
                      cax=None, ax=None, fig=None):
        """
        Plot a colorbar.

        Parameters
        ----------
        mappable : Image, ContourSet, etc.
            Image, ContourSet, etc to which the colorbar applied.  If None the
            last mappable object will be used.
        field : str
            Field to label colorbar with.
        label : str
            Colorbar label.  None will use a default value from the last field
            plotted.
        orient : str
            Colorbar orientation, either 'vertical' [default] or 'horizontal'.
        cax : Axis
            Axis onto which the colorbar will be drawn.  None is also valid.
        ax : Axes
            Axis onto which the colorbar will be drawn. None is also valid.
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
        # Find the axes locations to set colorbar
        box = make_axes_locatable(ax)

        if orient is None:
            orient = 'vertical'
        if orient == 'vertical':
            cax = box.append_axes("right", size="3%", pad=0.05)
        if orient == 'horizontal':
            cax = box.append_axes("bottom", size="6%", pad=0.50)
        cb = fig.colorbar(mappable, orientation=orient, ax=ax, cax=cax)
        cb.set_label(label)

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
        ax.set_xlabel('Horizontal distance from ' + self.origin + ' (km)')

    def label_yaxis_y(self, ax=None):
        """ Label the yaxis with the default label for y units. """
        ax = self._parse_ax(ax)
        ax.set_ylabel('Horizontal distance from ' + self.origin + ' (km)')

    def label_yaxis_z(self, ax=None):
        """ Label the yaxis with the default label for y units. """
        ax = self._parse_ax(ax)
        ax.set_ylabel('Altitude (km)')

    def label_xaxis_r(self, ax=None):
        """ Label the xaxis with the default label for r units. """
        ax = self._parse_ax(ax)
        ax.set_xlabel('Distance from ' + self.origin + ' (km)')

    def label_xaxis_rays(self, ax=None):
        """ Label the yaxis with the default label for rays. """
        ax = self._parse_ax(ax)
        ax.set_xlabel('Ray number (unitless)')

    def label_yaxis_field(self, field, ax=None):
        """ Label the yaxis with the default label for a field units. """
        ax = self._parse_ax(ax)
        ax.set_ylabel(self._get_colorbar_label(field))

    def set_aspect_ratio(self, aratio=None, ax=None):
        """ Set the aspect ratio for plot area. """
        ax = self._parse_ax(ax)
        if aratio is None:
            ax.set_aspect(0.75)
        else:
            ax.set_aspect(aratio)

    def _set_title(self, field, title, ax):
        """ Set the figure title using a default title. """
        if title is None:
            ax.set_title(self.generate_title(field))
        else:
            ax.set_title(title)

    def _set_ray_title(self, field, ray, title, ax):
        """ Set the figure title for a ray plot using a default title. """
        if title is None:
            ax.set_title(self.generate_ray_title(field, ray))
        else:
            ax.set_title(title)

    def _label_axes_sweep(self, axis_labels, ax):
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

    def _label_axes_ray(self, axis_labels, field, ax):
        """ Set the x and y axis labels for a ray plot. """
        x_label, y_label = axis_labels
        if x_label is None:
            self.label_xaxis_r(ax)
        else:
            ax.set_xlabel(x_label)
        if y_label is None:
            self.label_yaxis_field(field, ax)
        else:
            ax.set_ylabel(y_label)

    ##########################
    # name generator methods #
    ##########################

    def generate_filename(self, field, sweep, ext='png'):
        """
        Generate a filename for a plot.

        Generated filename has form:
            radar_name_field_time.ext

        Parameters
        ----------
        field : str
            Field plotted.
        sweep : int
            Sweep plotted.
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
        time_s = time_s.replace('-', ':')
        return '%s_%s_%s.%s' % (name_s, field_s, time_s, ext)

    def generate_title(self, field):
        """
        Generate a title for a plot.

        Parameters
        ----------
        field : str
            Field plotted.

        Returns
        -------
        title : str
            Plot title.

        """
        time_str = self.time_begin.isoformat() + 'Z'
        fixed_angle = self.fixed_angle
        l1 = "%s %.1f Deg. %s " % (self.radar_name, fixed_angle, time_str)
        field_name = self._generate_field_name(field)
        return l1 + '\n' + field_name

    def generate_ray_title(self, field, ray):
        """
        Generate a title for a ray plot.

        Parameters
        ----------
        field : str
            Field plotted.
        ray : int
            Ray plotted.

        Returns
        -------
        title : str
            Plot title.

        """
        time_str = self.time_begin.isoformat() + 'Z'
        l1 = "%s %s" % (self.radar_name, time_str)
        azim = self.azimuths[ray]
        elev = self.elevations[ray]
        l2 = "Ray: %i  Elevation: %.1f Azimuth: %.1f" % (ray, azim, elev)
        field_name = self._generate_field_name(field)
        return l1 + '\n' + l2 + '\n' + field_name

    def _generate_field_name(self, field):
        """ Return a nice field name for a particular field. """
        if 'standard_name' in self.fields[field]:
            field_name = self.fields[field]['standard_name']
        elif 'long_name' in self.fields[field]:
            field_name = self.fields[field]['long_name']
        else:
            field_name = str(field)
        field_name = field_name.replace('_', ' ')
        field_name = field_name[0].upper() + field_name[1:]
        return field_name

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

    def _get_data(self, field, mask_tuple):
        """ Retrieve and return data from a plot function. """
        data = self.fields[field]['data'][:]

        if mask_tuple is not None:  # mask data if mask_tuple provided
            mask_field, mask_value = mask_tuple
            mdata = self.fields[mask_field]['data'][:]
            data = np.ma.masked_where(mdata < mask_value, data)
        return data

    def _get_ray_data(self, field, ray, mask_tuple):
        """ Retrieve and return ray data from a plot function. """
        data = self.fields[field]['data'][ray]

        if mask_tuple is not None:
            mask_field, mask_value = mask_tuple
            mdata = self.fields[mask_field]['data'][ray]
            data = np.ma.masked_where(mdata < mask_value, data)
        return data

    def _get_x_z(self, field):
        """ Retrieve and return x and y coordinate in km. """
        return self.x/1000., self.z/1000.

    def _get_x_y(self, field):
        """ Retrieve and return x and y coordinate in km. """
        return self.x/1000.0, self.y/1000.0

    def _get_x_y_z(self, field):
        """ Retrieve and return x, y, and z coordinate in km. """
        x = self.x / 1000.0
        y = self.y / 1000.0
        z = self.z / 1000.0
        return x, y, z

    def _get_colorbar_label(self, field):
        """ Return a colorbar label for a given field. """
        last_field_dict = self.fields[field]
        if 'standard_name' in last_field_dict:
            standard_name = last_field_dict['standard_name']
        elif 'long_name' in last_field_dict:
            standard_name = last_field_dict['long_name']
        else:
            standard_name = field

        if 'units' in last_field_dict:
            units = last_field_dict['units']
        else:
            units = '?'
        return self._generate_colorbar_label(standard_name, units)
