"""
pyart.graph.radardisplay
=========================

Class for creating plots from Radar objects.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    RadarDisplay

"""

import warnings

import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import numpy as np
import netCDF4

from . import common
from ..exceptions import DeprecatedAttribute
from ..core.transforms import antenna_to_cartesian
from ..core.transforms import antenna_vectors_to_cartesian
from ..core.transforms import geographic_to_cartesian_aeqd
from ..util.datetime_utils import datetimes_from_radar


class RadarDisplay(object):
    """
    A display object for creating plots from data in a radar object.

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
    antenna_transition : array or None
        Antenna transition flag (1 in transition, 0 in transition) or None
        if no antenna transition.

    """

    def __init__(self, radar, shift=(0.0, 0.0)):
        """ Initialize the object. """
        # save radar object
        self._radar = radar

        # populate attributes from radar object
        self.fields = radar.fields
        self.scan_type = radar.scan_type
        self.ranges = radar.range['data']
        self.azimuths = radar.azimuth['data']
        self.elevations = radar.elevation['data']
        self.fixed_angle = radar.fixed_angle['data']
        if radar.antenna_transition is None:
            self.antenna_transition = None
        else:
            self.antenna_transition = radar.antenna_transition['data']

        # origin
        if shift != (0.0, 0.0):
            self.origin = 'origin'
        else:
            self.origin = 'radar'

        self.shift = shift

        # radar location in latitude and longitude
        if radar.latitude['data'].size == 1:
            lat = float(radar.latitude['data'])
            lon = float(radar.longitude['data'])
        else:
            # for moving platforms stores use the median location.
            # The RadarDisplay object does not give a proper
            # visualization for moving platform data as the origin
            # of each ray changes and needs to be calculated individually or
            # georeferences.  When that is not available the following
            # gives acceptable results.
            lat = np.median(radar.latitude['data'])
            lon = np.median(radar.longitude['data'])
            warnings.warn('RadarDisplay does not correct for moving platforms')
        self.loc = (lat, lon)

        # list to hold plots, plotted fields and plotted colorbars
        self.plots = []
        self.plot_vars = []
        self.cbs = []

    @property
    def starts(self):
        """ Deprecated starts attribute. """
        warnings.warn(
            "The 'starts' attribute has been deprecated and will be removed"
            "in future versions of Py-ART", category=DeprecatedAttribute)
        return self._radar.sweep_start_ray_index['data']

    @property
    def ends(self):
        """ Deprecated starts attribute. """
        warnings.warn(
            "The 'ends' attribute has been deprecated and will be removed"
            "in future versions of Py-ART", category=DeprecatedAttribute)
        return self._radar.sweep_end_ray_index['data']

    @property
    def time_begin(self):
        """ Depeciated datetime object describing first sweep time. """
        warnings.warn(
            "The 'time_begin' attribute has been deprecated and will be "
            "removed in future versions of Py-ART",
            category=DeprecatedAttribute)
        times = self._radar.time['data'][0]
        units = self._radar.time['units']
        calendar = self._radar.time['calendar']
        return netCDF4.num2date(times, units, calendar)

    @property
    def radar_name(self):
        """ Deprecated radar_name attribute. """
        warnings.warn(
            "The 'radar_name' attribute has been deprecated and will be "
            "removed in future versions of Py-ART",
            category=DeprecatedAttribute)
        if 'instrument_name' in self._radar.metadata:
            return self._radar.metadata['instrument_name']
        else:
            return ''

    @property
    def x(self):
        """ Deprecated x coordinate attribute. """
        warnings.warn(
            "The 'x' attribute has been deprecated and will be removed in "
            "future versions of Py-ART", category=DeprecatedAttribute)
        rg, azg = np.meshgrid(self.ranges, self.azimuths)
        rg, eleg = np.meshgrid(self.ranges, self.elevations)
        return antenna_to_cartesian(rg / 1000.0, azg, eleg)[0] + self.shift[0]

    @property
    def y(self):
        """ Deprecated y coordinate attribute. """
        warnings.warn(
            "The 'z' attribute has been deprecated and will be removed in "
            "future versions of Py-ART", category=DeprecatedAttribute)
        rg, azg = np.meshgrid(self.ranges, self.azimuths)
        rg, eleg = np.meshgrid(self.ranges, self.elevations)
        return antenna_to_cartesian(rg / 1000.0, azg, eleg)[1] + self.shift[1]

    @property
    def z(self):
        """ Deprecated z coordinate attribute. """
        warnings.warn(
            "The 'z' attribute has been deprecated and will be removed in "
            "future versions of Py-ART", category=DeprecatedAttribute)
        rg, azg = np.meshgrid(self.ranges, self.azimuths)
        rg, eleg = np.meshgrid(self.ranges, self.elevations)
        return antenna_to_cartesian(rg / 1000.0, azg, eleg)[2]

    ####################
    # Plotting methods #
    ####################

    def plot(self, field, sweep=0, **kwargs):
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
        sweep : int
            Sweep number to plot, not used for VPT scans.

        See Also
        --------
        plot_ppi : Plot a PPI scan
        plot_rhi : Plot a RHI scan
        plot_vpt : Plot a VPT scan

        """
        if self.scan_type == 'ppi':
            self.plot_ppi(field, sweep, **kwargs)
        elif self.scan_type == 'rhi':
            self.plot_rhi(field, sweep, **kwargs)
        elif self.scan_type == 'vpt':
            self.plot_vpt(field, **kwargs)
        else:
            raise ValueError('unknown scan_type % s' % (self.scan_type))
        return

    def plot_ray(self, field, ray, format_str='k-', mask_tuple=None,
                 ray_min=None, ray_max=None, mask_outside=False, title=None,
                 title_flag=True, axislabels=(None, None), gatefilter=None,
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
        gatefilter : GateFilter
            GateFilter instance. None will result in no gatefilter mask being
            applied to data.
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
        ax, fig = common.parse_ax_fig(ax, fig)

        # get the data and mask
        data = self._get_ray_data(field, ray, mask_tuple, gatefilter)

        # mask the data where outside the limits
        _mask_outside(mask_outside, data, ray_min, ray_max)

        # plot the data
        line, = ax.plot(self.ranges / 1000., data, format_str)

        if title_flag:
            self._set_ray_title(field, ray, title, ax)

        if axislabels_flag:
            self._label_axes_ray(axislabels, field, ax)

        # add plot and field to attribute lists
        self.plots.append(line)
        self.plot_vars.append(field)

    def plot_ppi(
            self, field, sweep=0, mask_tuple=None,
            vmin=None, vmax=None, norm=None, cmap=None, mask_outside=False,
            title=None, title_flag=True,
            axislabels=(None, None), axislabels_flag=True,
            colorbar_flag=True, colorbar_label=None,
            colorbar_orient='vertical', edges=True, gatefilter=None,
            filter_transitions=True, ax=None, fig=None,
            ticks=None, ticklabs=None, **kwargs):
        """
        Plot a PPI.

        Additional arguments are passed to Matplotlib's pcolormesh function.

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
        colorbar_orient : 'vertical' or 'horizontal'
            Colorbar orientation.
        ticks : array
            Colorbar custom tick label locations.
        ticklabs : array
                Colorbar custom tick labels.
        edges : bool
            True will interpolate and extrapolate the gate edges from the
            range, azimuth and elevations in the radar, treating these
            as specifying the center of each gate.  False treats these
            coordinates themselved as the gate edges, resulting in a plot
            in which the last gate in each ray and the entire last ray are not
            plotted.
        gatefilter : GateFilter
            GateFilter instance. None will result in no gatefilter mask being
            applied to data.
        filter_transitions : bool
            True to remove rays where the antenna was in transition between
            sweeps from the plot.  False will include these rays in the plot.
            No rays are filtered when the antenna_transition attribute of the
            underlying radar is not present.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.

        """
        # parse parameters
        ax, fig = common.parse_ax_fig(ax, fig)
        norm, vmin, vmax = common.parse_norm_vmin_vmax(
            norm, self._radar, field, vmin, vmax)
        cmap = common.parse_cmap(cmap, field)

        # get data for the plot
        data = self._get_data(
            field, sweep, mask_tuple, filter_transitions, gatefilter)
        x, y = self._get_x_y(sweep, edges, filter_transitions)

        # mask the data where outside the limits
        _mask_outside(mask_outside, data, vmin, vmax)

        # plot the data
        pm = ax.pcolormesh(
            x, y, data, vmin=vmin, vmax=vmax, cmap=cmap, norm=norm, **kwargs)

        if title_flag:
            self._set_title(field, sweep, title, ax)

        if axislabels_flag:
            self._label_axes_ppi(axislabels, ax)

        # add plot and field to lists
        self.plots.append(pm)
        self.plot_vars.append(field)

        if colorbar_flag:
            self.plot_colorbar(
                mappable=pm, label=colorbar_label, orient=colorbar_orient,
                field=field, ax=ax, fig=fig, ticks=ticks, ticklabs=ticklabs)

    def plot_rhi(
            self, field, sweep=0, mask_tuple=None,
            vmin=None, vmax=None, norm=None, cmap=None,
            mask_outside=False, title=None, title_flag=True,
            axislabels=(None, None), axislabels_flag=True,
            reverse_xaxis=None, colorbar_flag=True, colorbar_label=None,
            colorbar_orient='vertical', edges=True, gatefilter=None,
            filter_transitions=True, ax=None, fig=None,
            ticks=None, ticklabs=None, **kwargs):
        """
        Plot a RHI.

        Additional arguments are passed to Matplotlib's pcolormesh function.

        Parameters
        ----------
        field : str
            Field to plot.
        sweep : int,
            Sweep number to plot.

        Other Parameters
        ----------------
        mask_tuple : (str, float)
            2-Tuple containing the field name and value below which to mask
            field prior to plotting, for example to mask all data where
            NCP < 0.5 set mask to ['NCP', 0.5]. None performs no masking.
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
        reverse_xaxis : bool or None
            True to reverse the x-axis so the plot reads east to west, False
            to have east to west.  None (the default) will reverse the axis
            only when all the distances are negative.
        colorbar_flag : bool
            True to add a colorbar with label to the axis.  False leaves off
            the colorbar.
        colorbar_label : str
            Colorbar label, None will use a default label generated from the
            field information.
        colorbar_orient : 'vertical' or 'horizontal'
            Colorbar orientation.
        ticks : array
            Colorbar custom tick label locations.
        ticklabs : array
                Colorbar custom tick labels.
        edges : bool
            True will interpolate and extrapolate the gate edges from the
            range, azimuth and elevations in the radar, treating these
            as specifying the center of each gate.  False treats these
            coordinates themselved as the gate edges, resulting in a plot
            in which the last gate in each ray and the entire last ray are not
            not plotted.
        gatefilter : GateFilter
            GateFilter instance. None will result in no gatefilter mask being
            applied to data.
        filter_transitions : bool
            True to remove rays where the antenna was in transition between
            sweeps from the plot.  False will include these rays in the plot.
            No rays are filtered when the antenna_transition attribute of the
            underlying radar is not present.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.

        """
        # parse parameters
        ax, fig = common.parse_ax_fig(ax, fig)
        norm, vmin, vmax = common.parse_norm_vmin_vmax(
            norm, self._radar, field, vmin, vmax)
        cmap = common.parse_cmap(cmap, field)

        # get data for the plot
        data = self._get_data(
            field, sweep, mask_tuple, filter_transitions, gatefilter)
        x, y, z = self._get_x_y_z(sweep, edges, filter_transitions)

        # mask the data where outside the limits
        _mask_outside(mask_outside, data, vmin, vmax)

        # plot the data
        R = np.sqrt(x ** 2 + y ** 2) * np.sign(y)
        if reverse_xaxis is None:
            # reverse if all distances are nearly negative (allow up to 1 m)
            reverse_xaxis = np.all(R < 1.)
        if reverse_xaxis:
            R = -R
        pm = ax.pcolormesh(
            R, z, data, vmin=vmin, vmax=vmax, cmap=cmap, norm=norm, **kwargs)

        if title_flag:
            self._set_title(field, sweep, title, ax)

        if axislabels_flag:
            self._label_axes_rhi(axislabels, ax)

        # add plot and field to lists
        self.plots.append(pm)
        self.plot_vars.append(field)

        if colorbar_flag:
            self.plot_colorbar(
                mappable=pm, label=colorbar_label, orient=colorbar_orient,
                field=field, ax=ax, fig=fig, ticks=ticks, ticklabs=ticklabs)

    def plot_vpt(
            self, field, mask_tuple=None,
            vmin=None, vmax=None, norm=None, cmap=None, mask_outside=False,
            title=None, title_flag=True,
            axislabels=(None, None), axislabels_flag=True,
            colorbar_flag=True, colorbar_label=None,
            colorbar_orient='vertical', edges=True,
            filter_transitions=True, time_axis_flag=False,
            date_time_form=None, tz=None, ax=None, fig=None,
            ticks=None, ticklabs=None, **kwargs):
        """
        Plot a VPT scan.

        Additional arguments are passed to Matplotlib's pcolormesh function.

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
        ticks : array
            Colorbar custom tick label locations.
        ticklabs : array
                Colorbar custom tick labels.
        colorbar_orient : 'vertical' or 'horizontal'
            Colorbar orientation.
        edges : bool
            True will interpolate and extrapolate the gate edges from the
            range, azimuth and elevations in the radar, treating these
            as specifying the center of each gate.  False treats these
            coordinates themselved as the gate edges, resulting in a plot
            in which the last gate in each ray and the entire last ray are not
            not plotted.
        filter_transitions : bool
            True to remove rays where the antenna was in transition between
            sweeps from the plot.  False will include these rays in the plot.
            No rays are filtered when the antenna_transition attribute of the
            underlying radar is not present.
        time_axis_flag : bool
            True to plot the x-axis as time. False uses the index number.
            Default is False - index-based.
        date_time_form : str, optional
            Format of the time string for x-axis labels. Parameter is
            ignored if time_axis_flag is set to False.
        tz : str, optional
            Time zone info to use when creating axis labels (see datetime).
            Parameter is ignored if time_axis_flag is set to False.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.

        """
        # parse parameters
        ax, fig = common.parse_ax_fig(ax, fig)
        norm, vmin, vmax = common.parse_norm_vmin_vmax(
            norm, self._radar, field, vmin, vmax)
        cmap = common.parse_cmap(cmap, field)

        # get data for the plot
        data = self._get_vpt_data(field, mask_tuple, filter_transitions)
        if edges:
            y = np.empty((self.ranges.shape[0] + 1, ), dtype=self.ranges.dtype)
            y[1:-1] = (self.ranges[:-1] + self.ranges[1:]) / 2.
            y[0] = self.ranges[0] - (self.ranges[1] - self.ranges[0]) / 2.
            y[-1] = self.ranges[-1] - (self.ranges[-2] - self.ranges[-1]) / 2.
            y[y < 0] = 0    # do not allow range to become negative
            y = y / 1000.
            x = np.arange(data.shape[1] + 1)
        else:
            x = np.arange(data.shape[1])
            y = self.ranges / 1000.

        # set up the time axis
        if time_axis_flag:
            self._set_vpt_time_axis(ax, date_time_form=date_time_form, tz=tz)
            x = datetimes_from_radar(self._radar)

        # mask the data where outside the limits
        _mask_outside(mask_outside, data, vmin, vmax)

        # plot the data
        pm = ax.pcolormesh(
            x, y, data, vmin=vmin, vmax=vmax, cmap=cmap, norm=norm, **kwargs)

        if title_flag:
            self._set_vpt_title(field, title, ax)

        if axislabels_flag:
            self._label_axes_vpt(axislabels, time_axis_flag, ax)

        # add plot and field to lists
        self.plots.append(pm)
        self.plot_vars.append(field)

        if colorbar_flag:
            self.plot_colorbar(
                mappable=pm, label=colorbar_label, orient=colorbar_orient,
                field=field, ax=ax, fig=fig, ticks=ticks, ticklabs=ticklabs)

    def plot_azimuth_to_rhi(
            self, field, target_azimuth, mask_tuple=None,
            vmin=None, vmax=None, norm=None, cmap=None, mask_outside=False,
            title=None, title_flag=True,
            axislabels=(None, None), axislabels_flag=True,
            colorbar_flag=True, colorbar_label=None,
            colorbar_orient='vertical', edges=True, gatefilter=None,
            reverse_xaxis=None, filter_transitions=True,
            ax=None, fig=None, ticks=None, ticklabs=None, **kwargs):
        """
        Plot pseudo-RHI scan by extracting the vertical field associated
        with the given azimuth.

        Additional arguments are passed to Matplotlib's pcolormesh function.

        Parameters
        ----------
        field : str
            Field to plot.
        target_azimuth : integer
            Azimuthal angle in degrees where cross section will be taken.

        Other Parameters
        ----------------
        mask_tuple : (str, float)
            2-Tuple containing the field name and value below which to mask
            field prior to plotting, for example to mask all data where
            NCP < 0.5 set mask to ['NCP', 0.5]. None performs no masking.
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
        reverse_xaxis : bool or None
            True to reverse the x-axis so the plot reads east to west, False
            to have east to west.  None (the default) will reverse the axis
            only when all the distances are negative.
        colorbar_flag : bool
            True to add a colorbar with label to the axis.  False leaves off
            the colorbar.
        colorbar_label : str
            Colorbar label, None will use a default label generated from the
            field information.
        ticks : array
            Colorbar custom tick label locations.
        ticklabs : array
                Colorbar custom tick labels.
        colorbar_orient : 'vertical' or 'horizontal'
            Colorbar orientation.
        edges : bool
            True will interpolate and extrapolate the gate edges from the
            range, azimuth and elevations in the radar, treating these
            as specifying the center of each gate.  False treats these
            coordinates themselved as the gate edges, resulting in a plot
            in which the last gate in each ray and the entire last ray are not
            not plotted.
        gatefilter : GateFilter
            GateFilter instance. None will result in no gatefilter mask being
            applied to data.
        filter_transitions : bool
            True to remove rays where the antenna was in transition between
            sweeps from the plot.  False will include these rays in the plot.
            No rays are filtered when the antenna_transition attribute of the
            underlying radar is not present.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.

        """
        # parse parameters
        ax, fig = common.parse_ax_fig(ax, fig)
        norm, vmin, vmax = common.parse_norm_vmin_vmax(
            norm, self._radar, field, vmin, vmax)
        cmap = common.parse_cmap(cmap, field)

        data, x, y, z = self._get_azimuth_rhi_data_x_y_z(
            field, target_azimuth, edges, mask_tuple,
            filter_transitions, gatefilter)

        # mask the data where outside the limits
        _mask_outside(mask_outside, data, vmin, vmax)

        # plot the data
        R = np.sqrt(x ** 2 + y ** 2) * np.sign(y)
        if reverse_xaxis is None:
            # reverse if all distances (nearly, up to 1 m) negative.
            reverse_xaxis = np.all(R < 1.)
        if reverse_xaxis:
            R = -R
        pm = ax.pcolormesh(
            R, z, data, vmin=vmin, vmax=vmax, cmap=cmap, norm=norm, **kwargs)

        if title_flag:
            self._set_az_rhi_title(field, target_azimuth, title, ax)

        if axislabels_flag:
            self._label_axes_rhi(axislabels, ax)

        # add plot and field to lists
        self.plots.append(pm)
        self.plot_vars.append(field)

        if colorbar_flag:
            self.plot_colorbar(
                mappable=pm, label=colorbar_label, orient=colorbar_orient,
                field=field, ax=ax, fig=fig, ticks=ticks, ticklabs=ticklabs)

    def plot_range_rings(self, range_rings, ax=None, col='k', ls='-', lw=2):
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
        for range_ring_location_km in range_rings:
            self.plot_range_ring(
                range_ring_location_km, ax=ax, col=col, ls=ls, lw=lw)

    @staticmethod
    def plot_range_ring(
            range_ring_location_km, npts=100, ax=None, col='k', ls='-', lw=2):
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
        ax = common.parse_ax(ax)
        theta = np.linspace(0, 2 * np.pi, npts)
        r = np.ones([npts], dtype=np.float32) * range_ring_location_km
        x = r * np.sin(theta)
        y = r * np.cos(theta)
        ax.plot(x, y, c=col, ls=ls, lw=lw)

    @staticmethod
    def plot_grid_lines(ax=None, col='k', ls=':'):
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
        ax = common.parse_ax(ax)
        ax.grid(c=col, ls=ls)

    def plot_labels(
            self, labels, locations, symbols='r+', text_color='k', ax=None):
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
        ax = common.parse_ax(ax)

        if type(symbols) is str:
            symbols = [symbols] * len(labels)
        if len(labels) != len(locations):
            raise ValueError('length of labels and locations must match')
        if len(labels) != len(symbols):
            raise ValueError('length of labels and symbols must match')

        for loc, label, sym in zip(locations, labels, symbols):
            self.plot_label(label, loc, sym, text_color, ax)

    def plot_label(
            self, label, location, symbol='r+', text_color='k', ax=None):
        """
        Plot a single symbol and label at a given location.

        Transforms of the symbol location in latitude and longitude units to
        x and y plot units is performed using an azimuthal equidistance
        map projection centered at the radar.

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
        ax = common.parse_ax(ax)
        location_lat, location_lon = location
        radar_lat, radar_lon = self.loc
        location_x, location_y = geographic_to_cartesian_aeqd(
            location_lon, location_lat, radar_lon, radar_lat)
        location_x /= 1000.0
        location_y /= 1000.0
        ax.plot([location_x], [location_y], symbol)
        ax.text(location_x - 5.0, location_y, label, color=text_color)

    @staticmethod
    def plot_cross_hair(size, npts=100, ax=None):
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
        ax = common.parse_ax(ax)
        x = np.zeros(npts, dtype=np.float32)
        y = np.linspace(-size, size, npts)
        ax.plot(x, y, 'k-')  # verticle
        ax.plot(y, x, 'k-')  # horizontal

    def plot_colorbar(self, mappable=None, field=None, label=None,
                      orient='vertical', cax=None, ax=None, fig=None,
                      ticks=None, ticklabs=None):
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
        ticks : array
            Colorbar custom tick label locations.
        ticklabs : array
                Colorbar custom tick labels.
        """
        if fig is None:
            fig = plt.gcf()
        if mappable is None:
            mappable = self.plots[-1]
        if label is None:
            if field is None:
                field = self.plot_vars[-1]
            label = self._get_colorbar_label(field)

        cb = fig.colorbar(mappable, orientation=orient, ax=ax, cax=cax)
        if ticks is not None:
            cb.set_ticks(ticks)
        if ticklabs:
            cb.set_ticklabels(ticklabs)
        cb.set_label(label)

        self.cbs.append(cb)

    ##########################
    # Plot adjusting methods #
    ##########################

    @staticmethod
    def set_limits(xlim=None, ylim=None, ax=None):
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
        common.set_limits(xlim, ylim, ax)

    def label_xaxis_x(self, ax=None):
        """ Label the xaxis with the default label for x units. """
        ax = common.parse_ax(ax)
        ax.set_xlabel('East West distance from ' + self.origin + ' (km)')

    def label_yaxis_y(self, ax=None):
        """ Label the yaxis with the default label for y units. """
        ax = common.parse_ax(ax)
        ax.set_ylabel('North South distance from ' + self.origin + ' (km)')

    def label_xaxis_r(self, ax=None):
        """ Label the xaxis with the default label for r units. """
        ax = common.parse_ax(ax)
        ax.set_xlabel('Distance from ' + self.origin + ' (km)')

    def label_yaxis_z(self, ax=None):
        """ Label the yaxis with the default label for z units. """
        ax = common.parse_ax(ax)
        ax.set_ylabel('Distance Above ' + self.origin + '  (km)')

    @staticmethod
    def label_xaxis_rays(ax=None):
        """ Label the yaxis with the default label for rays. """
        ax = common.parse_ax(ax)
        ax.set_xlabel('Ray number (unitless)')

    @staticmethod
    def label_xaxis_time(ax=None):
        """ Label the yaxis with the default label for rays. """
        ax = common.parse_ax(ax)
        ax.set_xlabel('Time (HH:MM)')

    def label_yaxis_field(self, field, ax=None):
        """ Label the yaxis with the default label for a field units. """
        ax = common.parse_ax(ax)
        ax.set_ylabel(self._get_colorbar_label(field))

    @staticmethod
    def set_aspect_ratio(aspect_ratio=0.75, ax=None):
        """ Set the aspect ratio for plot area. """
        ax = common.parse_ax(ax)
        ax.set_aspect(aspect_ratio)

    def _set_title(self, field, sweep, title, ax):
        """ Set the figure title using a default title. """
        if title is None:
            ax.set_title(self.generate_title(field, sweep))
        else:
            ax.set_title(title)

    def _set_vpt_title(self, field, title, ax):
        """ Set the figure title using a default title. """
        if title is None:
            ax.set_title(self.generate_vpt_title(field))
        else:
            ax.set_title(title)

    def _set_ray_title(self, field, ray, title, ax):
        """ Set the figure title for a ray plot using a default title. """
        if title is None:
            ax.set_title(self.generate_ray_title(field, ray))
        else:
            ax.set_title(title)

    def _set_az_rhi_title(self, field, azimuth, title, ax):
        """ Set the figure title for a ray plot using a default title. """
        if title is None:
            ax.set_title(self.generate_az_rhi_title(field, azimuth))
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

    def _label_axes_vpt(self, axis_labels, time_axis_flag, ax):
        """ Set the x and y axis labels for a PPI plot. """
        x_label, y_label = axis_labels
        if x_label is None:
            if time_axis_flag:
                self.label_xaxis_time(ax)
            else:
                self.label_xaxis_rays(ax)
        else:
            ax.set_xlabel(x_label)
        if y_label is None:
            self.label_yaxis_z(ax)
        else:
            ax.set_ylabel(y_label)

    @staticmethod
    def _set_vpt_time_axis(ax, date_time_form=None, tz=None):
        """ Set the x axis as a time formatted axis.

        Parameters
        ----------
        ax : Matplotlib axis instance
            Axis to plot. None will use the current axis.
        date_time_form : str
            Format of the time string for x-axis labels.
        tz : str
            Time zone info to use when creating axis labels (see datetime).

        """
        if date_time_form is None:
            date_time_form = '%H:%M'

        # Set the date format
        date_Fmt = DateFormatter(date_time_form, tz=tz)
        ax.xaxis.set_major_formatter(date_Fmt)

        # Turn the tick marks outward
        ax.tick_params(which='both', direction='out')

    ##########################
    # name generator methods #
    ##########################

    def generate_filename(self, field, sweep, ext='png'):
        """
        Generate a filename for a plot.

        Generated filename has form:
            radar_name_field_sweep_time.ext

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
        return common.generate_filename(self._radar, field, sweep, ext)

    def generate_title(self, field, sweep):
        """
        Generate a title for a plot.

        Parameters
        ----------
        field : str
            Field plotted.
        sweep : int
            Sweep plotted.

        Returns
        -------
        title : str
            Plot title.

        """
        return common.generate_title(self._radar, field, sweep)

    def generate_vpt_title(self, field):
        """
        Generate a title for a VPT plot.

        Parameters
        ----------
        field : str
            Field plotted.

        Returns
        -------
        title : str
            Plot title.

        """
        return common.generate_vpt_title(self._radar, field)

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
        return common.generate_ray_title(self._radar, field, ray)

    def generate_az_rhi_title(self, field, azimuth):
        """
        Generate a title for a ray plot.

        Parameters
        ----------
        field : str
            Field plotted.
        azimuth : float
            Azimuth plotted.

        Returns
        -------
        title : str
            Plot title.

        """
        return common.generate_az_rhi_title(self._radar, field, azimuth)

    ###############
    # Get methods #
    ###############

    def _get_data(self, field, sweep, mask_tuple, filter_transitions,
                  gatefilter):
        """ Retrieve and return data from a plot function. """
        sweep_slice = self._radar.get_slice(sweep)
        data = self.fields[field]['data'][sweep_slice]

        # mask data if mask_tuple provided
        if mask_tuple is not None:
            mask_field, mask_value = mask_tuple
            mdata = self.fields[mask_field]['data'][sweep_slice]
            data = np.ma.masked_where(mdata < mask_value, data)

        # mask data if gatefilter provided
        if gatefilter is not None:
            mask_filter = gatefilter.gate_excluded[sweep_slice]
            data = np.ma.masked_array(data, mask_filter)

        # filter out antenna transitions
        if filter_transitions and self.antenna_transition is not None:
            in_trans = self.antenna_transition[sweep_slice]
            data = data[in_trans == 0]

        return data

    def _get_vpt_data(self, field, mask_tuple, filter_transitions):
        """ Retrieve and return vpt data from a plot function. """
        data = self.fields[field]['data']

        # mask data if mask_tuple provided
        if mask_tuple is not None:
            mask_field, mask_value = mask_tuple
            mdata = self.fields[mask_field]['data']
            data = np.ma.masked_where(mdata < mask_value, data)

        # filter out antenna transitions
        if filter_transitions and self.antenna_transition is not None:
            in_trans = self.antenna_transition
            data = data[in_trans == 0]

        return data.T

    def _get_ray_data(self, field, ray, mask_tuple, gatefilter):
        """ Retrieve and return ray data from a plot function. """
        data = self.fields[field]['data'][ray]

        if mask_tuple is not None:
            mask_field, mask_value = mask_tuple
            mdata = self.fields[mask_field]['data'][ray]
            data = np.ma.masked_where(mdata < mask_value, data)

        # mask data if gatefilter provided
        if gatefilter is not None:
            mask_filter = gatefilter.gate_excluded[ray]
            data = np.ma.masked_array(data, mask_filter)

        return data

    def _get_azimuth_rhi_data_x_y_z(self, field, target_azimuth,
                                    edges, mask_tuple,
                                    filter_transitions, gatefilter):
        """Retrieve and return pseudo-RHI data from a plot function. """
        # determine which rays from the ppi radar make up the pseudo RHI
        data = self.fields[field]['data']

        if mask_tuple is not None:
            mask_field, mask_value = mask_tuple
            mdata = self.fields[mask_field]['data']
            data = np.ma.masked_where(mdata < mask_value, data)

        # mask data if gatefilter provided
        if gatefilter is not None:
            mask_filter = gatefilter.gate_excluded
            data = np.ma.masked_array(data, mask_filter)

        # filter out antenna transitions
        if filter_transitions and self.antenna_transition is not None:
            in_trans = self.antenna_transition
            data = data[in_trans == 0]

        prhi_rays = []
        for sweep_slice in self._radar.iter_slice():
            sweep_azimuths = self.azimuths[sweep_slice]
            ray_number = np.argmin(np.abs(sweep_azimuths - target_azimuth))
            prhi_rays.append(ray_number + sweep_slice.start)

        azimuth = self.azimuths[prhi_rays]
        elevation = self.elevations[prhi_rays]

        data = data[prhi_rays]

        x, y, z = antenna_vectors_to_cartesian(
            self.ranges, azimuth, elevation, edges=edges)
        x = (x + self.shift[0]) / 1000.0
        y = (y + self.shift[1]) / 1000.0
        z = z / 1000.0
        return data, x, y, z

    def _get_x_z(self, sweep, edges, filter_transitions):
        """ Retrieve and return x and z coordinate in km. """
        x, _, z = self._get_x_y_z(sweep, edges, filter_transitions)
        return x, z

    def _get_x_y(self, sweep, edges, filter_transitions):
        """ Retrieve and return x and y coordinate in km. """
        x, y, _ = self._get_x_y_z(sweep, edges, filter_transitions)
        return x, y

    def _get_x_y_z(self, sweep, edges, filter_transitions):
        """ Retrieve and return x, y, and z coordinate in km. """
        x, y, z = self._radar.get_gate_x_y_z(
            sweep, edges=edges, filter_transitions=filter_transitions)
        # add shift and convert to km
        x = (x + self.shift[0]) / 1000.0
        y = (y + self.shift[1]) / 1000.0
        z = z / 1000.0
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
        return common.generate_colorbar_label(standard_name, units)


def _mask_outside(flag, data, v1, v2):
    """ Return the data masked outside of v1 and v2 when flag is True.  """
    if flag:
        data = np.ma.masked_invalid(data)
        data = np.ma.masked_outside(data, v1, v2)
    return data
