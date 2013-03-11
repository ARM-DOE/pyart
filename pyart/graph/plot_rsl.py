"""
pyart.graph.plot_rsl
====================

Routines for plotting radar data from files readable by RSL.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    RslDisplay
    _get_sweep_data

"""

import numpy as np

import pyart.io._rsl as _rsl
from .common import radar_coords_to_cart, dms_to_d

from ..io import rsl
from .radar_display import RadarDisplay


class RslDisplay(RadarDisplay):
    """
    A display object for creating plots from data in a RSL Radar objects.

    Parameters
    ----------
    rslradar : RSL Radar
        RSL Radar object to use for creating plots.

    Attributes
    ----------
    plots : list
        List of plots created.
    plot_vars : list
        List of fields plotted, order matches plot list.
    cbs : list
        List of colorbars created.
    loc : (float, float)
        Latitude and Longitude of radar in degrees.
    origin : 'Radar'
       'Radar'
    time_begin : datetime
        Time of first radar scan.
    radar_name : str
        Radar name.
    rslradar : RSL radar
        RSL radar object use to create plots.

    """

    def __init__(self, rslradar):
        """ Initialize the object. """
        self.plots = []
        self.plot_vars = []
        self.cbs = []

        latd = rslradar.contents.h.latd
        latm = rslradar.contents.h.latm
        lats = rslradar.contents.h.lats

        lond = rslradar.contents.h.lond
        lonm = rslradar.contents.h.lonm
        lons = rslradar.contents.h.lons

        lat = dms_to_d((latd, latm, lats))
        lon = dms_to_d((lond, lonm, lons))
        self.loc = [lat, lon]
        self.origin = 'radar'

        first_ray_header = rslradar.contents.volumes[0].sweeps[0].rays[0].h
        self.time_begin = rsl.ray_header_time_to_datetime(first_ray_header)

        self.radar_name = rslradar.contents.h.radar_name
        self.rslradar = rslradar

    # public methods which are overridden.

    def plot_ppi(self, field, tilt, mask_tuple=None, vmin=None, vmax=None,
                 cmap='jet', title=None, title_flag=True,
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
        x, y = self._get_x_y(field, tilt)   # added field

        # plot the data
        data = np.ma.masked_outside(data, vmin, vmax)   # added masking
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

    # added masking for data outside vmin, vmax
    def plot_rhi(self, field, tilt, mask_tuple=None, vmin=None, vmax=None,
                 cmap='jet', title=None, title_flag=True,
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
        x, y, z = self._get_x_y_z(field, tilt)  # added field

        # plot the data
        R = np.sqrt(x ** 2 + y ** 2) * np.sign(y)   # added masking
        data = np.ma.masked_outside(data, vmin, vmax)
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
        radar_name = self.radar_name
        volume_num = _rsl.fieldTypes().list.index(field)
        sweep = self.rslradar.contents.volumes[volume_num].sweep[tilt]
        fixed_angle = sweep.rays[0].azimuth

        time_str = self.time_begin.isoformat() + 'Z'
        l1 = "%s %.1f Deg. %s " % (radar_name, fixed_angle, time_str)
        return l1 + '\n' + field

    # private methods which are overridden

    def _parse_vmin_vmax(self, field, vmin, vmax):
        if vmin is None:
            if field in DEFAULT_VMIN_VMAX:
                vmin = DEFAULT_VMIN_VMAX[field][0]
            else:
                vmin = -16.0
        if vmax is None:
            if field in DEFAULT_VMIN_VMAX:
                vmax = DEFAULT_VMIN_VMAX[field][1]
            else:
                vmax = 64.0
        return vmin, vmax

    def _get_data(self, field, tilt, mask_tuple):
        """ Retrieve and return data from a plot function. """
        data = _get_sweep_data(self.rslradar, field, tilt)

        if mask_tuple is not None:
            mask_field, mask_value = mask_tuple
            mdata = _get_sweep_data(self.rslradar, mask_field, tilt)
            data = np.ma.masked_where(mdata < mask_value, data)
        return data

    def _get_x_y(self, field, tilt):
        """ Retrieve and return x and y coordinate in km. """
        volume_num = _rsl.fieldTypes().list.index(field)
        sweep = self.rslradar.contents.volumes[volume_num].sweeps[tilt]

        ranges = sweep.rays[0].dists / 1000.0
        elevs = [sweep.rays[i].h.elev for i in xrange(sweep.h.nrays)]
        rg, ele = np.meshgrid(ranges, elevs)
        azg = np.ones_like(rg) * sweep.h.azimuth
        x, y, z = radar_coords_to_cart(rg, azg, ele)
        return x / 1000.0, y / 1000.0

    def _get_x_y_z(self, field, tilt):
        """ Retrieve and return x, y, and z coordinate in km. """
        volume_num = _rsl.fieldTypes().list.index(field)
        sweep = self.rslradar.contents.volumes[volume_num].sweeps[tilt]

        ranges = sweep.rays[0].dists / 1000.0
        elevs = [sweep.rays[i].h.elev for i in xrange(sweep.h.nrays)]
        rg, ele = np.meshgrid(ranges, elevs)
        azg = np.ones_like(rg) * sweep.h.azimuth
        x, y, z = radar_coords_to_cart(rg, azg, ele)
        return x / 1000.0, y / 1000.0, z / 1000.0

    def _get_colorbar_label(self, field):
        """ Return a colobar label for a given field. """
        if field in DEFAULT_COLORBAR_LABELS:
            return DEFAULT_COLORBAR_LABELS[field]
        else:
            return "Unknown units"


def _get_sweep_data(rslradar, field, tilt):
    """ Extract and return the data for a given field and tilt (sweep). """
    volume_num = _rsl.fieldTypes().list.index(field)
    sweep = rslradar.contents.volumes[volume_num].sweeps[tilt]
    nrays = sweep.h.nrays
    nbins = sweep.rays[0].h.nbins
    data = np.zeros([nrays, nbins], 'float32') + 1.31072000e+05
    for raynum in xrange(nrays):
        ray_data = sweep.rays[raynum].data
        data[raynum, 0:len(ray_data)] = ray_data
    return data

DEFAULT_VMIN_VMAX = {
    'CZ': [-16., 64.],
    'DZ': [-16., 64.],
    'ZT': [-16.0, 64.0],
    'VE': [-40.0, 40.0],
    'VR': [-16, 16],
    'PH': [300, 350],
    'DR': [-4, 3]}

DEFAULT_COLORBAR_LABELS = {
    'CZ': 'Eq refl fact (dBz)',
    'DZ': 'Eq refl fact (dBz)',
    'ZT': 'Eq refl fact (dBz)',
    'PH': 'Phidp deg',
    'VE': 'Unfolded VR (m/s)',
    'DR': 'ZDR dB',
    'VR': 'Radial Velocity (+away) m/s'}
