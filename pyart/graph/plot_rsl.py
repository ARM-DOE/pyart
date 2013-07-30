"""
pyart.graph.plot_rsl
====================

Routines for plotting radar data from files readable by RSL.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    RslDisplay

"""

# Nothing from this module is imported into pyart.graph is RSL is not
# installed.
import numpy as np

from .common import radar_coords_to_cart, dms_to_d

from ..io import rsl
from .radar_display import RadarDisplay


# XXX currenly broken
class RslDisplay(RadarDisplay):
    """
    A display object for creating plots from data in a RslFile object.

    Parameters
    ----------
    rslfile : RslFile
        RslFile object to use for creating plots.

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
    rslfile : RslFile
        RSL radar object use to create plots.

    """

    def __init__(self, rslfile):
        """ Initialize the object. """
        self.plots = []
        self.plot_vars = []
        self.cbs = []

        # radar location in lat, lon
        latd = rslfile.latd
        latm = rslfile.latm
        lats = rslfile.lats

        lond = rslfile.lond
        lonm = rslfile.lonm
        lons = rslfile.lons

        lat = dms_to_d((latd, latm, lats))
        lon = dms_to_d((lond, lonm, lons))
        self.loc = [lat, lon]
        self.origin = 'radar'
        dt = rslfile.get_volume(0).get_sweep(0).get_ray(0).get_datetime()
        self.time_begin = dt

        self.radar_name = rslfile.get_radar_header()['radar_name']
        self.rslfile = rslfile

    # public methods which are overridden.

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
        volume_num = rsl.RSLNAME2VOLUMENUM[field]
        sweep = self.rslfile.get_volume(volume_num).get_sweep(tilt)
        fixed_angle = sweep.azimuth
        if fixed_angle == -999.0:   # ppi scan
            fixed_angle = sweep.elev

        time_str = self.time_begin.isoformat() + 'Z'
        l1 = "%s %.1f Deg. %s " % (self.radar_name, fixed_angle, time_str)
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
        volume_num = rsl.RSLNAME2VOLUMENUM[field]
        sweep = self.rslfile.get_volume(volume_num).get_sweep(tilt)
        data = sweep.get_data()

        if mask_tuple is not None:
            mask_field, mask_value = mask_tuple
            mask_num = rsl.RSLNAME2VOLUMENUM[mask_field]
            sweep = self.rslfile.get_volume(mask_num).get_sweep(tilt)
            mdata = sweep.get_data()
            data = np.ma.masked_where(mdata < mask_value, data)
        return data

    def _get_x_y(self, field, tilt):
        """ Retrieve and return x and y coordinate in km. """
        x, y, z = self._get_x_y_z(field, tilt)
        return x, y

    def _get_x_y_z(self, field, tilt):
        """ Retrieve and return x, y, and z coordinate in km. """
        volume_num = rsl.RSLNAME2VOLUMENUM[field]
        sweep = self.rslfile.get_volume(volume_num).get_sweep(tilt)
        ray = sweep.get_ray(0)
        ranges = ray.range_bin1 + ray.gate_size * np.arange(ray.nbins)
        ranges = ranges / 1000.0
        elevs = [sweep.get_ray(i).elev for i in xrange(sweep.nrays)]
        azimuths = [sweep.get_ray(i).azimuth for i in xrange(sweep.nrays)]
        rg, ele = np.meshgrid(ranges, elevs)
        rg, azg = np.meshgrid(ranges, azimuths)
        x, y, z = radar_coords_to_cart(rg, azg, ele)
        return x / 1000.0, y / 1000.0, z / 1000.0

    def _get_colorbar_label(self, field):
        """ Return a colobar label for a given field. """
        if field in DEFAULT_COLORBAR_LABELS:
            return DEFAULT_COLORBAR_LABELS[field]
        else:
            return "Unknown units"


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
