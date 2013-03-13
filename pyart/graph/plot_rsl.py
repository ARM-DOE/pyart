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

from .common import radar_coords_to_cart, dms_to_d

from ..io import rsl, _rsl
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
        volume_num = _rsl.fieldTypes().list.index(field)
        sweep = self.rslradar.contents.volumes[volume_num].sweeps[tilt]
        fixed_angle = sweep.h.azimuth
        if fixed_angle == -999.0:   # ppi scan
            fixed_angle = sweep.h.elev

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
        data = _get_sweep_data(self.rslradar, field, tilt)

        if mask_tuple is not None:
            mask_field, mask_value = mask_tuple
            mdata = _get_sweep_data(self.rslradar, mask_field, tilt)
            data = np.ma.masked_where(mdata < mask_value, data)
        return data

    def _get_x_y(self, field, tilt):
        """ Retrieve and return x and y coordinate in km. """
        x, y, z = self._get_x_y_z(field, tilt)
        return x, y

    def _get_x_y_z(self, field, tilt):
        """ Retrieve and return x, y, and z coordinate in km. """
        volume_num = _rsl.fieldTypes().list.index(field)
        sweep = self.rslradar.contents.volumes[volume_num].sweeps[tilt]
        ranges = sweep.rays[0].dists / 1000.0
        elevs = [sweep.rays[i].h.elev for i in xrange(sweep.h.nrays)]
        azimuths = [sweep.rays[i].h.azimuth for i in xrange(sweep.h.nrays)]
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
