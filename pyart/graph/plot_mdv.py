"""
pyart.graph.plot_mdv
====================

Routines for plotting radar data from MDV file.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    MdvDisplay

.. autosummary::
    :toctree: generated/

    _get_default_range

"""

import warnings

import numpy as np

from .radardisplay import RadarDisplay


class MdvDisplay(RadarDisplay):
    """
    A display object for creating plots from data in a MdvFile objects.

    This class has been deprecated and will be removed in the next release of
    Py-ART.

    Parameters
    ----------
    mdvfile : MdvFile
        MdvFile object to use for creating plots.

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
    mdvfile : MdvFile
        MdvFile object use to create plots.

    """

    def __init__(self, mdvfile):
        """ Initialize the object. """
        warnings.warn(
            ('This class has been deprecated and will be removed in the '
             'next release of Py-ART'), DeprecationWarning)
        self.plots = []
        self.plot_vars = []
        self.cbs = []
        lat = mdvfile.radar_info['latitude_deg']
        lon = mdvfile.radar_info['longitude_deg']
        self.loc = [lat, lon]
        self.origin = 'radar'
        self.time_begin = mdvfile.times['time_begin']
        self.radar_name = mdvfile.master_header['data_set_source']
        self.mdvfile = mdvfile

    # public methods which are overridden.

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
        az_deg, range_km, el_deg = self.mdvfile._calc_geometry()
        if self.mdvfile.projection == 'rhi':
            fixed_angle = az_deg[sweep]
        else:
            fixed_angle = el_deg[sweep]

        time_str = self.time_begin.isoformat() + 'Z'
        l1 = "%s %.1f Deg. %s " % (self.radar_name, fixed_angle, time_str)
        field_name = FANCY_NAMES[field]
        return l1 + '\n' + field_name

    # private methods which are overridden

    def _parse_vmin_vmax(self, field, vmin, vmax):
        if vmin is None:
            vmin = _get_default_range(self.mdvfile, field)[0]
        if vmax is None:
            vmax = _get_default_range(self.mdvfile, field)[1]
        return vmin, vmax

    def _get_data(self, field, sweep, mask_tuple, filter_transitions):
        """ Retrieve and return data from a plot function. """
        field_num = self.mdvfile.fields.index(field)
        data = self.mdvfile.read_a_field(field_num)[sweep]
        if mask_tuple is not None:
            mask_field, mask_value = mask_tuple
            mask_field_num = self.mdvfile.fields.index(mask_field)
            mdata = self.mdvfile.read_a_field(mask_field_num)[sweep]
            data = np.ma.masked_where(np.ma.less(mdata, mask_value), data)
        return data

    def _get_x_y(self, field, sweep, edges, filter_transitions):
        """ Retrieve and return x and y coordinate in km. """
        # TODO perform interpolating if edges is True
        carts = self.mdvfile._make_carts_dict()
        x = carts['x'][sweep] / 1000.0  # x coords in km
        y = carts['y'][sweep] / 1000.0  # y coords in km
        return x, y

    def _get_x_y_z(self, field, sweep, edges, filter_transitions):
        """ Retrieve and return x, y, and z coordinate in km. """
        # TODO perform interpolating if edges is True
        carts = self.mdvfile._make_carts_dict()
        x = carts['x'][sweep] / 1000.0  # x coords in km
        y = carts['y'][sweep] / 1000.0  # y coords in km
        z = carts['z'][sweep] / 1000.0  # y coords in km
        return x, y, z

    def _get_colorbar_label(self, field):
        """ Return a colobar label for a given field. """
        return UNIT_NAMES[field]


def _get_default_range(mdvfile, field):
    """ Return the default range for a field. """
    def_ranges = {
        'DBMHC': [-100, 0],
        'DBMVC': [-100, 0],
        'DBZ': [-16.0, 64.0],
        'DBZ_F': [-16.0, 64.0],
        'DBZVC': [-16.0, 64.0],
        'DBZVC_F': [-16.0, 64.0],
        'VEL': [-1.0*mdvfile.radar_info['unambig_vel_mps'],
                mdvfile.radar_info['unambig_vel_mps']],
        'VEL_F': [-1.0*mdvfile.radar_info['unambig_vel_mps'],
                  mdvfile.radar_info['unambig_vel_mps']],
        'WIDTH': [0.0, 10.0],
        'WIDTH_F': [0.0, 10.0],
        'ZDR': [-3, 6.0],
        'ZDR_F': [-3, 6.0],
        'RHOHV': [0.6, 1.0],
        'RHOHV_F': [0.6, 1.0],
        'PHIDP': [0, 180.0],
        'PHIDP_F': [-180.0, 180.0],
        'KDP': [-2, 6],
        'KDP_F': [-2, 6],
        'NCP': [0, 1],
        'NCP_F': [0, 1]}
    return def_ranges[field]


FANCY_NAMES = {
    'DBMHC': "Horizontal recieved power",
    'DBMVC': "Vertical recieved power",
    'DBZ': 'Horizontal equivalent reflectivity factor',
    'DBZ_F': 'Horizontal equivalent reflectivity factor',
    'DBZVC': 'Vertical equivalent reflectivity factor',
    'DBZVC_F': 'Vertical equivalent reflectivity factor',
    'VEL': "Radial velocity of scatterers (positive away)",
    'VEL_F': "Radial velocity of scatterers (positive away)",
    'WIDTH': "Spectral Width",
    'WIDTH_F': "Spectral Width",
    'ZDR': "Differential reflectivity",
    'ZDR_F': "Differential reflectivity",
    'RHOHV': 'Co-Polar correlation coefficient',
    'RHOHV_F': 'Co-Polar Correlation Coefficient',
    'PHIDP': "Differential propigation phase",
    'PHIDP_F': "Differential propigation phase",
    'KDP': "Specific differential phase",
    'KDP_F': "Specific differential phase",
    'NCP': "Normalized coherent power",
    'NCP_F': "Normalized coherent power"}


UNIT_NAMES = {
    'DBMHC': " H Rec Power (dBm)",
    'DBMVC': "V Rec Power (dBm)",
    'DBZ': "Hz Eq. Ref. Fac (dBz)",
    'DBZ_F': "Hz Eq. Ref. Fac (dBz)",
    'DBZVC': "V Eq. Ref. Fac (dBz)",
    'DBZVC_F': "V Eq. Ref. Fac (dBz)",
    'VEL': "Rad. Vel. (m/s, +away)",
    'VEL_F': "Rad. Vel. (m/s, +away)",
    'WIDTH': "Spec. Width (m/s)",
    'WIDTH_F': "Spec. Width (m/s)",
    'ZDR': "Dif Refl (dB)",
    'ZDR_F': "Dif Refl (dB)",
    'RHOHV': "Cor. Coef (frac)",
    'RHOHV_F': "Cor. Coef (frac)",
    'PHIDP': "Dif Phase (deg)",
    'PHIDP_F': "Dif Phase (deg)",
    'KDP': "Spec Dif Ph. (deg/km)",
    'KDP_F': "Spec Dif Ph. (deg/km)",
    'NCP': "Norm. Coh. Power (frac)",
    'NCP_F': "Norm. Coh. Power (frac)"}
