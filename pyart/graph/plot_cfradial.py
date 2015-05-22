"""
pyart.graph.plot_cfradial
=========================

Routines for plotting radar data from CF/Radial netCDF files.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    CFRadialDisplay

"""

import warnings

import numpy as np
import netCDF4

from .radardisplay import RadarDisplay
from .common import radar_coords_to_cart


class CFRadialDisplay(RadarDisplay):
    """
    A display object for creating plots from data in NetCDF4 Dataset objects.

    This class has been deprecated and will be removed in the next release of
    Py-ART.

    Parameters
    ----------
    dataset : Dataset
        NetCDF4 Dataset object to use for creating plots.
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
    ranges : array
        Gate ranges in meters.
    azimuths : array
        Azimuth angle in degrees.
    elevations : array
        Elevations in degrees.
    fixed_angle : array
        Scan angle in degrees.

    """
    # missing scan_type and fields attributes
    def __init__(self, dataset, shift=(0.0, 0.0)):
        """ Initialize the object. """
        warnings.warn(
            ('This class has been deprecated and will be removed in the '
             'next release of Py-ART'), DeprecationWarning)

        # populate attributes from dataset
        self.dataset = dataset
        self.ranges = dataset.variables['range']
        self.azimuths = dataset.variables['azimuth']
        self.elevations = dataset.variables['elevation']
        self.fixed_angle = dataset.variables['fixed_angle']
        self.antenna_transition = None
        if 'instrument_name' in dataset.ncattrs():
            self.radar_name = dataset.instrument_name
        else:
            self.radar_name = ''

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
        lon = dataset.variables['longitude'][:]
        lat = dataset.variables['latitude'][:]
        self.loc = [lat, lon]

        # datetime object describing first sweep time
        t = dataset.variables['time']
        self.time_begin = netCDF4.num2date(t[0], t.units, t.calendar)

        # sweep start and end indices
        self.starts = dataset.variables['sweep_start_ray_index'][:]
        self.ends = dataset.variables['sweep_end_ray_index'][:]

        # list to hold plots, plotted fields and plotted colorbars
        self.plots = []
        self.plot_vars = []
        self.cbs = []

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
        time_str = self.time_begin.isoformat() + 'Z'
        fixed_angle = self.fixed_angle[sweep]
        l1 = "%s %.1f Deg. %s " % (self.radar_name, fixed_angle, time_str)

        field_name = self.dataset.variables[field].standard_name
        field_name = field_name.replace('_', ' ')
        field_name = field_name[0].upper() + field_name[1:]
        return l1 + '\n' + field_name

    # private methods which are overridden

    def _parse_vmin_vmax(self, field, vmin, vmax):
        if vmin is None:
            vmin = self.dataset.variables[field].valid_min
        if vmax is None:
            vmax = self.dataset.variables[field].valid_max
        return vmin, vmax

    def _get_data(self, field, sweep, mask_tuple, filter_transitions):
        """ Retrieve and return data from a plot function. """
        start = self.starts[sweep]
        end = self.ends[sweep] + 1
        data = self.dataset.variables[field][start:end]

        if mask_tuple is not None:
            mask_field, mask_value = mask_tuple
            mdata = self.dataset.variables[mask_field][start:end]
            data = np.ma.masked_where(mdata < mask_value, data)
        return data

    def _get_colorbar_label(self, field):
        """ Return a colobar label for a given field. """
        last_field_dict = self.dataset.variables[field].ncattrs()
        if 'standard_name' in last_field_dict:
            standard_name = self.dataset.variables[field].standard_name
        else:
            standard_name = self.dataset.variables[field].long_name
        units = self.dataset.variables[field].units
        return self._generate_colorbar_label(standard_name, units)
