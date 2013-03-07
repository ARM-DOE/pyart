"""
pyart.graph.plot_mdv
====================

Routines for plotting radar data from MDV file.

.. autosummary::
    :toctree:: generated/

    MdvDisplay

"""

import getopt
import sys
import os

from pylab import *
import matplotlib.pyplot as plt
from numpy import ma, sin, cos, ones, zeros, pi, abs, sign
import numpy as np

from .common import dt_to_dict, corner_to_point, ax_radius
from .common_display import AbstractDisplay

class MdvDisplay(AbstractDisplay):
    """ 
    A class for creating plots from MdvFiles objects.
    """

    # XXX not generic
    def __init__(self, mdvfile):
        self.mdvfile = mdvfile
        info_dict = make_info(self.mdvfile, 'DBZ_F')
        lat, lon = info_dict['latitude_deg'], info_dict['longitude_deg']
        self._radar_location = [lat, lon]
        self._last_mappable = None

    # XXX not generic
    def plot_ppi(self, field, sweep, mask=None, vmin=None, vmax=None,
                 title_str=None, colorbar=False, fig=None, ax=None):
        """
        Plot a PPI.

        Parameters
        ----------
        field : str
            Field to plot.
        sweep : int,
            Sweep number to plot.
    
        Other Parameters
        ----------------
        mask : tuple
            2-Tuple containing the field name and value below which to mask
            field prior to plotting, for example to mask all data where 
            NCP < 0.5 set mask to ['NCP', 0.5]. None performs no masking.
        vmin : float
            Luminance minimum value, None for default value.
        vmax : float
            Luminance maximum value, None for default value.
        title_str : str
            String
        colorbar : bool
            True to add a colorbar with label to the axis.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.
        
        """
        # parse the parameters
        if ax is None:
            ax = plt.gca()
        if fig is None:
            fig = plt.gcf()

        info_dict = make_info(self.mdvfile, field)
        info_dict.update({'ele': self.mdvfile.el_deg[sweep]})
    
        if title_str is None:
            title_str = forminator()
    
        if vmin is None:
            vmin = get_default_range(self.mdvfile, field)[0]
        if vmax is None:
            vmax = get_default_range(self.mdvfile, field)[1]
     
        # extract coordinates in km
        x = self.mdvfile.carts['x'][sweep] / 1000.0  # x coords in km
        y = self.mdvfile.carts['y'][sweep] / 1000.0  # y coords in km
        z = self.mdvfile.carts['z'][sweep] / 1000.0  # z corrds in km
   
        # read in the data, mask if needed
        field_num = self.mdvfile.fields.index(field)
        data = self.mdvfile.read_a_field(field_num)[sweep]
        if mask is not None:
            mask_field, mask_value = mask
            mask_field_num = self.mdvfile.fields.index(mask_field)
            mdata = self.mdvfile.read_a_field(mask_field_num)[sweep]
            data = np.ma.masked_where(mdata < mask_value, data)

        # create the plot
        pm = ax.pcolormesh(x, y, data, vmin=vmin, vmax=vmax)
        self._last_mappable = pm
       
        # add labels and title
        ax.set_xlabel('x (km)')
        ax.set_ylabel('y (km)') 
        ax.set_title(title_str % info_dict)
 
        self._last_field = field
        self._last_field_dict = make_info(self.mdvfile, field)
        self._last_field_dict.update({'ele': self.mdvfile.el_deg[sweep]})
        
        if colorbar:
            self.plot_colorbar(pm, info_dict['name'], ax, fig)

    # XXX not generic
    def plot_rhi(self, field, sweep, mask=None, vmin=None, vmax=None,
                 title_str=None, colorbar=False, ax=None, fig=None):
        """ 
        Plot a RHI. 
        
        Parameters
        ----------
        field : str
            Field to plot.
        sweep : int,
            Sweep number to plot.
    
        Other Parameters
        ----------------
        mask : tuple
            2-Tuple containing the field name and value below which to mask
            field prior to plotting, for example to mask all data where 
            NCP < 0.5 set mask to ['NCP', 0.5]. None performs no masking.
        vmin : float
            Luminance minimum value, None for default value.
        vmax : float
            Luminance maximum value, None for default value.
        title_str : str
            String
        colorbar : bool
            True to add a colorbar with label to the axis.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.
 
        """
        # parse the parameters
        if ax is None:
            ax = plt.gca()
        if fig is None:
            fig = plt.gcf()

        info_dict = make_info(self.mdvfile, field)
        info_dict.update({'ele': self.mdvfile.az_deg[sweep]})

        if title_str is None:
            title_str = forminator()

        if vmin is None:
            vmin = get_default_range(self.mdvfile, field)[0]
        if vmax is None:
            vmax = get_default_range(self.mdvfile, field)[1]
 
        # extract coordinate in km
        x = self.mdvfile.carts['x'][sweep] / 1000.0
        y = self.mdvfile.carts['y'][sweep] / 1000.0
        z = self.mdvfile.carts['z'][sweep] / 1000.0
     
        # read in the data, mask if needed
        field_num = self.mdvfile.fields.index(field)
        data = self.mdvfile.read_a_field(field_num)[sweep]
        if mask is not None:
            mask_field, mask_value = mask
            mask_field_num = self.mdvfile.fields.index(mask_field)
            mdata = self.mdvfile.read_a_field(mask_field_num)[sweep]
            data = np.ma.masked_where(mdata < mask_value, data)

        # create the plot
        x_range = np.sign(y) * np.sqrt(y ** 2 + x ** 2)
        pm = ax.pcolormesh(x_range, z, data, vmin=vmin, vmax=vmax)
        
        self._last_mappable = pm
        
        # add labels and title
        ax.set_xlabel('Range (km)')
        ax.set_ylabel('Distance above radar (km)')
        ax.set_title(title_str % info_dict)

        self._last_field = field
        self._last_field_dict = make_info(self.mdvfile, field)
        self._last_field_dict.update({'ele': self.mdvfile.el_deg[sweep]})
        
        if colorbar:
            self.plot_colorbar(pm, info_dict['name'], ax, fig)


    # XXX not generic
    def plot_colorbar(self, mappable=None, label=None, ax=None, fig=None):
        """ 
        Plot a colorbar.
        
        Parameters
        ----------
        mappable : Image, ContourSet, etc.
            Image, ContourSet, etc to which the colorbar applied.  If None the
            last mappable object will be used.
        label :
            Colorbar label.  None will use a default value from the last field
            plotted.
        ax : Axis
            Axis to plot on.  None will use the current axis.
        fig : Figure
            Figure to place colorbar on.  None will use the current figure.
        
        """
        if ax is None:
            ax = plt.gca()
        if fig is None:
            fig = plt.gcf()
        if mappable is None:
            if self._last_mappable is None:
                raise ValueError("No mappable object found")
            mappable = self._last_mappable
            
        if label is None:
            label = self._last_field_dict['name']
        cb = fig.colorbar(mappable)
        cb.set_label(label)


def get_default_range(mdvfile, field):
    """ Return the default range for a field """
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
 

def forminator():
    """
    The forminator is used in future functions.. create a new forminator
    to format the titles of plots...
    """
    return '%(radar_name)s %(ele).1f Degree %(scan_type)s %(begin_year)04d-%(begin_month)02d-%(begin_day)02d %(begin_hour)02d:%(begin_minute)02d \n %(fancy_name)s '


def fancy_names():
    """
    Returns a dictionary for appending fancy names to a plot of MDV files...
    the moment names are typical of those from a TITAN setup
    """
    return {
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


def names_units():
    """
    Returns units for moments in an MDV file, moment names are typical of
    those in a file generated by TITAN
    """
    return {
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


def make_info(mdvobj, fld):
    """
    Concatinate metadata from the MDV file and info about the moment fld
    into one dictionary
    usage:
        dictionary = make_info(mdvobject,
                               string of the moment you want to append)

    """
    info = mdvobj.radar_info
    info.update({'scan_type': mdvobj.scan_type.upper()})
    info.update(dt_to_dict(mdvobj.times['time_begin'], pref='begin_'))
    info.update(dt_to_dict(mdvobj.times['time_end'], pref='end_'))
    name = names_units()[fld]
    units = dict([(mdvobj.fields[i], mdvobj.field_headers[i]['units'])
                  for i in range(len(mdvobj.fields))])[fld]
    fancy_name = fancy_names()[fld]
    info.update({'name': name, 'units': units, 'fancy_name': fancy_name})
    return info
