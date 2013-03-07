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


class MdvDisplay:
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
        if ax is None:
            ax = plt.gca()
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
        if ax is None:
            ax = plt.gca()
        theta = np.linspace(0, 2 * pi, npts)
        r = np.ones([npts], dtype=np.float32) * range_ring_location_km
        x = r * np.sin(theta)
        y = r * np.cos(theta)
        ax.plot(x, y, 'k-')

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
        if ax is None:
            ax = plt.gca()

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
        if ax is None:
            ax = plt.gca()
        loc_x, loc_y = corner_to_point(self._radar_location, location)
        loc_x /= 1000.0
        loc_y /= 1000.0
        ax.plot([loc_x], [loc_y], symbol)
        ax.text(loc_x - 5.0, loc_y , label, color=text_color)
 
    def plot_cross_hair(self, size, npts=100, ax=None):
        """
        Plot a cross-hair.
        
        Parameters
        ----------
        size : float
            Size of cross-hair in km.
        npts: int
            Number of points in the cross-hair, higher for better resolution.
        ax : Axis
            Axis to plot on.  None will use the current axis.
 
        """
        if ax is None:
            ax = plt.gca()
        x = np.zeros(npts, dtype=np.float32)
        y = np.linspace(-size, size, npts)
        ax.plot(x, y, 'k-')  # verticle
        ax.plot(y, x, 'k-')  # horizontal




def single_panel_ppi(mdvfile, sweep, fld, fig = None, ax = None, mask=None, 
                     ylim=None, xlim=None, colorbar=False, vmin=None, 
                     vmax=None, range_rings=None, cross=0,
                     labels=None, label_locations=None, label_symbols='r+', 
                     label_text_color = 'k', title_str = None):
    """
    Add a PPI plot from a MdvFile object to the current axis.

    Parameters
    ----------
    mdvfile : MdvFile 
        MdvFile object to plot data from.
    sweep : int,
        Sweep number to plot.
    fld : str
        Field to plot.
    
    Other Parameters
    ----------------
    ax : axis
        Axis to add PPI plot to.  None will use the current axis. 
    mask : tuple
        2-Tuple containing the field name and value below which to mask field 
        prior to plotting, for example to mask all data where NCP < 0.5 set 
        mask to ['NCP', 0.5]. None performs no masking.
    xlim : tuple, optional
        2-Tuple containing y-axis limits in km. None uses default limits.
    ylim : tuple, optional
        2-Tuple containing x-axis limits in km. None uses default limits.
    colorbar : bool
        True to add a colorbar with label to the axis.
    vmin : float
        Luminance minimum value, None for default value.
    vmax : float
        Luminance maximum value, None for default value.
    range_rings : list of floats
        List of distances in km to plot range rings, None for no rings.
    cross : float
        Length of cross hair, 0 for no cross hair.

    labels : list of str

    label_locations : list of 2-tuples
        
    label_symbols : list of str

    label_text_color : str

    title_str : str



    Returns
    -------

    """ 
    # parse the parameters
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()

    info_dict = make_info(mdvfile, fld)
    info_dict.update({'ele': mdvfile.el_deg[sweep]})
    
    if title_str is None:
        title_str = forminator()
    
    if vmin is None:
        vmin = get_default_range(mdvfile, fld)[0]
    if vmax is None:
        vmax = get_default_range(mdvfile, fld)[1]
     
    # extract coordinates in km
    x = mdvfile.carts['x'][sweep] / 1000.0  # x coords in km
    y = mdvfile.carts['y'][sweep] / 1000.0  # y coords in km
    z = mdvfile.carts['z'][sweep] / 1000.0  # z corrds in km
   
    # read in the data, mask if needed
    data = mdvfile.read_a_field(mdvfile.fields.index(fld))[sweep]
    if mask is not None:
        mask_field, mask_value = mask
        mask_field_num = mdvfile.fields.index(mask_field)
        mdata = mdvfile.read_a_field(mask_field_num)[sweep]
        data = np.ma.masked_where(mdata < mask_value, data)

    # create the plot
    pm = ax.pcolormesh(x, y, data, vmin=vmin, vmax=vmax)
   
    # set limits
    if ylim is not None:
        ax.set_ylim(ylim)
    if xlim is not None:
        ax.set_xlim(xlim)
    
    # add a colorbar 
    if colorbar:
        add_colorbar(fig, pm, info_dict['name'])

    # add range rings
    if range_rings is not None:
        for range_ring_location_km in range_rings:
            plot_range_ring(ax, range_ring_location_km)
   
    if cross:
        plot_cross_hair(ax, cross)

    # add labels and title
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)') 
    ax.set_title(title_str % info_dict)
   
    # add labels
    radar_location = [info_dict['latitude_deg'], info_dict['longitude_deg']]
    plot_labels(radar_location, labels, label_locations, label_symbols,
               label_text_color, ax=ax)

    return pm


def plot_labels(radar_location, labels=None, label_locations=None, 
               label_symbols='r+', label_text_color='k', ax=None):
    """ Plot symbols and labels at locations """
    if ax is None:
        ax = plt.gca()

    if labels is None:
        labels = []
    if label_locations is None:
        label_locations = []
    if type(label_symbols) is str:
        label_symbols = [label_symbols] * len(labels)
    if len(labels) != len(label_locations):
        raise ValueError('length of labels and label_locations must match')
    if len(labels) != len(label_symbols):
        raise ValueError('length of labels and symbols must match')

    for loc, label, sym in zip(label_locations, labels, label_symbols):
        loc_x, loc_y = corner_to_point(radar_location, loc)
        loc_x /= 1000.0
        loc_y /= 1000.0
        ax.plot([loc_x], [loc_y], sym)
        ax.text(loc_x - 5.0, loc_y , label, color=label_text_color)
    return


def add_colorbar(fig, mappable, label):
    """ Add a colorbar to a figure. """
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
 
def plot_range_ring(ax, range_ring_location_km):
    """
    plot a circle on the current plot at range rnge
    """
    npts = 100
    theta = linspace(0, 2 * pi, npts)
    r = ones([npts], dtype=float32) * range_ring_location_km
    x = r * sin(theta)
    y = r * cos(theta)
    ax.plot(x, y, 'k-')


def plot_cross_hair(ax, rnge):
    """
    plot a cross hair of length rnge
    """
    npts = 100
    #vert
    x = zeros(npts, dtype=float32)
    y = linspace(-rnge, rnge, npts)
    ax. plot(x, y, 'k-')
    #hor
    y = zeros(npts, dtype=float32)
    x = linspace(-rnge, rnge, npts)
    ax.plot(x, y, 'k-')


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


def single_panel_rhi(mdvobj, sweep, fld, **kwargs):
    """
    single_panel_rhi(mdv_object, sweep_number, fld, rangerings=rnge_list,
                     ylim=ylim_tup, xlim=xlim_tup, mask=mask_tup )

    mdv_object : py_mdv object for the radar
    sweep_number : int of the sweep number
    fld : String of the field (eg 'DBZ_F')
    rnge_list : a list of ranges that you want rings at (omit for none)
                eg: [100.0, 200.0]
    ylim_tup : a (2) tuple to set the ylimits
    xlim_tup : a (2) tuple to set the xlimits
    mask_tup : set for masking on another moment. eg to mask all data where
               NCP < 0.5 mask_tup=['NCP', 0.5]

    """
    locs = kwargs.get('locs', [])
    labels = kwargs.get('labels', [])
    def_ranges = {
        'DBMHC': [-100, 0],
        'DBMVC': [-100, 0],
        'DBZ': [-16.0, 64.0],
        'DBZ_F': [-16.0, 64.0],
        'DBZVC': [-16.0, 64.0],
        'DBZVC_F': [-16.0, 64.0],
        'VEL': [-1.0 * mdvobj.radar_info['unambig_vel_mps'],
                mdvobj.radar_info['unambig_vel_mps']],
        'VEL_F': [-1.0 * mdvobj.radar_info['unambig_vel_mps'],
                  mdvobj.radar_info['unambig_vel_mps']],
        'WIDTH': [0.0, 10.0],
        'WIDTH_F': [0.0, 10.0],
        'ZDR': [-3, 6.0],
        'ZDR_F': [-3, 6.0],
        'RHOHV': [0.6, 1.0],
        'RHOHV_F': [0.6, 1.0],
        'PHIDP': [0, 180.0],
        'PHIDP_F': [0, 180.0],
        'KDP': [-2, 6],
        'KDP_F': [-2, 6],
        'NCP': [0, 1],
        'NCP_F': [0, 1]}

    rges = kwargs.get('rges', def_ranges[fld])
    info_dict = make_info(mdvobj, fld)
    info_dict.update({'ele': mdvobj.az_deg[sweep]})
    tit_str = kwargs.get('tit_str', forminator())
    my_title = tit_str % info_dict
    x = mdvobj.carts['x'][sweep, :, :, ]
    y = mdvobj.carts['y'][sweep, :, :, ]
    z = mdvobj.carts['z'][sweep, :, :, ]
    if 'mask' in kwargs.keys():
        mask = mdvobj.read_a_field(mdvobj.fields.index(
            kwargs['mask'][0]))[sweep, :, :]
        data = ma.masked_where(
            mask < kwargs['mask'][1],
            mdvobj.read_a_field(mdvobj.fields.index(fld))[sweep, :, :])
    else:
        data = mdvobj.read_a_field(mdvobj.fields.index(fld))[sweep, :, :]
    pcolormesh(sign(y) * sqrt(y ** 2 + x ** 2) / 1000.0, z / 1000.0,
               data, vmin=rges[0], vmax=rges[1])
    if 'rangerings' in kwargs.keys():
        for rge in kwargs['rangerings']:
            plot_rangering(rge)
    xlabel(kwargs.get('xlab', 'Range (km)'))
    ylabel(kwargs.get('ylab', 'Distance above radar (km)'))
    if 'ylim' in kwargs.keys():
        ylim(kwargs['ylim'])
    if 'xlim' in kwargs.keys():
        xlim(kwargs['xlim'])
    cb = colorbar()
    cb.set_label(info_dict['name'])
    title(my_title)
