"""
pyart.graph.common_display
==========================

Abstracr display object common

"""

import matplotlib.pyplot as plt
import numpy as np

from .common import corner_to_point


class AbstractDisplay:
    """
    AbstractDisplay is the base class for other Display classes.   
    """

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
        loc_x, loc_y = corner_to_point(self.loc, location)
        loc_x /= 1000.0
        loc_y /= 1000.0
        ax.plot([loc_x], [loc_y], symbol)
        ax.text(loc_x - 5.0, loc_y , label, color=text_color)
 
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
        if ax is None:
            ax = plt.gca()
        x = np.zeros(npts, dtype=np.float32)
        y = np.linspace(-size, size, npts)
        ax.plot(x, y, 'k-')  # verticle
        ax.plot(y, x, 'k-')  # horizontal


