"""
pyart.graph.radar_display
=========================

General class for creating plots from Radar objects.

.. autosummary::
    :toctree:: generated/

"""

from pylab import gca, pcolormesh, colorbar, meshgrid, plot, get_cmap, text
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import num2date

from .common import corner_to_point, ax_radius, dt_to_dict 
from .common import radar_coords_to_cart as _radar_coords_to_cart

from .common_display import AbstractDisplay

def radar_coords_to_cart(rng, az, ele, debug=False):
    """
    rng in meters
    """
    return _radar_coords_to_cart(rng / 1000.0, az, ele, debug=False)


class RadarDisplay(AbstractDisplay):
    """
    Class for display radar data stored in a python radar object

    Example code:
    myradar = radar.Radar(rsl_object,
        add_meta={'instrument_name': testfilename.split('/')[-2]})
    f = figure()
    my_display = radar_display(myradar)
    my_display.plot_ppi('reflectivity_horizontal', 1)
    my_display.add_cb()
    savefig('test.dir/' +
            my_display.generate_filename('reflectivity_horizontal', 1))
    close(f)

    """

    def __init__(self, radar, **kwargs):
        # initialize the object with all the information needed to make
        # neat code downstream
        self.fields = radar.fields  # data
        self.fixed_angle = radar.sweep_info['fixed_angle']['data']
        self.scan_type = radar.scan_type
        self.ranges = radar.range['data']
        self.azimuths = radar.azimuth['data']
        self.elevations = radar.elevation['data']

        if 'shift' in kwargs:
            self.origin = 'Origin'
        else:
            self.origin = 'Radar'
        self.shift = kwargs.get('shift', [0.0, 0.0])
        rg, azg = meshgrid(self.ranges, self.azimuths)
        rg, eleg = meshgrid(self.ranges, self.elevations)

        # appending carts
        self.x, self.y, self.z = radar_coords_to_cart(rg, azg, eleg)
        self.x = self.x + self.shift[0]
        self.y = self.y + self.shift[1]

        #append a datetime object
        self.time_begin = num2date(radar.time['data'][0], radar.time['units'],
                                   calendar=radar.time['calendar'])
        self.radar_name = radar.metadata['instrument_name']
        self.plots = []
        self.plot_vars = []
        self.cbs = []
        self.starts = radar.sweep_info['sweep_start_ray_index']['data']
        self.ends = radar.sweep_info['sweep_end_ray_index']['data']
        self.loc = [radar.location['latitude']['data'],
                    radar.location['longitude']['data']]

    def generate_filename(self, var, tilt):
        infodict = {'name': self.radar_name, 'tilt': tilt, 'var': var}
        infodict.update(dt_to_dict(self.time_begin, pref='begin_'))
        return '%(name)s_%(var)s_%(tilt)02d_%(begin_year)04d%(begin_month)02d%(begin_day)02d%(begin_hour)02d%(begin_minute)02d.png' % infodict

    def generate_title(self, var, tilt):
        infodict = {'name': self.radar_name, 'tilt': self.fixed_angle[tilt],
                    'var': var.replace('_', ' ')}
        infodict.update(dt_to_dict(self.time_begin, pref='begin_'))
        return '%(name)s %(var)s %(tilt).1f deg %(begin_year)04d%(begin_month)02d%(begin_day)02d%(begin_hour)02d%(begin_minute)02d' % infodict

    def append_x(self, ax=None):
        if ax is None:
            ax = plt.gca()
        ax.set_xlabel('East West distance from ' + self.origin + ' (km)')

    def append_y(self, ax=None):
        if ax is None:
            ax = plt.gca()
        ax.set_ylabel('North South distance from ' + self.origin + ' (km)')

    def append_r(self, ax=None):
        if ax is None:
            ax = plt.gca()
        ax.set_xlabel('Distance from ' + self.origin + ' (km)')

    def append_z(self, ax=None):
        if ax is None:
            ax = plt.gca()
        ax.set_ylabel('Distance Above ' + self.origin + '  (km)')

    def plot_ppi(self, field, tilt, mask=None, vmin=None, vmax=None, 
            cmap='jet', title_and_labels=False, colorbar=False, 
            ax=None, fig=None):
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
        mask : tuple
            2-Tuple containing the field name and value below which to mask
            field prior to plotting, for example to mask all data where 
            NCP < 0.5 set mask to ['NCP', 0.5]. None performs no masking.
        vmin : float
            Luminance minimum value, None for default value.
        vmax : float
            Luminance maximum value, None for default value.
        cmap : str
            Matplotlib colormap name.   
        title_and_labels : bool
            True to set default figure titles and axis labels.
        colorbar : bool
            True to add a colorbar with label to the axis.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.
        
        """
        if ax is None:
            ax = plt.gca()
        if fig is None:
            fig = plt.gcf()
        start_index = self.starts[tilt]
        end_index = self.ends[tilt]
      
        field_dict = self.fields[field]

        if vmin is None:
            if 'valid_min' in field_dict:
                vmin = field_dict['valid_min']
            else:
                vmin = -6
        if vmax is None:
            if 'valid_max' in field_dict:
                vmax = field_dict['valid_max']
            else:
                vmax = 100

        x = self.x[start_index:end_index] / 1000.0
        y = self.y[start_index:end_index] / 1000.0
        data = field_dict['data'][start_index:end_index, :]

        if mask is not None:
            mask_field, mask_value = mask
            mdata = self.fields[mask_field]['data'][start_index:end_index]
            data = np.ma.masked_where(mdata < mask_value, data)

        pm = ax.pcolormesh(x, y, data, vmin=vmin, vmax=vmax, cmap=cmap) 
        
        if title_and_labels:
            self.append_x(ax)
            self.append_y(ax)
            ax.set_title(self.generate_title(field, tilt))

        self.plots.append(pm)
        self.plot_vars.append(field)
        
        if colorbar:
            self.plot_colorbar(mappable=pm, fig=fig) 

    def plot_rhi(self, field, tilt, mask=None, vmin=None, vmax=None, 
            cmap='jet', title_and_labels=False, colorbar=False, 
            ax=None, fig=None):
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
        mask : tuple
            2-Tuple containing the field name and value below which to mask
            field prior to plotting, for example to mask all data where 
            NCP < 0.5 set mask to ['NCP', 0.5]. None performs no masking.
        vmin : float
            Luminance minimum value, None for default value.
        vmax : float
            Luminance maximum value, None for default value.
        cmap : str
            Matplotlib colormap name.   
        title_and_labels : bool
            True to set default figure titles and axis labels.
        colorbar : bool
            True to add a colorbar with label to the axis.
        ax : Axis
            Axis to plot on. None will use the current axis.
        fig : Figure
            Figure to add the colorbar to. None will use the current figure.
        
        """ 
        if ax is None:
            ax = plt.gca()
        if fig is None:
            fig = plt.gcf()
        start_index = self.starts[tilt]
        end_index = self.ends[tilt]
      
        field_dict = self.fields[field]

        if vmin is None:
            if 'valid_min' in field_dict:
                vmin = field_dict['valid_min']
            else:
                vmin = -6
        if vmax is None:
            if 'valid_max' in field_dict:
                vmax = field_dict['valid_max']
            else:
                vmax = 100

        x = self.x[start_index:end_index] / 1000.0
        y = self.y[start_index:end_index] / 1000.0
        z = self.z[start_index:end_index] / 1000.0
        data = field_dict['data'][start_index:end_index, :]

        if mask is not None:
            mask_field, mask_value = mask
            mdata = self.fields[mask_field]['data'][start_index:end_index]
            data = np.ma.masked_where(mdata < mask_value, data)

        R = np.sqrt(x ** 2 + y ** 2) * np.sign(y) 
        pm = ax.pcolormesh(R, z, data, vmin=vmin, vmax=vmax, cmap=cmap)
        
        if title_and_labels:
            self.append_r(ax)
            self.append_z(ax)
            ax.set_title(self.generate_title(field, tilt))
        
        self.plots.append(pm)
        self.plot_vars.append(field)

        if colorbar:
            self.plot_colorbar(mappable=pm, fig=fig)

    def labelator(self, standard_name, units):
        """ Return a label (str) for the given field name and units. """
        return standard_name.replace('_', ' ') + ' (' + units + ')'

    def plot_colorbar(self, mappable=None, label=None, cax=None, fig=None):
        """
        Plot a colorbar

        Parameters
        ----------
        mappable : Image, ContourSet, etc.
            Image, ContourSet, etc to which the colorbar applied.  If None the
            last mappable object will be used.
        label :
            Colorbar label.  None will use a default value from the last field
            plotted.
        cax : Axis
            Axis onto which the colorbar will be drawn.  None is also valid.
        fig : Figure
            Figure to place colorbar on.  None will use the current figure.
  
        """
        if fig is None:
            fig = plt.gcf()
        if mappable is None:
            mappable = self.plots[-1]
        if label is None:
            last_field_dict = self.fields[self.plot_vars[-1]]
            if 'standard_name' in last_field_dict:
                standard_name = last_field_dict['standard_name']
            else:
                standard_name = last_field_dict['long_name']
            units = last_field_dict['units']
            label = self.labelator(standard_name, units)

        cb = fig.colorbar(mappable, cax=cax)
        cb.set_label(label)
        self.cbs.append(cb)
