""" A class to make nice plots from the radar class

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL Scott Collis or Argonne National Laboratory BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-------------------------------------------------------------------------------


USE
---
for example: 

REQUIREMENTS
------------
numpy
netCDF4
radar
matplotlib

HISTORY
-------
Oct 19 2012: Started development
"""


#import os
#import sys
#pyart_dir=os.environ.get('PYART_DIR',os.environ['HOME']+'/python')
#sys.path.append(pyart_dir)
#from pyart.io import radar, py4dd
from pylab import gca, pcolormesh, colorbar, meshgrid, plot
import numpy as np
from netCDF4 import num2date

def corner_to_point(corner, point):
	print corner, point
        Re=6371.0*1000.0
        Rc=ax_radius(point[0], units='degrees')
        #print Rc/Re
        y=((point[0]-corner[0])/360.0)*np.pi*2.0*Re
        x=((point[1]-corner[1])/360.0)*np.pi*2.0*Rc
        return x,y

def ax_radius(lat, units='radians'):
        #Determine the radius of a circle of constant longitude at a certain
        #Latitude
	Re=6371.0*1000.0
        if units=='degrees':
                const=np.pi/180.0
        else:
                const=1.0
        R=Re*np.sin(np.pi/2.0 - np.abs(lat*const))
        return R



def plot_x(rnge):
        npts=100
        #vert
        x=zeros(npts, dtype=float32)
        y=linspace(-rnge, rnge, npts)
        plot(x,y,'k-')
        #hor
        y=zeros(npts, dtype=float32)
        x=linspace(-rnge, rnge, npts)
        plot(x,y,'k-')



def dt_to_dict(dt, **kwargs):
        pref=kwargs.get('pref', '')
        return dict( [(pref+key, getattr(dt, key)) for key in    ['year', 'month', 'day', 'hour', 'minute', 'second']])

def radar_coords_to_cart(rng, az, ele, debug=False):
    """
    Asumes standard atmosphere, ie R=4Re/3
    Note that this v
    """
    Re=6371.0*1000.0
    #h=(r^2 + (4Re/3)^2 + 2r(4Re/3)sin(ele))^1/2 -4Re/3
    #s=4Re/3arcsin(rcos(ele)/(4Re/3+h))
    p_r=4.0*Re/3.0
    rm=rng
    z=(rm**2 + p_r**2 + 2.0*rm*p_r*np.sin(ele*np.pi/180.0))**0.5 -p_r
    #arc length
    s=p_r*np.arcsin(rm*np.cos(ele*np.pi/180.)/(p_r+z))
    if debug: print "Z=", z, "s=", s
    y=s*np.cos(az*np.pi/180.0)
    x=s*np.sin(az*np.pi/180.0)
    return x,y,z


class radar_display:
	"""Class for display radar data stored in a python radar object
	Example code:
	myradar= radar.Radar(rsl_object, add_meta={'instrument_name': testfilename.split('/')[-2]})
	f=figure()
	my_display=radar_display(myradar)
	my_display.plot_ppi('reflectivity_horizontal', 1)
	my_display.add_cb()
	savefig('/home/sc8/results/display_tests/'+my_display.generate_filename('reflectivity_horizontal', 1))
	close(f)
	"""
	def __init__(self, radar, **kwargs):
		#initialize the object with all the information needed to make neat code downstream
		self.fields=radar.fields #data
		self.scan_type=radar.scan_type
		self.ranges=radar.range['data']
		self.azimuths=radar.azimuth['data']
		self.elevations=radar.elevation['data']
		rg,azg=meshgrid(self.ranges,self.azimuths)
		rg,eleg=meshgrid(self.ranges,self.elevations)
		self.x,self.y,self.z=radar_coords_to_cart(rg,azg, eleg) #appending carts
		self.time_begin=num2date(radar.time['data'][0], radar.time['units'], calendar=radar.time['calendar']) #append a datetime object
		self.radar_name=radar.metadata['instrument_name']
		self.plots=[]
		self.plot_vars=[]
		self.cbs=[]
		self.starts, self.ends=radar.sweep_info['sweep_start_ray_index']['data'], radar.sweep_info['sweep_end_ray_index']['data']
		self.loc=[radar.location['latitude']['data'], radar.location['longitude']['data']]
	def generate_filename(self, var, tilt):
		infodict={'name':self.radar_name,'tilt':tilt, 'var':var}
		infodict.update(dt_to_dict(self.time_begin, pref='begin_'))
		return '%(name)s_%(var)s_%(tilt)02d_%(begin_year)04d%(begin_month)02d%(begin_day)02d%(begin_hour)02d%(begin_minute)02d.png' %infodict
	def generate_title(self, var, tilt):
		infodict={'name':self.radar_name,'tilt':self.elevations[tilt], 'var':var.replace('_', ' ')}
		infodict.update(dt_to_dict(self.time_begin, pref='begin_'))
		return '%(name)s %(var)s %(tilt).1f deg %(begin_year)04d%(begin_month)02d%(begin_day)02d%(begin_hour)02d%(begin_minute)02d' %infodict
	def append_x(self, **kwargs):
		kwargs.get('axis', gca()).set_xlabel('East West distance from radar (km)')
	def append_y(self, **kwargs):
		kwargs.get('axis', gca()).set_ylabel('North South distance from radar (km)')
	def plot_ppi(self, var, tilt, **kwargs):
		start_index=self.starts[tilt]
		end_index=self.ends[tilt]
		this_plot=pcolormesh(self.x[start_index:end_index, :]/1000.0, self.y[start_index:end_index, :]/1000.0,
			self.fields[var]['data'][start_index:end_index, :], #note we assume a masked array here.. if you want you can always mask the data field
			vmin=kwargs.get('vmin', self.fields[var]['valid_min']), vmax=kwargs.get('vmin', self.fields[var]['valid_max']))
		self.plots.append(this_plot)
		self.plot_vars.append(var)
	def labelator(self, standard_name, units):
		return standard_name.replace('_', ' ')+ ' (' + units +')'
	def add_cb(self, **kwargs):
		"""
		Adds a colorbar to a plot,
		defaults to the last created.. 
		"""
		if 'target_axis' in kwargs.keys():
			this_cb=colorbar(cax=kwargs['target_axis'], mappable=kwargs.get('plot', self.plots[-1]))
		else:
			this_cb=colorbar(mappable=kwargs.get('plot', self.plots[-1]))
		this_cb.set_label(kwargs.get('label', self.labelator(self.fields[self.plot_vars[-1]]['standard_name'], self.fields[self.plot_vars[-1]]['units'])))
		self.cbs.append(this_cb)
	def plot_rangering(self, rnge, **kwargs):
		"""
		Plots range rings on your PPI
		"""
		npts=kwargs.get('npts', 100)
		use_axis=kwargs.get('axis', gca())
		theta=np.linspace(0,2*np.pi, npts)
		r=np.ones([npts])*rnge
		x=r*np.sin(theta)
		y=r*np.cos(theta)
		use_axis.plot(x,y,'k-')
	def plot_x(self, rnge, **kwargs):
		"""
		Puts cross hairs on the ppi plot
		"""
		npts=kwargs.get('npts', 100)
		use_axis=kwargs.get('axis', gca())
		#vert
		x=np.zeros(npts)
		y=np.linspace(-rnge, rnge, npts)
		use_axis.plot(x,y,'k-')
		#hor
		y=np.zeros(npts)
		x=np.linspace(-rnge, rnge, npts)
		use_axis.plot(x,y,'k-')
	def plot_locs(self, locs, labels, **kwargs):
		use_axis=kwargs.get('axis', gca())
		for i in range(len(locs)):
			carts=corner_to_point(self.loc, locs[i])
			use_axis.plot([carts[0]/1000.0, carts[0]/1000.0], [carts[1]/1000.0, carts[1]/1000.0], kwargs.get('sym', ['r+']*len(locs))[i])
			use_axis.text(carts[0]/1000.0-5.0, carts[1]/1000.0, labels[i], color=kwargs.get('tc', 'k'))









