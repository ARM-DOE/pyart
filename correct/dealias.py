""" Bindings to the University of Washington 4DD code through the NASA TRMM RSL
code

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Argonne National Laboratory nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL Scott Collis BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-------------------------------------------------------------------------------

Scott Collis

Argonne National Laboratory, 
Environmental Sciences Division 
Computing Environment and Life Sciences Directorate
ARM Climate Research Facility

Brookhaven National Library,
Atmospheric Sciences Department


United States Department of Energy

USE
---


REQUIREMENTS
------------
numpy
scipy

HISTORY
-------
2012-10-25 Start of development 
Scott Collis scollis.acrf@gmail.com
"""
__author__ = "Scott Collis <scollis@anl.gov>"
__version__ = "0.1"
import os
import sys
pyart_dir=os.environ.get('PYART_DIR',os.environ['HOME']+'/python')
sys.path.append(pyart_dir)
from pyart.io import  radar
from pyart.io import py4dd, rsl_utils
import netCDF4
import datetime
from numpy import where, isnan, ma


class dealiaser:
	""" Create a first guess field from a sounding and use this to attempt to dealias and return
		dopler velocities.
		Usage radar_field_object=get_dealias_field_sounding(radar_obj, sounding_heights, sounding_wind_speeds, sounding_wind_direction)
			radar_obj: a Py-RART radar object, needs to have nyquist defined in the inst_params and have at least reflectivity_horizontal and 
			mean_doppler_velocity
			sounding_heights, sounding_wind_speeds, sounding_wind_direction: Numpy arrays in meters, degrees and m/s
			datetime_sounding: A datetime object for the mean time of the profiler
			NOTE: Due to limitations in the C code do not call with numpy arrays over 900 elements long
	"""
	def __init__(self, myradar, 
				sonding_heights, sounding_wind_speeds, 
				sounding_wind_direction, datetime_sounding):
		self.radar=myradar
		self.rsl_radar=rsl_utils.radar_to_rsl(myradar, {'reflectivity_horizontal': 'DZ', 'mean_doppler_velocity':'VR'})
		dayofyear=(datetime.datetime(datetime_sounding.year, datetime_sounding.month, datetime_sounding.day)-datetime.datetime(datetime_sounding.year, 01, 01)).days+1
		juldate=(datetime_sounding.year-int(("%(d)d" %{'d':datetime_sounding.year})[0:2])*100)*1000+dayofyear
		self.fulljuldate=juldate*10000 + datetime_sounding.hour*100 + datetime_sounding.minute
		self.sounding_heights, self.sounding_wind_speeds, self.sounding_wind_direction= sonding_heights, sounding_wind_speeds, sounding_wind_direction
	def __call__(self, prep=1, low_dbz=-10, filt=1, rsl_badval=131072, fill_value=-9999.0):
		self.my_new_volume, self.sonde_volume=py4dd.dealias_radar_array(
			self.rsl_radar,None,self.sounding_heights, self.sounding_wind_speeds,
			self.sounding_wind_direction,self.fulljuldate, prep=prep, LOWDBZ=low_dbz, filt=filt)
		dealiased_data=radar.create_cube_array_lim(self.my_new_volume[0], 
			self.my_new_volume.contents.h.nsweeps, 
			self.my_new_volume.contents.sweeps[0].h.nrays)
		dealiased_data[where(isnan(dealiased_data))]=fill_value
		dealiased_data[where(dealiased_data == rsl_badval)]=fill_value
		meta=self.radar.get_mdv_meta(self.rsl_radar, 'VEL_COR') #fetch metadata
		dealiased_fielddict={'data':ma.masked_equal(dealiased_data,-9999.0).reshape(dealiased_data.shape[0]*dealiased_data.shape[1], dealiased_data.shape[2])} 
		dealiased_fielddict.update(meta)
		return dealiased_fielddict

def find_time_in_interp_sonde(interp_sonde, radar_datetime):
	"""Take a netcdf4 object pointing to an ARM Interpsonde file and get the correct time"""
	sonde_datetimes=netCDF4.num2date(interp_sonde.variables['time'][:], interp_sonde.variables['time'].units)
	selected=sorted (sonde_datetimes, key=lambda x: abs (x-radar_datetime))[0]
	my_index=list(sonde_datetimes).index(selected)
	return interp_sonde.variables['height'][:], interp_sonde.variables['wspd'][my_index,:], interp_sonde.variables['wdir'][my_index,:]

