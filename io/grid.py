""" A general gridded data object for radar data

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
datetime
pyart (mapping)


HISTORY
-------
2012: First work started
Dec 19 2012: Started with the blank class scollis
Jan 22 2013: added functions for gridding from radar objects
"""

__author__ = "Scott Collis <scollis@anl.gov>"
__version__ = "1.0"
import os
import sys
pyart_dir=os.environ.get('PYART_DIR',os.environ['HOME']+'/python')
sys.path.append(pyart_dir)

import pyart.map.ballsy as ballsy
import numpy as np
from pyart.graph import radar_display

#now to grid
def dms_to_d(dms):
	#Convert degrees minutes seconds tuple to decimal degrees
    return dms[0]+(dms[1] + dms[2]/60.0)/60.0


def grid2(radars, **kwargs):
    """map tuples of Py-Radar objects to a cartesian grid..
    Usage: (xr,yr,zr), (nx,ny,nz), grids2=grid2((radar1, radar2, radarn,  nx=81, ny=81, nz=69, zr=(0.,17000.), yr=(-20000., 20000.), xr=(-30000., 20000.), origin=(lat, lon, height_m))
    """
    nx,ny,nz=kwargs.get('nxyz', (81,81,69))
    xr,yr,zr=kwargs.get('xyzr', ((-30000., 20000), (-20000., 20000.), (0., 17000.)))
    toa=kwargs.get('toa', 17000.0)
    cf_lat, cf_lon, cf_alt=kwargs.get('origin', (dms_to_d([36.0, 36.0, 18.35]), -1.0*dms_to_d( [97.0, 29.0, 10.69]), 320.0 /1000.0))#need to just set to the radar..
    #initialize blank arrays to be filled by the radar gate coordinates 
    x=np.array([]); y=np.array([]); z=np.array([])
    #parameters to be mapped, reflectivity should always be first!
    parms=kwargs.get('params',['corrected_reflectivity_horizontal', 'reflectivity_horizontal', 'rain_rate_A','recalculated_diff_phase'])
    data=dict([(parms[i], np.array([])) for i in range(len(parms))])
    #Loop over radars and append them to the lists to be used to call the mapping object
    for myradar in radars:
    	#calculate radar offset to the origin
        displacement=radar_display.corner_to_point([cf_lat, cf_lon], [myradar.location['latitude']['data'], myradar.location['longitude']['data']])
        #meshgrid up azimuths, ranges and elevations
        rg,azg=np.meshgrid(myradar.range['data'],myradar.azimuth['data'])
        rg,eleg=np.meshgrid(myradar.range['data'],myradar.elevation['data'])
        #Calculate cartesian locations of gates
        xdash,ydash,zdash=radar_display.radar_coords_to_cart(rg,azg, eleg) 
        del rg, azg #housekeeping
        #only use gates below toa
        within_sensible=np.where(zdash.flatten() < toa)[0]
        x_disp=displacement[0]; y_disp=displacement[1]
        #append geolocation.. 
        x=np.append( x, (xdash+x_disp).flatten()[within_sensible])
        y=np.append(y, (ydash+y_disp).flatten()[within_sensible])
        z=np.append(z, zdash.flatten()[within_sensible])
        #append gate variables.. 
        for var in parms:
            data[var]=np.append(data[var],myradar.fields[var]['data'].flatten()[within_sensible])
    #find NaNs and crazy reflectivities
    where_the_data_is_good=np.where(np.logical_and(np.isfinite(data[parms[0]]), data[parms[0]] < 100.0))[0]
    #Create a meshgrid(cube) to allow calculation of radii of influence
    zg,yg,xg= np.mgrid[zr[0]:zr[1]:nz*1j, yr[0]:yr[1]:ny*1j, xr[0]:xr[1]:nx*1j ]
    #Virtual beam width and beam spacing
    nb=1.5
    bsp=1.0
    #Query radius of influence, flattened
    qrf=(zg/20.0 + np.sqrt(yg**2 +xg**2)*np.tan(nb*bsp*np.pi/180.0) + 500.0).flatten() 
    #flattened query points
    ask=np.column_stack((xg.flatten(),yg.flatten(),zg.flatten()))
    #fill values
    badval=-9999.0
    #mask and flatten the data
    masked_data=np.ma.masked_array(np.column_stack([data[key][where_the_data_is_good] for key in parms]), np.column_stack([data[key][where_the_data_is_good] for key in parms]) == badval)
    #populate the ball tree
    mapping_obj = ballsy.BallsyMapper(np.column_stack((x[where_the_data_is_good],y[where_the_data_is_good],z[where_the_data_is_good])) ,masked_data , debug=True)
    #house keeping
    del data
    #query the tree and get the flattened interpolation
    interpol = mapping_obj(ask, qrf, debug=True, func='Barnes')
    grids={}
    #reshape and store the grids in a dictionary 
    for i in range(len(parms)):
        grids.update({parms[i]:interpol[:,i].reshape((nz,ny,nx))})
    return (xr,yr,zr), (nx,ny,nz), grids



class pyGrid:
	def __init__(self, *args, **kwargs):
		if len(args)==0:
			#initialize an empty pyGrid object
			self.fields={}
			self.metadata={}
			self.axes={}
		elif 'count' in dir(args[0]):
			#a tuple of radar objects
			#grid the data
			(xr,yr,zr), (nx,ny,nz), grids=grid2(args[0], **kwargs)
			#create the fields
			self.fields=grids
			#move metadata from the radar to the grid
			for fld in grids.keys():
				for meta in args[0][0].fields[fld].keys():
					if meta!='data':
						self.fields.update({meta:args[0][0].fields[fld][meta]})
			#create some axes
			x_array=np.linspace(xr[0],xr[1],nx)
			y_array=np.linspace(yr[0],yr[1],ny)
			z_array=np.linspace(zr[0],zr[1],nz)
			xaxis={'data':x_array, 
				'long_name':'x-coordinate in Cartesian system',
				'axis':'X', 'units':'m'}
			yaxis={'data':y_array, 
				'long_name':'y-coordinate in Cartesian system',
				'axis':'Y', 'units':'m'}
			zaxis={'data':z_array, 
				'long_name':'z-coordinate in Cartesian system',
				'axis':'Z', 'units':'m', 'positive':'up'}
			self.axes={'x_disp':xaxis, 'y_disp':yaxis, 'z_disp':zaxis}
		else:
			print("foo")
			#TBI from various grid sources, WRF etc..
