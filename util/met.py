""" general meteorological calculations useful to other pyart modules

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

Adapted by Scott Collis, Argonne National Laboratory, 
Computing Environment and Life Sciences Directorate
ARM Climate Research Facility

United States Department of Energy

USE
---


REQUIREMENTS
------------
numpy
HISTORY
-------
2011-09-01 Start of development 
Scott Collis scollis.acrf@gmail.com
"""
import os
import numpy as np
import netCDF4
from pylab import datestr2num, date2num, num2date
import heapq

#def construct_ident_field(grid, mask, x,y, **kwargs):
	##send in a particular level, not whole grid
	##set default classifications for stataform and convective
	#badval=kwargs.get( 'badval',1.31072000e+05)
	#radius=kwargs.get('radius', 11.0*1000.0)
	#loud=kwargs.get('loud', False)
	#strat_class=kwargs.get('strat_class', 2)
	#conv_class=kwargs.get('conv_class', 1)
	#masked_class=kwargs.get('masked_class', 0)
	##set thresholding
	#conv_thresh=kwargs.get('conv_thresh', 42.43) #start with Steiner's definition
	#ZN=kwargs.get('ZN', [10.0,180.0])
	##some basic calcs, WE ASSUME EQUISPACED GRID
	#dx=x[1]-x[0]
	#dy=y[1]-y[0]
	#npx=round(radius/dx)
	#npy=round(radius/dy)
	#print npx
	#print npy
	#class_grid=zeros(grid.shape, dtype=float)
	##loop over pixels
	#for i in range(grid.shape[0]):
		#for j in range(grid.shape[1]):
			##is this a convective core?
			#if grid[i,j] == badval:
				#class_grid[i,j]=masked_class
			#elif (grid[i,j] >=conv_thresh) and ((i-npx/2 >0) and ( j-npy/2 >0) and (i+npx/2 < grid.shape[0]) and (j+npy/2 < grid.shape[1])): 
				##print 'bopop'
				## convective greater than a certain distance from the edge of the grid
				##grab a subgrid
				#subgrid=grid[(i-npx/2):(i+npx/2),( j-npy/2):( j+npy/2)]
				#submask=mask[(i-npx/2):(i+npx/2),( j-npy/2):( j+npy/2)]
				#subclass=class_grid[(i-npx/2):(i+npx/2),( j-npy/2):( j+npy/2)]
				#background_z=(subgrid*submask).sum()/submask.sum()
				#if isnan(background_z):
					#background_z=-100.0
				#if background_z < 0.0:
					#delta_z = ZN[0]
				#elif (background_z >= 0.0) and (background_z < conv_thresh):
					#delta_z= ZN[0]-(background_z**2)/ZN[1]
				#elif background_z >= conv_thresh:
					#delta_z=0.0
				#for l in range(subgrid.shape[0]):
					#for m in range(subgrid.shape[1]):
						#if subgrid[l,m]==badval:
							#subclass[l,m]=masked_class
						#elif subgrid[l,m] >= background_z+delta_z: #a convetive pixel within the convective core
							#subclass[l,m]=conv_class
						#elif subclass[l,m]!=conv_class: #dont make pixels that are already convective strataform
							#subclass[l,m]=strat_class
				#class_grid[(i-npx/2):(i+npx/2),( j-npy/2):( j+npy/2)]=subclass
			#elif (grid[i,j] >=conv_thresh) and not((i-npx/2 >0) and ( j-npy/2 >0) and (i+npx/2 < grid.shape[0]) and (j+npy/2 < grid.shape[1])):
				#class_grid[i,j]=conv_class
			#elif  (grid[i,j] <conv_thresh) and (class_grid[i,j]!=conv_class):
				#class_grid[i,j]=strat_class
	#return class_grid

def nth_smallest(n, iter):
    return heapq.nsmallest(n, iter)[-1]


def get_best_sounding(target, sdir, minl, maxl):
	sondes=os.listdir(sdir)
	sondes.sort()
	offsets=[np.abs(datestr2num(s[18:33].replace('.', ' '))-date2num(target)) for s in sondes]
	cont=True
	n=1
	while cont:
		test_sonde=sondes[offsets.index(nth_smallest(n,offsets))]
		ncf_obj=netCDF4.Dataset(sdir+test_sonde,'r')
		if ncf_obj.variables['alt'][:].min() < minl and ncf_obj.variables['alt'][:].max() > maxl:
			cont=False
			chosen_sonde=test_sonde
		ncf_obj.close()
		n=n+1
	return chosen_sonde


