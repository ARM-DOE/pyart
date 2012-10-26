#!/usr/bin/python
	
""" Example command line utility for plotting two panels from a radar file

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
 Scott Collis, Argonne National Laboratory, 
Environmental Sciences Division
ARM Climate Research Facility

United States Department of Energy

USE
---


REQUIREMENTS
------------
numpy
HISTORY
-------
2012-10-26 Start of development 
Scott Collis scollis.acrf@gmail.com
"""
import matplotlib
matplotlib.use('agg')
import netCDF4
import sys
import os
pyart_dir=os.environ.get('PYART_DIR',os.environ['HOME']+'/python')
sys.path.append(pyart_dir)
from pyart.graph import rayplot, radar_display
from pyart.io import radar, py4dd, py_mdv
from pylab import *
import copy
if __name__ == "__main__":
	filename=sys.argv[1]
	outdir=sys.argv[2]
	tilt=int(sys.argv[3])
	print "plotting ", filename, " to ", outdir
	if ".mdv" in filename:
		my_object=py_mdv.read_mdv(filename, debug=True)
		myradar=radar.Radar(my_object)
	elif ".nc" in filename:
		my_object=netCDF4.Dataset(filename)
		myradar=radar.Radar(my_object)
		my_object.close()
	else:
		py4dd.RSL_radar_verbose_on()
		my_object = py4dd.RSL_anyformat_to_radar(filename)
	my_display=radar_display.radar_display(myradar)
	myradar.instrument_name="Radar"
	f=figure(figsize=[12,5])
	subplot(1,2,1)
	my_display.plot_ppi('mean_doppler_velocity',1)
	my_display.add_cb()
	pl.subplot(1,2,2)
	my_display.plot_ppi('corrected_	mean_doppler_velocity',1)
	my_display.add_cb()	
	figname=my_display.generate_filename('six_panel', tilt)
	savefig(outdir+'/'+figname.replace(' ','_'))
	close(f)
