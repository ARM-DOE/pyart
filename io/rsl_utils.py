""" Utilities for converting to RSL radar objects

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
Computing Environment and Life Sciences Directorate
ARM Climate Research Facility

United States Department of Energy

USE
---

REQUIREMENTS
------------

Ctypes, py4dd or pyrsl

HISTORY
-------
2011-07-26 Start of development 
Scott Collis scollis.acrf@gmail.com
0.1: basic functionality

"""
import ctypes
import os
import sys
pyart_dir=os.environ.get('PYART_DIR',os.environ['HOME']+'/python')
sys.path.append(pyart_dir+'/pyart/io/')
try:
	import py4dd as py4dd
except ImportError:
	print "No Py-ART in pythonpath, please either add to pythonpath or set the PYART_BASE environment variable"

def mdv_to_rsl(myfile, trans):
	radar = py4dd.RSL_new_radar(40)
	for mdv_fld in trans.keys():
		mdv_data=myfile.read_a_field(myfile.fields.index(mdv_fld))
		rsl_index=getattr(py4dd.fieldTypes(), trans[mdv_fld])
		print "Transfering ", mdv_fld, " to ", trans[mdv_fld], " which has an RSL index of ", rsl_index
		nsweeps, nrays, ngates=mdv_data.shape
		print nsweeps, nrays, ngates
		vol = py4dd.RSL_new_volume(nsweeps)
		radar.contents.v[rsl_index] = vol
		vol.contents.h.field_type = trans[mdv_fld]
		vol.contents.h.nsweeps = nsweeps
		vol.contents.h.calibr_const = 0
		field_f, field_invf = py4dd.conversion_functions[trans[mdv_fld]]
		vol.contents.h.f    = field_f
		vol.contents.h.invf = field_invf
		for i_s in range(nsweeps):
			print "Sweep: ", i_s, "of ", nsweeps
			sweep=py4dd.RSL_new_sweep(nrays)
			vol.contents.sweep[i_s] = sweep
			sweep.contents.h.field_type  = trans[mdv_fld]
			sweep.contents.h.sweep_num   = i_s  + 1 # one-indexed
			sweep.contents.h.elev        = myfile.el_deg[i_s]
			sweep.contents.h.beam_width  = 1.0 #change this
			sweep.contents.h.nrays       = nrays
			sweep.contents.h.f           = field_f
			sweep.contents.h.invf        = field_invf
			sweep_time = myfile.times['time_begin']
			if myfile.scan_type=='rhi':
				sweep_az=[myfile.az_deg[0]]*nrays
			else:
				sweep_az = myfile.az_deg
			#print sweep_az
			sweep_beam_width = 1.0
			sweep_gate_width = myfile.field_headers[1]['grid_dx']*1000
			sweep_nyquist    = myfile.radar_info['unambig_vel_mps']
			for i_r in range(nrays):
				ray=py4dd.RSL_new_ray(ngates)
				sweep.contents.ray[i_r] = ray
				ray.contents.h.nbins     = ngates
				ray.contents.h.year      = sweep_time.year
				ray.contents.h.month     = sweep_time.month
				ray.contents.h.day       = sweep_time.day
				ray.contents.h.hour      = sweep_time.hour
				ray.contents.h.minute    = sweep_time.minute
				ray.contents.h.sec       = sweep_time.second
				ray.contents.h.azimuth   = sweep_az[i_r]
				ray.contents.h.ray_num   = i_r + 1 # one-indexed
				ray.contents.h.elev      = myfile.el_deg[i_s]
				ray.contents.h.elev_num  = i_s + 1 # one-indexed
				ray.contents.h.range_bin1= int(myfile.radar_info['start_range_km']*1000.0)
				ray.contents.h.gate_size = int(myfile.field_headers[1]['grid_dx']*1000.0)
				ray.contents.h.fix_angle = myfile.el_deg[i_s]
				ray.contents.h.lat       = myfile.radar_info['latitude_deg']
				ray.contents.h.lon       = myfile.radar_info['longitude_deg']
				ray.contents.h.alt       = int(myfile.radar_info['altitude_km']*1000.0)
				ray.contents.h.beam_width= 1.0
				ray.contents.h.nyq_vel   = myfile.radar_info['unambig_vel_mps']
				ray.contents.h.f         = field_f
				ray.contents.h.invf      = field_invf
				for i_g in range(ngates):
					float_value = ctypes.c_float(mdv_data[i_s,i_r,i_g])
					range_value = ray.contents.h.invf(float_value)
					ray.contents.range[i_g] = range_value
	return radar























