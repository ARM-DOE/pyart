""" Calculating attenuation from polarimetric radars

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

Code adapted from PAPER by Scott Giangrande et al

Adapted by Scott Collis and Scott Giangrande, 

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
0.5 Scott Collis: Added filtering to the reflectivity to remove second trips etc...
"""
__author__ = "Scott Collis <scollis@anl.gov>"
__version__ = "0.5"
import os
import sys
pyart_dir=os.environ.get('PYART_DIR',os.environ['HOME']+'/python')
sys.path.append(pyart_dir)
from pyart.correct import phase_proc
import numpy as np
import scipy
import copy


def calculate_attenuation(radar,z_offset, **kwargs):
	doc=kwargs.get('doc', 15)
	fzl=kwargs.get('fzl', 4000.0)
	rhv_min=kwargs.get('rhv_min',0.8)
	ncp_min=kwargs.get('ncp_min',0.5)
	a_coef=kwargs.get('a_coef',0.06)
	beta=kwargs.get('beta',0.8)
	debug=kwargs.get('debug',False)
	is_cor=radar.fields['copol_coeff']['data'] > rhv_min
	is_coh= radar.fields['norm_coherent_power']['data'] > ncp_min
	is_good = np.logical_and(is_cor, is_cor)
	good_indeces=np.where(is_good)
	refl=ma.masked_where(copy.deepcopy(radar.fields['reflectivity_horizontal']['data'])+z_offset, np.logical_not(is_good))
	ref_init_correct=refl +radar.fields['proc_dp_phase_shift']['data']*a_coef
	npts_good=len(good_indeces[0])
	dr=(radar.range['data'][1]-radar.range['data'][0])/1000.0
	specific_atten=np.zeros(radar.fields['reflectivity_horizontal']['data'].shape)
	atten=np.zeros(radar.fields['reflectivity_horizontal']['data'].shape)
	for sweep in range(len(radar.sweep_info['sweep_start_ray_index']['data'])):
		if debug: print "Doing ", sweep
		end_gate, start_ray, end_ray=phase_proc.det_process_range(radar,sweep,fzl, doc=doc)
		for i in range(start_ray,end_ray+1):
			is_good= np.logical_and(is_cor[i,0:end_gate], is_cor[i,0:end_gate])
			good_indeces=np.where(is_good)
			maximum_phidp=np.median(radar.fields['proc_dp_phase_shift']['data'][i,0:end_gate][good_indeces[0][-6:]])
			smoothed_reflectivity=phase_proc.smooth_and_trim(ref_init_correct[i,0:end_gate], window_len=5)
			reflectivity_in_linear_units=10.0**(0.1*beta*smoothed_reflectivity)
			I_indef=scipy.integrate.cumtrapz(0.46*beta*dr*reflectivity_in_linear_units[::-1])
			I_indef=np.append(I_indef, I_indef[-1])[::-1]
			self_cons_number=10.0**(0.1*beta*a_coef*maximum_phidp)-1.0
			#qq = 10.0^(0.1*beta*al*fmax(j))-1  ; does not require an 'a' coefficient because Z already corrected
			specific_atten[i,0:end_gate]=reflectivity_in_linear_units*self_cons_number/(I_indef[0] + self_cons_number*(I_indef))
			atten[i,:-1]=scipy.integrate.cumtrapz(specific_atten[i,:])*dr*2.0
			atten[i,-1]=atten[i,-2]
	spec_at=copy.deepcopy(radar.fields['diff_phase'])
	spec_at['data']=specific_atten
	spec_at['valid_min']=0.0
	spec_at['valid_max']=1.0
	spec_at['standard_name']='specific_attenuation'
	spec_at['long_name']='specific_attenuation'
	spec_at['least_significant_digit']=4
	spec_at['units']='dB/km'
	cor_z=copy.deepcopy(radar.fields['reflectivity_horizontal'])
	cor_z['data']=ref_init_correct#atten+cor_z['data']+z_offset
	cor_z['least_significant_digit']=2
	cor_z['valid_max']=80.0
	return spec_at, cor_z

