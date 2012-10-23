""" Colored ray plotter

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

United States Department of Energy

USE
---


REQUIREMENTS
------------
numpy
matplotlib
-------
2012-09-18 Start of development 
Scott Collis scollis.acrf@gmail.com
"""
__author__ = "Scott Collis <scollis@anl.gov>"
__version__ = "0.1"

from pylab import figure, plot, title, xlabel, ylabel, subplot, xlim, ylim
import numpy as np

class rayplot:
	def __init__(self, radar, ncp_min=0.5, rhv_min=0.8, **kwargs):
		self.radar=radar
		self.ncp_min=ncp_min
		self.rhv_min=rhv_min
		self.is_coh=radar.fields['norm_coherent_power']['data'] > self.ncp_min
		self.is_cor=radar.fields['copol_coeff']['data'] > self.ncp_min
		self.is_good=np.logical_and(self.is_coh, self.is_cor)
		self.good_idx=np.where(self.is_good)
		self.bad_idx=np.where(np.logical_not(self.is_good))
		self.distance=radar.range['data']
		self.figsize=kwargs.get('figsize', [15,10])
	def __call__(self, ray, vars, limits=[0,-1] ):
		nfig=len(vars)
		self.fig=figure(figsize=self.figsize)
		ax_list=[]
		plot_list={}
		i=1
		ray_number=self.radar.sweep_info['sweep_start_ray_index']['data'][ray[0]]+ray[1]
		for var in vars:
			cur_ax=subplot(nfig, 1, i)
			cur_good=plot(self.distance[limits[0]:limits[1]], np.ma.masked_where(np.logical_not(self.is_good), self.radar.fields[var]['data'])[ray_number, limits[0]:limits[1]])
			cur_bad=plot(self.distance[limits[0]:limits[1]], np.ma.masked_where(self.is_good, self.radar.fields[var]['data'])[ray_number, limits[0]:limits[1]], 'r-')
			ylabel(self.radar.fields[var]['units'])
			title(self.radar.fields[var]['standard_name'].replace('_',' '))
			cur_ax.yaxis.grid(color='gray', linestyle='dashed')
			if i!= nfig: 
				cur_ax.get_xaxis().set_ticks([])
			else:
				xlabel(self.radar.range['standard_name']+' ('+self.radar.range['units']+')')
			if var =='copol_coeff':
				ylim([.5,1.0])
			plot_list.update({var+'_good':cur_good})
			plot_list.update({var+'_bad':cur_bad})
			ax_list.append(cur_ax)
			i=i+1


