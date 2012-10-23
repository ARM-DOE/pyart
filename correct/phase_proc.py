""" Utilities working phase data

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
pyGLPK
HISTORY
-------
2012-09-18 Start of development 
Scott Collis scollis.acrf@gmail.com
"""
__author__ = "Scott Collis <scollis@anl.gov>"
__version__ = "0.1"

from time import time
import numpy as np
from numpy import ma
import copy
import glpk

def snr(line, **kwargs):
	wl=kwargs.get('wl', 11)
	signal=smooth_and_trim(line, window_len=wl)
	noise=smooth_and_trim(np.sqrt((line-signal)**2), window_len=wl)
	return abs(signal)/noise



def unwrap_masked(lon, centered=False, copy=True):
    """
    Unwrap a sequence of longitudes or headings in degrees.

    Optionally center it as close to zero as possible

    By default, return a copy; if *copy* is False, avoid a
    copy when possible.

    Returns a masked array only of the input is a masked array.
    """
    masked_input = ma.isMaskedArray(lon)
    if masked_input:
        fill_value = lon.fill_value
        # masked_invalid loses the original fill_value (ma bug, 2011/01/20)
    lon = np.ma.masked_invalid(lon).astype(float)
    if lon.ndim != 1:
        raise ValueError("Only 1-D sequences are supported")
    if lon.shape[0] < 2:
        return lon
    x = lon.compressed()
    if len(x) < 2:
        return lon
    w = np.zeros(x.shape[0]-1, int)
    ld = np.diff(x)
    np.putmask(w, ld > 180, -1)
    np.putmask(w, ld < -180, 1)
    x[1:] += (w.cumsum() * 360.0)
    if centered:
        x -= 360 * np.round(x.mean() / 360.0)
    if lon.mask is ma.nomask:
        lon[:] = x
    else:
        lon[~lon.mask] = x
    if masked_input:
        lon.fill_value = fill_value
        return lon
    else:
        return lon.filled(np.nan)


def smooth_and_trim(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    elif window == 'sg_smooth':
    	w=np.array([0.1, .25, .3, .25, .1])
    else:
        w=eval('np.'+window+'(window_len)')
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[window_len/2:len(x)+window_len/2]


def get_phidp_unf(radar, **kwargs):
	#['norm_coherent_power', 'reflectivity_horizontal', 'dp_phase_shift', 'doppler_spectral_width', 'diff_reflectivity', 'mean_doppler_velocity', 'copol_coeff', 'diff_phase']
	ncp_lev=kwargs.get('ncp_lev', 0.4)
	rhohv_lev=kwargs.get('rhohv_lev', 0.6)
	debug=kwargs.get('debug', False)
	ncpts=kwargs.get('ncpts', 20)
	doc=kwargs.get('doc', -10)
	#print radar.fields['dp_phase_shift']['data'][:, 0:2]
	if doc!=None:
		my_phidp=radar.fields['dp_phase_shift']['data'][:,0:doc]
		my_rhv=radar.fields['copol_coeff']['data'][:,0:doc]
		my_ncp=radar.fields['norm_coherent_power']['data'][:,0:doc]
		my_z=radar.fields['reflectivity_horizontal']['data'][:,0:doc]
	else:
		my_phidp=radar.fields['dp_phase_shift']['data']
		my_rhv=radar.fields['copol_coeff']['data']
		my_ncp=radar.fields['norm_coherent_power']['data']
		my_z=radar.fields['reflectivity_horizontal']['data']
	t=time()
	system_zero=-135.0#det_sys_phase(myfile, -141.097702627)
	cordata=np.zeros(my_rhv.shape, dtype=float)
	for radial in range(my_rhv.shape[0]):
			my_snr=snr(my_z[radial,:])
			notmeteo=np.logical_or(np.logical_or(my_ncp[radial,:] < ncp_lev, my_rhv[radial,:] < rhohv_lev), my_snr < 10.0)
			x_ma=ma.masked_where(notmeteo, my_phidp[radial,:])
			try:
				ma.notmasked_contiguous(x_ma)
				for slc in ma.notmasked_contiguous(x_ma):
					if slc.stop-slc.start < ncpts or slc.start < ncpts: #so trying to get rid of clutter and small things that should not add to phidp anyway
						x_ma.mask[slc.start-1:slc.stop+1]=True
				c=0
			except TypeError:#non sequence, no valid regions
				#print "No Valid regions"
				#sys.stderr.write(':NVR:')
				c=1 #ie do nothing
				x_ma.mask[:]=True
			except AttributeError:
				sys.stderr.write('No Valid Regions, ATTERR \n ')
				sys.stderr.write(myfile.times['time_end'].isoformat()+'\n')
				#print x_ma
				#print x_ma.mask
				c=1 #also do nothing
				x_ma.mask=True
			unwrapped=unwrap_masked(x_ma, centered=False)
			#system_zero=unwrapped[np.where(np.logical_not(reallynotmeteo))][0:30].mean()
			system_max=unwrapped[np.where(np.logical_not(notmeteo))][-10:-1].mean()-system_zero #end so no clutter expected
			unwrapped_fixed=np.zeros(len(x_ma), dtype=float)
			based=unwrapped-system_zero
			based[0]=0.0
			notmeteo[0]=False
			based[-1]=system_max
			notmeteo[-1]=False
			unwrapped_fixed[np.where(np.logical_not(based.mask))[0]]=based[np.where(np.logical_not(based.mask))[0]]
			#f=interp1d(np.where(np.logical_not(based.mask))[0], based[np.where(np.logical_not(based.mask))[0]], kind='cubic')
			#unwrapped_fixed[np.where(based.mask)[0]]=f(np.where(based.mask)[0])#, np.where(np.logical_not(based.mask))[0], based[np.where(np.logical_not(based.mask))[0]])
			if len(based[np.where(np.logical_not(based.mask))[0]]) > 11:
				unwrapped_fixed[np.where(based.mask)[0]]=np.interp(np.where(based.mask)[0], np.where(np.logical_not(based.mask))[0], smooth_and_trim(based[np.where(np.logical_not(based.mask))[0]]))
			else:
				unwrapped_fixed[np.where(based.mask)[0]]=np.interp(np.where(based.mask)[0], np.where(np.logical_not(based.mask))[0], based[np.where(np.logical_not(based.mask))[0]])
			if c!=1:
				cordata[radial, :]=unwrapped_fixed
			else:
				cordata[radial, :]=np.zeros(my_rhv.shape[1])
	if debug: print "Exec time: ", time()-t
	return cordata


# class LP_phase_proc:
# 	def __init__(self, radar, min_phidp=0.01, min_ncp=0.5, min_rhv=0.8, fzl=4000.0):
# 		self.radar=radar
# 		self.min_phidp, self.min_ncp, self.min_rhv=(min_phidp, min_rhv, min_ncp)
# 		self.is_low_z=radar.fields['reflectivity_horizontal']['data'] <10.0
# 		self.is_high_z=radar.fields['reflectivity_horizontal']['data'] >53.0
# 		z_mod=copy.deepcopy(radar.fields['reflectivity_horizontal']['data'])
# 		z_mod[np.where(self.is_high_z)]=53.0
# 		z_mod[np.where(self.is_low_z)]=10.0
# 		self.z_mod=z_mod
# 		self.not_coherent=radar.fields['norm_coherent_power']['data'] < min_ncp
# 		self.not_correlated=radar.fields['copol_coeff']['data'] < min_rhv
# 		phidp_mod=copy.deepcopy(radar.fields['unf_dp_phase_shift']['data'])
# 		self.phidp_neg=phidp_mod < min_phidp
# 		phidp_mod[np.where(self.phidp_neg)]=min_phidp
# 		self.phidp_mod=phidp_mod
# 		self.A_matrix=self.construct_A_matrix(self.radar)
# 		self.set_up_problem()
# 		phidp_proc=np.zeros(self.radar.fields['unf_dp_phase_shift']['data'].shape)
# 		is_garbage=np.logical_or(self.not_coherent , self.not_correlated)
# 		m=np.where(is_garbage)
# 		weights=np.zeros(radar.fields['reflectivity_horizontal']['data'].shape)
# 		weights[m]=1.0
# 		self.weights=weights
# 	def __call__(self):
# 		for ray_number in range(phidp_proc.shape[1]):
# 			phidp_proc[i,:]=self.process_ray()
# 		return phidp_proc
# 	def make_vectors(self):
# 		n_gates=len(radar.range['data'])
# 		filter_length=len(self.filter)
# 		side_pad=(filter_length-1)/2
# 		top_of_B_vectors=np.bmat([-self.phidp_mod,self.phidp_mod]).transpose()
# 		self_consistency_constraint=(((10.0**(0.1*self.z_mod))**0.914)/50000.0)[:, side_pad:-side_pad]
# 		self_consistency_constraint[np.where(self_consistency_constraint <0.0)]=0.0
# 		self_consistency_constraint[0:side_pad]=list_corrl[0:side_pad]
# 		data_edges=np.bmat([self.phidp_mod[0:side_pad], np.zeros(n_gates-filter_length+1),
# 			self.phidp_mod[-side_pad:]])
# 		ii=filter_length-1
# 		jj=data_edges.shape[1]-1
# 		list_corrl=np.zeros(jj-ii+1)
# 		for count in range(len(list_corrl)):
# 			list_corrl[count]= -1.0*(array(self.filter)*np.squeeze(np.asarray(data_edges))[count:count+ii+1]).sum()
# 		self_consistency_constraint[-side_pad:]=list_corrl[-side_pad:]
# 		b_array=np.bmat([np.squeeze(np.asarray(top_of_B_vector)),np.squeeze(np.asarray(self_consistency_constraint))]).transpose()
# 		n_gates=len(self.radar.range['data'])
# 		objective_weights=np.bmat([self.weights, np.zeros(n)])
# 		for i in range(2*n_gates+n_gates-4):
#     		self.problem.rows[i].bounds = b_array[i], None
#     		lp.cols.add(2*n)
# 		for i in range(2*n):
#     		lp.cols[i].bounds = 0.0, None
#     		lp.obj[i]=nw[0,i]
# 		for cur_row in range(2*n+n-4):
#     		lp.rows[cur_row].matrix=list(np.squeeze(np.asarray(A_full[cur_row, :])))
# 		lp.simplex(msg_lev=glpk.LPX.MSG_ON, meth=glpk.LPX.PRIMAL, it_lim=6000, presolve=True)
# 		mysoln=zeros(n)
# 		for i in range(n):
#     		mysoln[i]=lp.cols[i+n].primal
# 	def process_ray(self, rn):
# 		n_gates=len(radar.range['data'])
# 		filter_length=len(self.filter)
# 		side_pad=(filter_length-1)/2
# 		top_of_B_vector=np.bmat([-self.phidp_mod[rn,:],self.phidp_mod[rn,:]]).transpose()
# 		self_consistency_constraint=(((10.0**(0.1*self.z_mod[rn,:]))**0.914)/50000.0)[side_pad:-side_pad]
# 		self_consistency_constraint[np.where(self_consistency_constraint <0.0)]=0.0
# 		self_consistency_constraint[0:side_pad]=list_corrl[0:side_pad]
# 		data_edges=np.bmat([self.phidp_mod[0:side_pad], np.zeros(n_gates-filter_length+1),
# 			self.phidp_mod[-side_pad:]])
# 		ii=filter_length-1
# 		jj=data_edges.shape[1]-1
# 		list_corrl=np.zeros(jj-ii+1)
# 		for count in range(len(list_corrl)):
#     		list_corrl[count]= -1.0*(array(self.filter)*np.squeeze(np.asarray(data_edges))[count:count+ii+1]).sum()
# 		self_consistency_constraint[-side_pad:]=list_corrl[-side_pad:]
# 		b_array=np.bmat([np.squeeze(np.asarray(top_of_B_vector)),np.squeeze(np.asarray(self_consistency_constraint))]).transpose()
# 		n_gates=len(self.radar.range['data'])
# 		objective_weights=np.bmat([self.weights, np.zeros(n)])
# 		for i in range(2*n_gates+n_gates-4):
#     		self.problem.rows[i].bounds = b_array[i], None
#     		lp.cols.add(2*n)
# 		for i in range(2*n):
#     		lp.cols[i].bounds = 0.0, None
#     		lp.obj[i]=nw[0,i]
# 		for cur_row in range(2*n+n-4):
#     		lp.rows[cur_row].matrix=list(np.squeeze(np.asarray(A_full[cur_row, :])))
# 		lp.simplex(msg_lev=glpk.LPX.MSG_ON, meth=glpk.LPX.PRIMAL, it_lim=6000, presolve=True)
# 		mysoln=zeros(n)
# 		for i in range(n):
#     		mysoln[i]=lp.cols[i+n].primal
# 	def set_up_lp_problem(self):
# 		n_gates=len(self.radar.range['data'])
# 		self.problem=glpk.LPX()
# 		self.problem.name= 'LP_MIN'
# 		self.obj.maximize= False
# 		self.problem.rows.add(2*n_gates+n_gates-4)         # Append rows 
# 		glpk.env.term_on = True
# 		self.problem.cols.add(2*n_gates)
# 		for cur_row in range(2*n_gates+n_gates-4):
#     		self.problem.rows[cur_row].matrix=list(np.squeeze(np.asarray(self.A_matrix	[cur_row, :])))
# 	def construct_A_matrix(self, radar):
# 		n_gates=len(radar.range['data'])
# 		Identity=eye(n_gates)
# 		St_Gorlv_differential_5pts=[-.2, -.1, 0, .1, .2]
# 		self.filter=St_Gorlv_differential_5pts
# 		filter_length=len(self.filter)
# 		M_matrix_middle=np.diag(ones(n_gates-filter_length+1), k=0)*0.0
# 		posn=np.linspace(-1.0*(filter_length-1)/2, (filter_length-1)/2, l)
# 		for diag in range(filter_length):
#     		M_matrix_middle=M_matrix_middle+np.diag(ones(n_gates-filter_length+1 - np.abs(posn[diag])), k=posn[diag])*self.filter[diag]
# 		side_pad=(filter_length-1)/2
# 		M_matrix=np.bmat([np.zeros([n_gates-filter_length+1, side_pad], dtype=float),
# 			 M_matrix_middle, 
# 			 np.zeros([n_gates-filter_length+1, side_pad], dtype=float)])
# 		Z_matrix=zeros([n-l+1, n])
# 		retrurn np.bmat([[Identity, -1.0*Identity], [Identity, Identity], [Z_matrix, M_matrix]])
