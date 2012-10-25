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

def fzl_index(fzl, ranges, elevation, radar_height):
	Re=6371.0*1000.0
	p_r=4.0*Re/3.0
	z=radar_height+(ranges**2 + p_r**2 + 2.0*ranges*p_r*np.sin(elevation*np.pi/180.0))**0.5 -p_r
	return np.where(z < fzl)[0].max()

def det_process_range(radar, sweep, fzl, doc=10):
	fzl_i=fzl_index(4000.0, radar.range['data'], 
		radar.sweep_info['fixed_angle']['data'][sweep],
		radar.location['altitude']['data'])
	if fzl_i > len(radar.range['data'])-doc:
		gate_end=len(radar.range['data'])-doc
	else:
		gate_end=fzl_i
	ray_start=radar.sweep_info['sweep_start_ray_index']['data'][sweep]
	ray_end=radar.sweep_info['sweep_end_ray_index']['data'][sweep]
	return gate_end, ray_start, ray_end


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
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman', 'sg_smooth']:
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

def construct_A_matrix(n_gates, filt):
    """Construct the A A matrix which is the block matrix given by:
    $$\bf{A}=\begin{bmatrix} \bf{I} & \bf{-I} \\\\ \bf{-I}& \bf{I} \\\\ \bf{Z} & \bf{M} \end{bmatrix} $$
	where $\bf{I}$ is the identity matrix, $\bf{Z}$ is a matrix of zeros and $\bf{M}$ contains our differential constraints. Each block is of shape n_gates by n_gates making shape($\bf{A}$)=(3*n, 2*n).
	Note that $\bf{M}$ contains some side padding to deal with edge issues
	"""
    Identity=np.eye(n_gates)
    filter_length=len(filt)
    M_matrix_middle=np.diag(np.ones(n_gates-filter_length+1), k=0)*0.0
    posn=np.linspace(-1.0*(filter_length-1)/2, (filter_length-1)/2, filter_length)
    for diag in range(filter_length):
        M_matrix_middle=M_matrix_middle+np.diag(np.ones(n_gates-filter_length+1 - np.abs(posn[diag])), k=posn[diag])*filt[diag]
    side_pad=(filter_length-1)/2
    M_matrix=np.bmat([np.zeros([n_gates-filter_length+1, side_pad], dtype=float),
        M_matrix_middle, 
        np.zeros([n_gates-filter_length+1, side_pad], dtype=float)])
    Z_matrix=np.zeros([n_gates-filter_length+1, n_gates])
    return np.bmat([[Identity, -1.0*Identity], [Identity, Identity], [Z_matrix, M_matrix]])

def construct_B_vectors(phidp_mod, z_mod, filt, **kwargs):
	n_gates=phidp_mod.shape[1]
	n_rays=phidp_mod.shape[0]
	filter_length=len(filt)
	side_pad=(filter_length-1)/2
	top_of_B_vectors=np.bmat([[-phidp_mod,phidp_mod]])
	data_edges=np.bmat([phidp_mod[:, 0:side_pad], np.zeros([n_rays, n_gates-filter_length+1]), phidp_mod[:,-side_pad:]])
	ii=filter_length-1
	jj=data_edges.shape[1]-1
	list_corrl=np.zeros([n_rays, jj-ii+1])
	for count in range(list_corrl.shape[1]):
	    list_corrl[:, count]= -1.0*(np.array(filt)*(np.asarray(data_edges))[:, count:count+ii+1]).sum(axis=1)
	sct=(((10.0**(0.1*z_mod))**kwargs.get('coef',0.914))/kwargs.get('dweight',60000.0))[:, side_pad:-side_pad]
	sct[np.where(sct <0.0)]=0.0
	sct[:, 0:side_pad]=list_corrl[:, 0:side_pad]
	sct[:, -side_pad:]=list_corrl[:, -side_pad:]
	B_vectors=np.bmat([[top_of_B_vectors,sct]])
	return B_vectors


def LP_solver(A_Matrix, B_vectors, weights, it_lim=7000, presolve=True, verbose=True):
	if verbose: 
		message_state=glpk.LPX.MSG_ON
	else:
		message_state=glpk.LPX.MSG_OFF
	n_gates=weights.shape[1]/2
	n_rays=B_vectors.shape[0]
	mysoln=np.zeros([n_rays, n_gates])
	lp = glpk.LPX()        # Create empty problem instance
	lp.name = 'LP_MIN'     # Assign symbolic name to problem
	lp.obj.maximize = False # Set this as a maximization problem
	lp.rows.add(2*n_gates+n_gates-4)# Append rows 
	lp.cols.add(2*n_gates)
	glpk.env.term_on = True
	for cur_row in range(2*n_gates+n_gates-4):
		lp.rows[cur_row].matrix=list(np.squeeze(np.asarray(A_Matrix[cur_row, :])))
	for i in range(2*n_gates):
		lp.cols[i].bounds = 0.0, None
	for raynum in range(n_rays):
		this_soln=np.zeros(n_gates)
		for i in range(2*n_gates+n_gates-4):
			lp.rows[i].bounds = B_vectors[raynum, i], None
		for i in range(2*n_gates):
			lp.obj[i]=weights[raynum,i]
		lp.simplex(msg_lev=message_state, meth=glpk.LPX.PRIMAL, it_lim=it_lim, presolve=presolve)
		for i in range(n_gates):
			this_soln[i]=lp.cols[i+n_gates].primal
		mysoln[raynum, :]=smooth_and_trim(this_soln, window_len=5, window='sg_smooth')
	return mysoln

class phase_proc:
	def __init__(self,radar, offset, **kwargs):
		debug=kwargs.get('debug', False)
		if debug: print('populating')
		low_z=kwargs.get('low_z', 10.0)
		high_z=kwargs.get('high_z', 53.0)
		self.min_phidp=kwargs.get('min_phidp',0.01)
		self.min_ncp=kwargs.get('min_ncp', 0.5)
		self.min_rhv=kwargs.get('min_rhv',0.8)
		self.fzl=kwargs.get('fzl':4000.0)
		refl=copy.deepcopy(radar.fields['reflectivity_horizontal']['data'])+offset
		is_low_z=(refl) <low_z
		is_high_z=(refl) >high_z
		refl[np.where(is_high_z)]=high_z
		refl[np.where(is_low_z)]=low_z
		self.z_mod=refl
		self.not_coherent=radar.fields['norm_coherent_power']['data'] < min_ncp
		self.not_correlated=myradar.fields['copol_coeff']['data'] < min_rhv
		if debug: print('Unfolding')
		my_unf=get_phidp_unf(myradar, ncp_lev=min_ncp, rhohv_lev=min_rhv, ncpts=2, doc=None)
		my_new_ph=copy.deepcopy(radar.fields['dp_phase_shift'])
		my_unf[:,-1]=my_unf[:,-2]
		my_new_ph['data']=my_unf
		radar.fields.update({'unf_dp_phase_shift':my_new_ph})
		phidp_mod=copy.deepcopy(radar.fields['unf_dp_phase_shift']['data'])
		phidp_neg=phidp_mod < min_phidp
		phidp_mod[np.where(phidp_neg)]=min_phidp
		self.phidp_mod=phidp_mod
		self.radar=radar
	__call__(self, debug=False):
		proc_ph=copy.deepcopy(self.radar.fields['unf_dp_phase_shift'])
		proc_ph['data']=self.phidp_mod
		St_Gorlv_differential_5pts=[-.2, -.1, 0, .1, .2]
		for sweep in range(len(self.radar.sweep_info['sweep_start_ray_index']['data'])):
    		if debug:print "Doing ", sweep
    		end_gate, start_ray, end_ray=det_process_range(self.radar,sweep,self.fzl, doc=15)
    		start_gate=0
    		A_Matrix=construct_A_matrix(len(self.radar.range['data'][start_gate:end_gate]),St_Gorlv_differential_5pts )
    		B_vectors=construct_B_vectors(self.phidp_mod[start_ray:end_ray,start_gate:end_gate], self.z_mod[start_ray:end_ray,start_gate:end_gate], St_Gorlv_differential_5pts)
    		weights=np.ones(self.phidp_mod['data'][start_ray:end_ray,start_gate:end_gate].shape)
    		nw=np.bmat([weights, np.zeros(weights.shape)])
    		mysoln=phase_proc.LP_solver(A_Matrix, B_vectors, nw, it_lim=7000, presolve=True, verbose=debug)
    		proc_ph['data'][start_ray:end_ray,start_gate:end_gate]=mysoln
    		last_gates=proc_ph['data'][start_ray:end_ray,-16]
    		proc_ph['data'][start_ray:end_ray,-16:]=np.meshgrid(ones([16]), last_gates)[1]
			proc_ph['valid_min']=0.0
			proc_ph['valid_max']=400.0
			self.radar.fields.update({'proc_dp_phase_shift':proc_ph})
		return self.radar