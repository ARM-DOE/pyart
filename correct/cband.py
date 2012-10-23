""" Utilities working with polarimentric moments and attenuation correction of ZH and ZDR at C-band

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

Code is adapted with permission from Scott Giangrande's IDL code

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
2011-08-09 Start of development 
Scott Collis scollis.acrf@gmail.com
"""
__author__ = "Scott Collis <scollis@anl.gov>"
__version__ = "0.1"

import numpy as np
import numpy.ma as ma
import matplotlib as mpl
mpl.use('Agg')
import sys
import os
sys.path.append('/home/titan5/python/')
import pyart.io.py_mdv as py_mdv
from time import time
from itertools import groupby
import scipy.integrate
from scipy.interpolate import interp1d

def fdif(data):
	"""REALLY basic finite difference""" 
	de=np.insert(data, 0, data[0]-(data[1]-data[0]))
	dd=data-de[0:len(data)]
	return dd


#First for completeness and self containing we have the numpy cookbook example smooth algorithm 

def smooth(x,window_len=11,window='hanning'):
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
    else:
        w=eval('np.'+window+'(window_len)')
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y


def KDP(phidp, dx, window_len):
    kdp=np.zeros(phidp.shape, dtype=float)
    myshape=kdp.shape
    for swn in range(myshape[0]):
        #print "Sweep ", swn+1
        for rayn in range(myshape[1]):
            #print "ray ", rayn+1
            kdp[swn, rayn, :]=sobel(phidp[swn, rayn,:], window_len=window_len)/dx
    return kdp

def sobel(x,window_len=11):
    """Sobel differential filter for calculating KDP
    output:
        differential signal (Unscaled for gate spacing
	    example:
	"""
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    w=2.0*np.arange(window_len)/(window_len-1.0) -1.0
    #print w
    w=w/(abs(w).sum())
    y=np.convolve(w,s,mode='valid')
    return -1.0*y[window_len/2:len(x)+window_len/2]/(window_len/3.0)

def snr(line, **kwargs):
	wl=kwargs.get('wl', 11)
	signal=smooth_and_trim(line, window_len=wl)
	noise=smooth_and_trim(np.sqrt((line-signal)**2), window_len=wl)
	return abs(signal)/noise

def noise(line, **kwargs):
	wl=kwargs.get('wl', 11)
	signal=smooth_and_trim(line, window_len=wl)
	noise=np.sqrt((line-signal)**2)
	return noise

 


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



def get_phidp_unf_sg(myfile, **kwargs):
	ncp_lev=kwargs.get('ncp_lev', 0.4)
	rhohv_lev=kwargs.get('rhohv_lev', 0.6)
	debug=kwargs.get('debug', False)
	ncpts=kwargs.get('ncpts', 20)
	doc=kwargs.get('doc', -10)
	my_phidp=myfile.read_a_field(myfile.fields.index('PHIDP_F'))[:,:,0:doc]
	my_rhv=myfile.read_a_field(myfile.fields.index('RHOHV_F'))[:,:,0:doc]
	my_ncp=myfile.read_a_field(myfile.fields.index('NCP_F'))[:,:,0:doc]
	my_z=myfile.read_a_field(myfile.fields.index('DBZ_F'))[:,:,0:doc]
	t=time()
	system_zero=det_sys_phase(myfile, -141.097702627)
	cordata=np.zeros(my_rhv.shape, dtype=float)
	for sweep in range(my_rhv.shape[0]):
		if debug: print "sweep ::  ", sweep
		for radial in range(my_rhv.shape[1]):
			my_snr=snr(my_z[sweep,radial,:])
			notmeteo=np.logical_or(np.logical_or(my_ncp[sweep,radial,:] < ncp_lev, my_rhv[sweep,radial,:] < rhohv_lev), my_snr < 10.0)
			#reallynotmeteo=np.logical_or(myfile.NCP_F[sweep,radial,:] < 0.5, myfile.RHOHV_F[sweep,radial,:] < 0.95)
			x_ma=ma.masked_where(notmeteo, my_phidp[sweep,radial,:])
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
				cordata[sweep, radial, :]=unwrapped_fixed
			else:
				cordata[sweep, radial, :]=np.zeros(my_rhv.shape[2])
	if debug: print "Exec time: ", time()-t
	return cordata


def det_sys_phase_sg(myfile, fg, **kwargs):
	print "dooooing"
	ncp_lev=kwargs.get('ncp_lev', 0.4)
	rhohv_lev=kwargs.get('rhohv_lev', 0.6)
	print rhohv_lev, ncp_lev
	good=False
	n=0
	phases=[]
	mncp=myfile.NCP_F[:,:,30:]
	mrhv= myfile.RHOHV_F[:,:,30:]
	#mncp[:,:,0:30]=0.0
	#mrhv[:,:,0:30]=0.0
	for sweep in [1]:
		for radial in range(myfile.RHOHV_F.shape[1]):
			meteo=np.logical_and(mncp[sweep, radial,:] > ncp_lev, mrhv[sweep,radial,:] > rhohv_lev)
			mpts=np.where(meteo)
			#print len(mpts),  mpts[-1]
			if len(mpts[0]) > 25:
				good=True
				msmth_phidp=smooth_and_trim(myfile.PHIDP_F[sweep,radial,30:], 9)
				phases.append(msmth_phidp[mpts].min())
	print phases[0:30]
	if not(good): sys_phase=fg
	print fg
	return np.median(phases[0:30])

def det_sys_phase(myfile, fg, **kwargs):
	print "dooooing"
	ncp_lev=kwargs.get('ncp_lev', 0.4)
	rhohv_lev=kwargs.get('rhohv_lev', 0.6)
	print rhohv_lev, ncp_lev
	good=False
	n=0
	phases=[]
	mncp=myfile.NCP_F
	mrhv= myfile.RHOHV_F
	mncp[:,:,0:30]=0.0
	mrhv[:,:,0:30]=0.0
	for sweep in [1]:
		for radial in range(myfile.RHOHV_F.shape[1]):
			meteo=np.logical_and(mncp[sweep, radial,:] > ncp_lev, mrhv[sweep,radial,:] > rhohv_lev)
			mpts=np.where(meteo)
			#print len(mpts),  mpts[-1]
			if len(mpts[0]) > 25:
				good=True
				msmth_phidp=smooth_and_trim(myfile.PHIDP_F[sweep,radial,:], 20)
				phases.append(msmth_phidp[mpts].min())
	phases.sort()
	print phases[0:30]
	if not(good): sys_phase=fg
	print fg
	return np.median(phases[0:30])



def append_phidp_unf(myfile, **kwargs):
	ncp_lev=kwargs.get('ncp_lev', 0.4)
	rhohv_lev=kwargs.get('rhohv_lev', 0.6)
	debug=kwargs.get('debug', False)
	ncpts=kwargs('ncpts', 20)
	d=myfile.read_a_field(myfile.fields.index('PHIDP_F'))
	d=myfile.read_a_field(myfile.fields.index('RHOHV_F'))
	d=myfile.read_a_field(myfile.fields.index('NCP_F'))
	d=myfile.read_a_field(myfile.fields.index('DBZ_F'))
	t=time()
	system_zero=det_sys_phase(myfile, -141.097702627)
	cordata=np.zeros(myfile.RHOHV_F.shape, dtype=float)
	for sweep in range(myfile.RHOHV_F.shape[0]):
		if debug: print "sweep ::  ", sweep
		for radial in range(myfile.RHOHV_F.shape[1]):
			my_snr=snr(myfile.DBZ_F[sweep,radial,:])
			notmeteo=np.logical_or(np.logical_or(myfile.NCP_F[sweep,radial,:] < ncp_lev, myfile.RHOHV_F[sweep,radial,:] < rhohv_lev), my_snr < 10.0)
			#reallynotmeteo=np.logical_or(myfile.NCP_F[sweep,radial,:] < 0.5, myfile.RHOHV_F[sweep,radial,:] < 0.95)
			x_ma=ma.masked_where(notmeteo, myfile.PHIDP_F[sweep,radial,:])
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
				cordata[sweep, radial, :]=unwrapped_fixed
			else:
				cordata[sweep, radial, :]=np.zeros(myfile.RHOHV_F.shape[2])
	if debug: print "Exec time: ", time()-t
	old_len=len(myfile.field_headers)
	hdr=myfile.field_headers[myfile.fields.index('PHIDP_F')].copy()
	newf='PHIDP_UNF'
	hdr['field_name']=newf
	myfile.field_headers.append(hdr)
	setattr(myfile, newf, cordata)
	myfile.fields.append(newf)
	return 1


def cont_regions(logic_array):
	"""Returns contigious regions where the array is True"""
	result = []
	for k, group in groupby(enumerate(np.flatnonzero(logic_array)), lambda (i,x):i-x):
		tmp = np.array([g[1] for g in group], int)
		result.append( slice(tmp[0], tmp[-1]))
	return result

def volume_filter(volume, wlen):
	#Hanning filter on a whole volume
	funf=np.zeros(volume.shape, dtype=float)
	for sweep in range(volume.shape[0]):
		for ray in range(volume.shape[1]):
			funf[sweep, ray, :]=smooth(volume[sweep,ray,:], window_len=wlen)[wlen/2:-wlen/2+1]
	return funf

#ok now for the real deal, the attenuation correction. This code is put together with extensive help from Scott Giangrande and is based on work found in 
#"Polarimetric Attenuation Correction in Heavy Rain at C Band", JAMC, Gu etal
#The basic idea is that there is regions of heavy rain (45+dBz) you get increased phidp gradients... Now you want to actually optimize alpha and beta terms
#by doing a self consistency check with the integrated reflectivity but you want to exclude these hot cells from this calculation.


def att_corr(dbz, zdr, phidp_unf, rhohv, ncp, z, a0=0.06, b0=0.01, beta=0.8, cell_width=16, min_zdr=0.1, n_guess=30.0, fzl=4000.0, debug=False, **kwargs):
	"""Correct reflectivity and ZDR for attenuation doing a self consistency check with integral of reflectivity filtering out hot
	cell regions.
	att_corr(dbz, zdr, phidp_unf, rhohv, ncp, a0=0.06, b0=0.01, beta=0.8, cell_width=16, min_zdr=0.1, n_guess=30.0, flz=4000.0, .... )
	Input:
	dbz: Uncorrected reflectivity (dBz) shape= [sweeps,azimuths,ranges]
	zdr: differential reflectivty (dB) shape= [sweeps,azimuths,ranges]
	phidp_unf: Unfolded phidp, filtered, phidp(0)=0 (deg)
	rhohv: corel coef (ratio)
	ncp: Norm Coh Power (ratio)
	Optargs:
	a0: first guess for a term
	b0: " for b term
	beta: The exponent in the Ah = a*Zh^beta relation
	cell width: Minimal width of the heavy rain cell (gates)
	min_zdr: Zdr is not allowed to drop below this threshold
	n_guess: number of a's used to find optimal guess
	fzl: height (in m) of the top freezing level.. we do not correct from here... 
	Returns:
	Corrected_Z, Corrected_zdr shape= [sweeps,azimuths,ranges]
	"""
	if debug: print "Doing whole of volume calculations"
	israin=kwargs.get('israin',45.0)
	ismeteo=kwargs.get('ismeteo',0.8)
	iscoh=kwargs.get('iscoh',0.5)
	dbin=kwargs.get('dbin', 119.916982949)
	nsweeps, nrays, ngates=dbz.shape
	#apply first pass correction
	badpts=np.isnan(dbz)#First we need to identify NaNs in the reflectivity
	dbz[badpts]=-100.0 #Ok, this is a bit dodge, but I replace NaNs with values which will basically not contribute.. We will put the NaNs back in at the end
	c_dbz_1st=dbz+phidp_unf*a0
	#create all our boolean fields
	has_rain=c_dbz_1st >=israin #hot spots
	has_met=rhohv >=ismeteo #Correlated returns
	has_coh=ncp >= iscoh #Coherent returns for rejecting 2nd trips...
	has_liquid= z < fzl
	ac_dbz=np.copy(dbz) #this is the array that will get filled.. we initialize to the uncorrected dbz so that i regions with no rhohv > ismeteo and for areas above fzl it def. to this
	zb=10.0**(dbz*0.1*beta)
	a_guesses=-0.05+0.01*np.arange(n_guess)
	if debug: print "Done, now doing ray by ray corrections"
	#ok.. so hopefully we have completed all the full arraywise ops.. this clears overhead so now we go beam by beam
	for sweep in range(nsweeps):
		print "Doing sweep ", sweep +1, "of ", nsweeps 
		for ray in range(nrays):
			liquid=np.where(has_liquid[sweep, ray,:])[0] #first determine elements above the FZL
			if len(np.where(has_met[sweep, ray,liquid])[0]) > 5: #so at least 5 gates with correlated returns to even bother doing anything, no else for this one, the ray will just default to the 1st pass attenuation corrected dbz
				hot_cells=cont_regions(np.logical_and(dbz[sweep, ray, liquid] > 45.0, has_met[sweep, ray, liquid]))
				long_hot_cells=[]
				delta_phidp=[]
				for i in range(len(hot_cells)):
					if hot_cells[i].stop-hot_cells[i].start > cell_width: 
						long_hot_cells.append(i)
						delta_phidp.append(phidp_unf[sweep,ray,hot_cells[i].stop]-phidp_unf[sweep,ray,hot_cells[i].start])
				total_hot_cell_dphidp=np.array(delta_phidp).sum()
				if len(long_hot_cells) > 0 and total_hot_cell_dphidp > 10.0: #We really do not care about hot spots unless the cumulative dPhiDP is gt 10 degrees
					#Ok, in this area we pick a self consistent a0 using the integral of attenuation corrected reflectivity plus phidp
					#print "We have ",len(long_hot_cells), " significant hot cells "," sum phidp change: ", array(delta_phidp).sum()
					rain_cells=cont_regions(has_met[sweep, ray, liquid])
					long_rain_cells=[]
					for i in range(len(rain_cells)):
						if rain_cells[i].stop - rain_cells[i].start > cell_width: long_rain_cells.append(rain_cells[i].start)
					long_rain_cells.sort()
					rmin=np.array(long_rain_cells).min()# gate number at the start if the first cell bigger that cell_width gates
					rmax=np.where(has_met[sweep, ray, liquid])[0].max() #last correlated return under the FZL
					I_indef=scipy.integrate.cumtrapz(0.46*beta*dbin*zb[sweep,ray, rmin:rmax])#indefinate integral of the reflectivity
					I_indef=np.append(I_indef, I_indef[-1]) #integral returns 1 less element... 
					#print I_indef.shape, zb[sweep,ray, rmin:rmax].shape
					#note that I0 from Scott G's code is simply the last element of the indef integral or I_indef[-1]
					fdpmax=max(phidp_unf[sweep,ray, (rmax-6):rmax].mean(),0.0)
					#print fdpmax
					sumkdp=0.5*a0*(fdpmax  - total_hot_cell_dphidp)
					spc_atten_guesses=[]#will be filled with 30 guesses at specific attenuation
					difz=[]
					wsd=np.where(~has_rain[sweep, ray, rmin:rmax])[0] #45dbz return free gates
					for i in range(len(a_guesses)):
						qq=(10.0)**(0.1*beta*(a0*fdpmax + a_guesses[i]*total_hot_cell_dphidp))-1.0
						this_spec_at=zb[sweep,ray, rmin:rmax]*qq/(I_indef[-1] +qq*I_indef)
						spc_atten_guesses.append(this_spec_at)
						difz.append(abs(dbin*np.array(this_spec_at[wsd]).sum() - sumkdp))
					#print "Rmin: ", rmin, " Rmax: ", rmax
					my_index=difz.index(np.array(difz).min())
					my_spc_atten=spc_atten_guesses[my_index]
					abs_atten=scipy.integrate.cumtrapz(my_spc_atten)*2.0*dbin
					#print abs_atten.shape
					abs_atten_plus=np.append(abs_atten, abs_atten[-1])
					#print my_index
					#print a_guesses[my_index]
					ac_dbz[sweep, ray, rmin:rmax]= ac_dbz[sweep, ray, rmin:rmax] + abs_atten_plus# zhl(i,l) = zz(i)+2.0*bin_spacing*total(ahl(0:i,l))
				else: #Just use the basic correction 
					ac_dbz[sweep, ray, liquid]=c_dbz_1st[sweep, ray, liquid]
			#no meteo returns so leave the ray as is
		#done this ray
	#done this sweep
	ac_dbz[badpts]=np.nan #put the bad points back in
	return ac_dbz




