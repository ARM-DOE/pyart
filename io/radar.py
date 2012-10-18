""" A general central radial scanning (or dwelling) instrument class

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
xsapr = py4dd.RSL_anyformat_to_radar(filename)
myradar=radar.Radar(xsapr)

works with rsl and MDV objects at the moment

REQUIREMENTS
------------
numpy
netCDF4
py4DD
datetime

HISTORY
-------
2012: First work started
Oct 18 2012:Updated sys.path,append to work on multi-system installs

"""



import sys
import os
from numpy import tile,array, isnan, where, ma, linspace, arange, zeros, float32
from netCDF4 import date2num
pyart_dir=os.environ.get('PYART_DIR',os.environ['HOME']+'/python')
sys.path.append(pyart_dir+'/pyart/io/')
import py4dd
from datetime import datetime


def dms_to_d(dms):
    return dms[0]+(dms[1] + dms[2]/60.0)/60.0

def csapr_standard_names():
	prop_names={'DBZ_F':'reflectivity_horizontal', 'VEL_F':'mean_doppler_velocity', 
	'WIDTH_F':'doppler_spectral_width','ZDR_F':'diff_reflectivity', 'RHOHV_F':'copol_coeff',
	'NCP_F':'norm_coherent_power', 'KDP_F':'diff_phase', 'PHIDP_F':'dp_phase_shift', 'VEL_COR':'corrected_mean_doppler_velocity', 
	'PHIDP_UNF':'unfolded_dp_phase_shift', 'KDP_SOB':'recalculated_diff_phase',
	'DBZ_AC':'attenuation_corrected_reflectivity_horizontal' }
	return prop_names

def dt_to_dict(dt, **kwargs):
	pref=kwargs.get('pref', '')
	return dict( [(pref+key, getattr(dt, key)) for key in    ['year', 'month', 'day', 'hour', 'minute', 'second']])

def defaut_mdv_metadata_map():
	"""produce the default mappings from mdv name space to cf-radial name space"""
	mdm={'instrument_name':'data_set_source', 'source':'data_set_info'}
	return mdm

def get_avail_moments(volumes):
	av=[]
	for i in range(len(volumes)):
		if volumes[i]!=None:
			av.append(py4dd.fieldTypes().list[i])
	return av

def create_cube_array(volume):
	ppi=zeros([volume.h.nsweeps, volume.sweeps[0].h.nrays, volume.sweeps[0].rays[0].h.nbins], dtype=float32)+1.31072000e+05
	for levnum in range(volume.h.nsweeps):
		for raynum in range(volume.sweeps[0].h.nrays):
			data=volume.sweeps[levnum].rays[raynum].data
			ppi[levnum,raynum, 0:len(data)]=data
	return ppi

def create_cube_array_lim(volume, nsweeps, nrays):
	ppi=zeros([nsweeps, nrays, volume.sweeps[0].rays[0].h.nbins], dtype=float32)+1.31072000e+05
	for levnum in range(nsweeps):
		for raynum in range(nrays):
			data=volume.sweeps[levnum].rays[raynum].data
			ppi[levnum,raynum, 0:len(data)]=data
	return ppi


class Radar:
	"""
	A class for storing antenna coordinate radar data which will interact nicely with CF-Radial files and other pyart code
	"""
	def __init__(self, radarobj):
		#check format of file object
		if 'vlevel_headers' in dir(radarobj):#this is a mdv file
			self.mdv2rad(radarobj)
		elif 'contents' in dir(radarobj):#this is a ctype, possibly a rsl object
			if 'h' in dir(radarobj.contents): #yep a rsl object
				self.rsl2rad(radarobj)
	def ray_header_time_to_dict(self, h):
		return {'year':h.year, 'month': h.month, 'day':h.day,'hour':h.hour, 'minute':h.minute, 'second':h.sec}
	def extract_rsl_pointing(self, volume):
		#this needs to be moved into C
		print self.nsweeps, self.nrays
		azimuth=zeros([self.nsweeps, self.nrays], dtype=float)
		elevation=zeros([self.nsweeps, self.nrays], dtype=float)
		for i in range(self.nsweeps):
			print i
			print len(volume.sweeps[i].rays)
			azimuth[i,:]=[volume.sweeps[i].rays[k].h.azimuth for k in range(self.nrays)]
			elevation[i,:]=[volume.sweeps[i].rays[k].h.elev for k in range(self.nrays)]
		return azimuth, elevation
	#time=zeros([self.nsweeps, self.nrays], dtype=float)
		#self.tu="seconds since %(year)d-%(month)02d-%(day)02d %(hour)02d:%(minute)02d:%(second)02d.0" % self.ray_header_time_to_dict(volume.sweeps[0].rays[0].h)
		#self.cal="gregorian"
			#time[i,:]=[date2num(datetime(volume.sweeps[i].rays[k].h.year, 
				#volume.sweeps[i].rays[k].h.month, 
				#volume.sweeps[i].rays[k].h.day, 
				#volume.sweeps[i].rays[k].h.hour, 
				#volume.sweeps[i].rays[k].h.minute,
				#volume.sweeps[i].rays[k].h.sec),self.tu, self.cal) for k in range(self.nrays)]
	def prtmode(self,h):
		if h.prf2 !=h.prf:
			mode='dual'
		else:
			mode='fixed'
		return mode
	def rsl2rad(self, radarobj, **kwargs):
		#We only want to transfer fields that we have valid names for... 
		#An issue that needs to be resolved is that this code likes all sweeps to have the same number of rays.. so for now we take min(nrays) across sweeps and drop rays out side of this... this is an "easy" issue to resolve caused by the fact I have been treating things as cubes and then flattening them
		name_transfer={'ZT':'DBZ', 'VR':'VEL_F', 'DR':'ZDR', 'KD':'KDP_F', 'SQ':'NCP_F', 'PH':'PHIDP_F', 'VE':'VEL_COR','RH':'RHOHV_F', 'DZ':'DBZ_F', 'SW':'WIDTH', 'ZD':'ZDR_F'}
		available_data=get_avail_moments(radarobj.contents.volumes)
		print available_data
		fields=[name_transfer[key] for key in available_data]
		print fields
		todo_fields=set(fields)&set(csapr_standard_names().keys())
		print todo_fields
		flat_dict={}
		sample_volume=radarobj.contents.volumes[py4dd.fieldTypes().list.index(available_data[0])]
		#determine the min number of rays
		self.nsweeps=sample_volume.h.nsweeps
		rays=array([sample_volume.sweeps[i].h.nrays for i in range(self.nsweeps)])
		print rays.max(), rays.min()
		self.nrays=rays.min()#sample_volume.sweeps[0].h.nrays
		self.ngates=sample_volume.sweeps[0].rays[0].h.nbins
		print "hoo", sample_volume.sweeps[0].h.azimuth
		if sample_volume.sweeps[0].h.azimuth == -999.0:
			self.scan_type='ppi'
			self.naz=self.nrays
			self.nele=self.nsweeps
		else:
			self.scan_type='rhi'
			self.naz=self.nsweeps
			self.nele=self.nrays
		azimuth,elevation=self.extract_rsl_pointing(sample_volume)
		self.range={'data':sample_volume.sweeps[0].rays[0].dists, 'units':'meters', 'standard_name':'projection_range_coordinate', 'long_name':'range_to_measurement_volume', 'comment':'Coordinate variable for range. Range to center of each bin.', 'spacing_is_constant':'true', 'meters_to_center_of_first_gate':sample_volume.sweeps[0].rays[0].h.range_bin1, 'meters_between_gates':sample_volume.sweeps[0].rays[0].h.gate_size} #the range array which describes the range of all beams (note in this 
		self.azimuth={'data':azimuth.flatten(), 'units':'degrees', 'comment':'Azimuth of antenna relative to true north', 'long_name':'azimuth_angle_from_true_north', 'standard_name':'beam_azimuth_angle'} #The flat azimuth array which describes the azimuth of each beam
		self.elevation={'data':elevation.flatten(), 'units':'degrees', 'standard_name':'beam_elevation_angle', 'comment':'Elevation of antenna relative to the horizontal plane', 'long_name':'elevation_angle_from_horizontal_plane'}
		self.tu="seconds since %(year)d-%(month)02d-%(day)02d %(hour)02d:%(minute)02d:%(second)02d.0" % self.ray_header_time_to_dict(sample_volume.sweeps[0].rays[0].h)
		self.cal="gregorian"
		time_end=date2num(datetime(sample_volume.sweeps[-1].rays[-1].h.year, 
				sample_volume.sweeps[-1].rays[-1].h.month, 
				sample_volume.sweeps[-1].rays[-1].h.day, 
				sample_volume.sweeps[-1].rays[-1].h.hour, 
				sample_volume.sweeps[-1].rays[-1].h.minute,
				sample_volume.sweeps[-1].rays[-1].h.sec),self.tu, self.cal)	
		time_array=linspace(0, time_end, self.nrays*self.nsweeps)
		self.time={'data':time_array, 'units':self.tu, 'calendar':self.cal, 'comment':'Coordinate variable for time. Time at the center of each ray, in fractional seconds since the global variable time_coverage_start', 'standard_name':'time', 'long_name':'time in seconds since volume start'}
		for field in todo_fields: #create a dictionary tree for all data fields
			print "Doing ",field
			rsl_field=[key for key, value in name_transfer.iteritems() if value == field][0]
			print "Corresponds to ",  rsl_field
			data=create_cube_array_lim(radarobj.contents.volumes[py4dd.fieldTypes().list.index(rsl_field)], self.nsweeps, self.nrays)#radarobj.read_a_field(radarobj.fields.index(field)) #grab data from MDV object
			data[where(isnan(data))]=-9999.0
			data[where(data == 131072)]=-9999.0
			meta=self.get_mdv_meta(radarobj, field) #fetch metadata
			fielddict={'data':ma.masked_equal(data,-9999.0).reshape(data.shape[0]*data.shape[1], data.shape[2])} 
			fielddict.update(meta)
			flat_dict.update({csapr_standard_names()[field]:fielddict})
		self.fields=flat_dict
		if self.scan_type=='ppi':
			self.nsweeps=self.nele
			sweep_number={'data':range(self.nsweeps), 'units':'count', 'long_name':'sweep_number'}
			sweep_mode={'data':self.nsweeps*['azimuth_surveillance    '], 'long_name':'sweep_mode', 'units':'uniteless', 'comment':'Options are:"sector","coplane",rhi","vertical_pointing","idle","azimuth_surveillance","elevation_surveillance","sunscan","pointing","manual_ppi","manual_rhi"'}
			fixed_angle={'data':array([sample_volume.sweeps[i].h.elev for i in range(self.nsweeps)]), 'long_name':'target_angle_for_sweep', 'units':'degrees', 'standard_name':'target_fixed_angle'}
			sweep_start_ray_index={'data':arange(0,len(self.time['data']), self.naz), 'long_name':'index of first ray in sweep, 0-based', 'units':'count'}
			sweep_end_ray_index={'data':arange(self.naz-1,len(self.time['data']), self.naz), 'long_name':'index of last ray in sweep, 0-based', 'units':'count'}
		elif self.scan_type=='rhi':
			print "WASSL!"
			self.nsweeps=self.naz
			sweep_number={'data':range(self.nsweeps), 'units':'count', 'long_name':'sweep_number'}
			sweep_mode={'data':self.nsweeps*['rhi                     '], 'long_name':'sweep_mode', 'units':'uniteless', 'comment':'Options are:"sector","coplane",rhi","vertical_pointing","idle","azimuth_surveillance","elevation_surveillance","sunscan","pointing","manual_ppi","manual_rhi"'}
			fixed_angle={'data':array([sample_volume.sweeps[i].h.azimuth for i in range(self.nsweeps)]), 'long_name':'target_angle_for_sweep', 'units':'degrees', 'standard_name':'target_fixed_angle'}
			sweep_start_ray_index={'data':arange(0,len(self.time['data']), self.nele), 'long_name':'index of first ray in sweep, 0-based', 'units':'count'}
			sweep_end_ray_index={'data':arange(self.nele-1,len(self.time['data']), self.nele), 'long_name':'index of last ray in sweep, 0-based', 'units':'count'}
		self.sweep_info={'sweep_number':sweep_number, 'sweep_mode':sweep_mode, 'fixed_angle':fixed_angle, 'sweep_start_ray_index':sweep_start_ray_index, 'sweep_end_ray_index':sweep_end_ray_index}
		self.sweep_mode=array([self.scan_type]*self.nsweeps)
		self.sweep_number=linspace(0,self.nsweeps-1, self.nsweeps)
		metadata={'original_container':'rsl'}
		self.metadata=metadata
		#now for location variables
		radar_loc=[ dms_to_d((radarobj.contents.h.latd, radarobj.contents.h.latm, radarobj.contents.h.lats)) , dms_to_d((radarobj.contents.h.lond, radarobj.contents.h.lonm, radarobj.contents.h.lons)) ]
		lat={'data':radar_loc[0], 'standard_name':'Latitude', 'units':'degrees_north'}
		lon={'data':radar_loc[1], 'standard_name':'Longitude', 'units':'degrees_east'}
		elv={'data':radarobj.contents.h.height, 'standard_name':'Altitude', 'units':'meters'}
		self.location={'latitude':lat, 'longitude':lon, 'altitude':elv}
		#now for instrument parameters.. sorry but I am just going to brute force this!
		#prt mode: Need to fix this.. assumes dual if two prts 
		inst_params={'prt_mode':{'data':array([self.prtmode(sample_volume.sweeps[i].rays[0].h) for i in range(self.nsweeps)]), 'comments':'Pulsing mode Options are: "fixed", "staggered", "dual" Assumed "fixed" if missing.'}, 
		'prt':{'data':array([sample_volume.sweeps[i].rays[0].prf for i in range(self.nsweeps)]), 'units':'seconds', 'comments':"Pulse repetition time.For staggered prt, also see prt_ratio."},
		'unambiguous_range':{'data':array([sample_volume.sweeps[i].rays[0].unam_rng*1000.0 for i in range(self.nsweeps)]), 'units':'meters', 'comment':'Unambiguous range'}}
		self.inst_params=inst_params
	
	def mdv2rad(self, radarobj):
		#We only want to transfer fields that we have valid names for... 
		valid_fields=csapr_standard_names()
		todo_fields=set(radarobj.fields)&set(csapr_standard_names().keys()) #the intersection of the set of available and valid names
		flat_dict={} #a holder to be updated
		self.naz=len(radarobj.az_deg)
		self.nele=len(radarobj.el_deg)
		self.ngates=len(radarobj.range_km)
		self.azimuth={'data':tile(radarobj.az_deg, self.nele), 'units':'degrees', 'comment':'Azimuth of antenna relative to true north', 'long_name':'azimuth_angle_from_true_north', 'standard_name':'beam_azimuth_angle'} #The flat azimuth array which describes the azimuth of each beam
		self.range={'data':array(radarobj.range_km), 'units':'meters', 'standard_name':'projection_range_coordinate', 'long_name':'range_to_measurement_volume', 'comment':'Coordinate variable for range. Range to center of each bin.', 'spacing_is_constant':'true', 'meters_to_center_of_first_gate':(radarobj.range_km[0])/1000.0, 'meters_between_gates':(radarobj.range_km[1]- radarobj.range_km[0])/1000.0} #the range array which describes the range of all beams (note in this 
		self.elevation={'data':array(radarobj.el_deg).repeat( self.naz), 'units':'degrees', 'standard_name':'beam_elevation_angle', 'comment':'Elevation of antenna relative to the horizontal plane', 'long_name':'elevation_angle_from_horizontal_plane'}
		#append time
		tu="seconds since %(year)d-%(month)02d-%(day)02d %(hour)02d:%(minute)02d:%(second)02d.0" % dt_to_dict(radarobj.times['time_begin'])
		cal="gregorian"
		time_array=linspace(date2num(radarobj.times['time_begin'], tu, cal), date2num(radarobj.times['time_end'], tu, cal),self.naz*self.nele )
		self.time={'data':time_array, 'units':tu, 'calendar':cal, 'comment':'Coordinate variable for time. Time at the center of each ray, in fractional seconds since the global variable time_coverage_start', 'standard_name':'time', 'long_name':'time in seconds since volume start'}
		for field in todo_fields: #create a dictionary tree for all data fields
			print "Doing ",field 
			data=radarobj.read_a_field(radarobj.fields.index(field)) #grab data from MDV object
			data[where(isnan(data))]=-9999.0
			data[where(data == 131072)]=-9999.0
			meta=self.get_mdv_meta(radarobj, field) #fetch metadata
			fielddict={'data':ma.masked_equal(data,-9999.0).reshape(data.shape[0]*data.shape[1], data.shape[2])} 
			fielddict.update(meta)
			flat_dict.update({csapr_standard_names()[field]:fielddict})
		self.fields=flat_dict
		#sweep based stuff now:
		if radarobj.scan_type=='ppi':
			self.nsweeps=self.nele
			sweep_number={'data':range(self.nsweeps), 'units':'count', 'long_name':'sweep_number'}
			sweep_mode={'data':self.nsweeps*['azimuth_surveillance    '], 'long_name':'sweep_mode', 'units':'uniteless', 'comment':'Options are:"sector","coplane",rhi","vertical_pointing","idle","azimuth_surveillance","elevation_surveillance","sunscan","pointing","manual_ppi","manual_rhi"'}
			fixed_angle={'data':array(radarobj.el_deg), 'long_name':'target_angle_for_sweep', 'units':'degrees', 'standard_name':'target_fixed_angle'}
			sweep_start_ray_index={'data':arange(0,len(self.time['data']), self.naz), 'long_name':'index of first ray in sweep, 0-based', 'units':'count'}
			sweep_end_ray_index={'data':arange(self.naz-1,len(self.time['data']), self.naz), 'long_name':'index of last ray in sweep, 0-based', 'units':'count'}
		elif radarobj.scan_type=='rhi':
			self.nsweeps=self.naz
			sweep_number={'data':range(self.nsweeps), 'units':'count', 'long_name':'sweep_number'}
			sweep_mode={'data':self.nsweeps*['rhi                     '], 'long_name':'sweep_mode', 'units':'uniteless', 'comment':'Options are:"sector","coplane",rhi","vertical_pointing","idle","azimuth_surveillance","elevation_surveillance","sunscan","pointing","manual_ppi","manual_rhi"'}
			fixed_angle={'data':array(radarobj.az_deg), 'long_name':'target_angle_for_sweep', 'units':'degrees', 'standard_name':'target_fixed_angle'}
			sweep_start_ray_index={'data':arange(0,len(self.time['data']), self.nele), 'long_name':'index of first ray in sweep, 0-based', 'units':'count'}
			sweep_end_ray_index={'data':arange(self.nele-1,len(self.time['data']), self.nele), 'long_name':'index of last ray in sweep, 0-based', 'units':'count'}
		elif radarobj.scan_type=='vpr':
			self.nsweeps=1
		self.sweep_info={'sweep_number':sweep_number, 'sweep_mode':sweep_mode, 'fixed_angle':fixed_angle, 'sweep_start_ray_index':sweep_start_ray_index, 'sweep_end_ray_index':sweep_end_ray_index}
		self.sweep_mode=array([radarobj.scan_type]*self.nsweeps)
		self.sweep_number=linspace(0,self.nsweeps-1, self.nsweeps)
		mapme=defaut_mdv_metadata_map()
		metadata={}
		masterdata=dict([(key, radarobj.master_header[mapme[key]]) for key in mapme.keys()])
		metadata.update(masterdata)
		self.metadata=metadata
		#now for location variables
		lat={'data':radarobj.radar_info['latitude_deg'], 'standard_name':'Latitude', 'units':'degrees_north'}
		lon={'data':radarobj.radar_info['longitude_deg'], 'standard_name':'Longitude', 'units':'degrees_east'}
		elv={'data':radarobj.radar_info['altitude_km']*1000.0, 'standard_name':'Altitude', 'units':'meters'}
		self.location={'latitude':lat, 'longitude':lon, 'altitude':elv}
		#now for instrument parameters.. sorry but I am just going to brute force this!
		#prt mode: Need to fix this.. assumes dual if two prts 
		if radarobj.radar_info['prt2_s']==0.0:
			prt_mode='fixed'
		else:
			prt_mode='dual'
		inst_params={'prt_mode':{'data':array([prt_mode]*self.nele), 'comments':'Pulsing mode Options are: "fixed", "staggered", "dual" Assumed "fixed" if missing.'}, 
		'prt':{'data':array([radarobj.radar_info['prt_s']]*self.nele), 'units':'seconds', 'comments':"Pulse repetition time.For staggered prt, also see prt_ratio."},
		'unambiguous_range':{'data':array([radarobj.radar_info['unambig_range_km']*1000.0]*self.naz*self.nele), 'units':'meters', 'comment':'Unambiguous range'}}
		self.inst_params=inst_params
		
	def get_mdv_meta(self, radarobj, field):
		debug=True
		print "go"
		csapr_names=csapr_standard_names()
		moment_fixes={'DBZ_F':{'units':'dBZ', 'standard_name':'equivalent_reflectivity_factor','long_name':'equivalent_reflectivity_factor', 'valid_max':80.0, 'valid_min':-45.0, 'least_significant_digit':2},
		'VEL_F':{'units':'m/s', 'standard_name':'radial_velocity_of_scatterers_away_from_instrument','long_name':'radial_velocity_of_scatterers_away_from_instrument', 'valid_max':95.0, 'valid_min':-95.0, 'least_significant_digit':2},
		'KDP_F':{'units':'degrees/km', 'standard_name':'specific_differential_phase_hv', 'long_name':'specific_differential_phase_hv','valid_max':20.0, 'valid_min':-10.0, 'least_significant_digit':2},
		'ZDR_F':{'units':'dB', 'standard_name':'log_differential_reflectivity_hv','long_name':'log_differential_reflectivity_hv', 'valid_max':8.0, 'valid_min':-6.0, 'least_significant_digit':3},
		'RHOHV_F':{'units':'ratio', 'standard_name':'cross_correlation_ratio_hv','long_name':'cross_correlation_ratio_hv', 'valid_max':1.0, 'valid_min':0.0, 'least_significant_digit':5},
		'NCP_F':{'units':'ratio', 'standard_name':'signal_quality', 'long_name':'signal_quality','valid_max':1.0, 'valid_min':0.0, 'comment':'Also know as Normalized Coherent Power', 'least_significant_digit':5},
		'WIDTH_F':{'units':'m/s', 'standard_name':'spectrum_width','long_name':'spectrum_width', 'valid_max':45.0, 'valid_min':0.0, 'least_significant_digit':2},
		'PHIDP_F':{'units':'degrees', 'standard_name':'differential_phase_hv','long_name':'differential_phase_hv', 'valid_max':180.0, 'valid_min':-180.0, 'least_significant_digit':2},
		'VEL_COR':{'units':'m/s', 'standard_name':'radial_velocity_of_scatterers_away_from_instrument','long_name':'radial_velocity_of_scatterers_away_from_instrument', 'valid_max':45.0, 'valid_min':-45.0, 'least_significant_digit':2},
		'PHIDP_UNF':{'units':'degrees', 'standard_name':'differential_phase_hv','long_name':'differential_phase_hv', 'valid_max':480.0, 'valid_min':0.0, 'least_significant_digit':2},
		'DBZ_AC':{'units':'dBZ', 'standard_name':'equivalent_reflectivity_factor', 'long_name':'equivalent_reflectivity_factor','valid_max':80.0, 'valid_min':-45.0, 'least_significant_digit':2},
		'KDP_SOB':{'units':'degrees/km', 'standard_name':'specific_differential_phase_hv','long_name':'specific_differential_phase_hv', 'valid_max':20.0, 'valid_min':-1.0, 'least_significant_digit':3}}
		return moment_fixes[field]

