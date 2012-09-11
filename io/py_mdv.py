""" Utilities for reading of MDV data into numpy objects

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

Code is adapted from Nitin Bharadwaj's Matlab code 

Scott Collis, Argonne National Laboratory, 
Computing Environment and Life Sciences Directorate
ARM Climate Research Facility

United States Department of Energy

USE
---
Currently this code works only on polar MDV files which are gzipped.. future versions will be expanded to load cartesian and deal with non-gzipped files

This will be updates with some examples

REQUIREMENTS
------------
Needs Numpy, datetime, gzip

HISTORY
-------
2011-05-23 Start of development 
Scott Collis scollis.acrf@gmail.com
0.1: basic functionality
2011-05-24
1.0: Added helper functions to calculate cartesian co-ordinates, forced the init function to get all headers on execute

"""
__author__ = "Scott Collis <scollis@anl.gov>"
__version__ = "1.0"

import struct
import gzip
import StringIO
import zlib
from numpy import array, float32, putmask, NaN, floor, arange, zeros, cos, sin, pi, arcsin, meshgrid, ones
import datetime

def radar_coords_to_cart(rng, az, ele, debug=False):
    """
    Asumes standard atmosphere, ie R=4Re/3
    """
    Re=6371.0*1000.0
    #h=(r^2 + (4Re/3)^2 + 2r(4Re/3)sin(ele))^1/2 -4Re/3
    #s=4Re/3arcsin(rcos(ele)/(4Re/3+h))
    p_r=4.0*Re/3.0
    rm=rng*1000.0
    z=(rm**2 + p_r**2 + 2.0*rm*p_r*sin(ele*pi/180.0))**0.5 -p_r
    #arc length
    s=p_r*arcsin(rm*cos(ele*pi/180.)/(p_r+z))
    if debug: print "Z=", z, "s=", s
    y=s*cos(az*pi/180.0)
    x=s*sin(az*pi/180.0)
    return x,y,z



class read_mdv:
	def __init__(self, filename, **kwargs):
		debug=kwargs.get('debug', False) #Noise on or off
		#The class initializes by loading the binary data
		if debug: print "Opening ", filename
		self.fileptr=open(filename, 'rb') #open and add the binary file to the class
		#Some information about the MDV file structure
		self.MDV_CHUNK_INFO_LEN=480
		self.MDV_INFO_LEN=512
		self.MDV_LONG_FIELD_LEN=64
		self.MDV_MAX_PROJ_PARAMS=8
		self.MDV_MAX_VLEVELS=122
		self.MDV_NAME_LEN=128
		self.MDV_SHORT_FIELD_LEN=16
		self.MDV_TRANSFORM_LEN=16
		self.MDV_UNITS_LEN=16
		self.MDV_N_COORD_LABELS=3
		self.MDV_COORD_UNITS_LEN=32
		# (x,y) in degrees. Simple latitude-longitude grid.
		# Also known as the Simple Cylindrical or Platte Carree projection.
		self.PROJ_LATLON = 0
		# (x,y) in km. Lambert Conformal Conic projection.
		self.PROJ_LAMBERT_CONF = 3
		# (x,y) in km. Polar Stereographic projection.
		self.PROJ_POLAR_STEREO = 5
		# Cartesian, (x,y) in km. This is a simple line-of-sight projection used for single
		# radar sites. The formal name is Oblique Lambert Azimuthal projection.
		self.PROJ_FLAT = 8
		# radar data in native Plan Position Indicator (PPI) coordinates of
		# range, azimuth angle and elevation angle. x is radial range (km), y is azimuth angle (deg), z is
		# elev angle (deg).
		self.PROJ_POLAR_RADAR = 9
		# (x,y) in km. Oblique Stereographic projection.
		self.PROJ_OBLIQUE_STEREO = 12
		# radar data in native Range Height Indicator (RHI) coordinates.
		# x is radial range (km), y is elev angle (deg), z is az angle (deg).
		self.PROJ_RHI_RADAR = 13
		#  ***************** COMPRESSION *******************
		self.COMPRESSION_NONE = 0# no compression
		self.COMPRESSION_ZLIB = 3# Lempel-Ziv
		self.COMPRESSION_BZIP = 4# bzip2
		self.COMPRESSION_GZIP = 5# Lempel-Ziv in gzip format
		
		#  ***************** COMPRESSION CODE *******************
		self.TA_NOT_COMPRESSED =791621423
		self.GZIP_COMPRESSED = 4160223223
		
		#  ***************** TRANSFORM *******************
		self.DATA_TRANSFORM_NONE = 0# None 
		self.DATA_TRANSFORM_LOG = 1# Natural log
		
		#  ***************** BIT ENCODING *******************
		self.ENCODING_INT8 = 1 # unsigned 8 bit integer
		self.ENCODING_INT16 = 2 # unsigned 16 bit integer
		self.ENCODING_FLOAT32 = 5 # 32 bit IEEE floating point
		
		#  ***************** CHUNK HEADER and DATA *******************
		self.CHUNK_DSRADAR_PARAMS = 3
		self.CHUNK_DSRADAR_ELEVATIONS=10
		self.CHUNK_DSRADAR_CALIB=7
		self.DS_LABEL_LEN=40
		self.NCHAR_DS_RADAR_PARAMS= 2*self.DS_LABEL_LEN
		self.DS_RADAR_CALIB_NAME_LEN=16
		self.DS_RADAR_CALIB_MISSING=-9999.0
		if debug: print "Getting master header"
		self.get_master_header()
		if debug: print "getting field headers"
		self.get_field_headers(self.master_header['nfields'])
		if debug: print "getting vlevel headers"
		self.get_vlevel_headers(self.master_header['nfields'])
		if debug: print "getting chunk headers"
		self.get_chunk_headers(self.master_header['nchunks'])
		if debug: print "Reading Chunks"
		self.read_chunks()
		if debug: print "Calculating Radar coordinates"
		self.calc_geometry()
		if debug: print "Making usable time objects"
		self.append_time()
		if debug: print "Calculating cartesian coordinates"
		self.append_carts()
		if debug: print "indexing fields"
		self.index_fields()


	
	def read_list_from_file(self, items, forms, endian, **kwargs):
		debug=kwargs.get('debug', False) #Noise on or off
		my_dict={}
		for i in range(len(forms)):
			nbytes=struct.calcsize(forms[i])
			if debug: print "reading ", nbytes, " of header data"
			bin_data=self.fileptr.read(nbytes)
			if 'c' in forms[i]:
				my_dict.update({items[i]: "".join(struct.unpack(endian+forms[i], bin_data)).split('\x00')[0]})
			else:
				my_dict.update({items[i]: struct.unpack(endian+forms[i], bin_data)})
				if len(my_dict[items[i]])==1:
					my_dict[items[i]]=my_dict[items[i]][0]
		return my_dict
	
	def close(self):
		#close the mdv file
		self.fileptr.close()
	
	def get_master_header(self):
		endian='>'
		items=["record_len1","struct_id","revision_number","time_gen","user_time","time_begin","time_end","time_centroid","time_expire","num_data_times","index_number","data_dimension","data_collection_type","user_data","native_vlevel_type","vlevel_type","vlevel_included","grid_orientation","data_ordering","nfields","ngates","nrays","nsweeps","nchunks","field_hdr_offset","vlevel_hdr_offset","chunk_hdr_offset","field_grids_differ","user_data_si328","time_written","unused_si325","user_data_fl326","sensor_lon","sensor_lat","sensor_alt","unused_fl3212","data_set_info","data_set_name","data_set_source","record_len2"]
		forms=['i','i','i','i','i','i','i','i','i','i','i','i','i','i','i','i','i','i','i','i','i','i','i','i','i','i','i', 'i', '8i', 'i', '5i', '6f', 'f', 'f', 'f', '12f', '%(n)dc' %{'n':self.MDV_INFO_LEN}, '%(n)dc' %{'n':self.MDV_NAME_LEN}, '%(n)dc' %{'n':self.MDV_NAME_LEN}, 'i']
		#print len(items)
		#print len(forms)
		#print forms
		my_dict=self.read_list_from_file(items, forms, '>')
		self.master_header=my_dict
	
	def get_field_headers(self, nfields):
		items=["record_len1", "struct_id",  "field_code",  "user_time1",  "forecast_delta", "user_time2",  "user_time3",   "forecast_time",   "user_time4",  "ngates",   "nrays",   "nsweeps",   "proj_type",   "encoding_type",   "data_element_nbytes",   "field_data_offset",    "volume_size",   "user_data_si32", "compression_type",  "transform_type",    "scaling_type",     "native_vlevel_type",   "vlevel_type",  "dz_constant",   "data_dimension",   "zoom_clipped",   "zoom_no_overlap",   "unused_si32",   "proj_origin_lat",   "proj_origin_lon",    "proj_param",    "vert_reference",    "grid_dx",     "grid_dy",    "grid_dz",    "grid_minx",   "grid_miny",    "grid_minz",    "scale",     "bias",      "bad_data_value",    "missing_data_value",   "proj_rotation",     "user_data_fl32",  "min_value",    "max_value",    "min_value_orig_vol",      "max_value_orig_vol",   "unused_fl32",  "field_name_long",  "field_name",  "units",   "transform",  "unused_char",  "record_len2"]
		forms=['i']*17 + ['10i'] + ['i']*9+['4i']+['f', 'f', '%(n)df' %{'n':self.MDV_MAX_PROJ_PARAMS}]+['f']*12+['4f']+['f']*5+['%(n)dc' %{'n':self.MDV_LONG_FIELD_LEN}, '%(n)dc' %{'n':self.MDV_SHORT_FIELD_LEN},'%(n)dc' %{'n':self.MDV_UNITS_LEN},'%(n)dc' %{'n':self.MDV_TRANSFORM_LEN},'%(n)dc' %{'n':self.MDV_UNITS_LEN}, 'i']
		fld_headers=[]
		for i in range(nfields):
			my_dict=self.read_list_from_file(items, forms, '>')
			fld_headers.append(my_dict)
		self.field_headers=fld_headers

	def get_vlevel_headers(self, nfields):
		items=["record_len1","struct_id","type","unused_si32","level","unused_fl32","record_len2"]
		forms=['i', 'i','%(n)di' %{'n':self.MDV_MAX_VLEVELS}, '4i', '%(n)df' %{'n':self.MDV_MAX_VLEVELS}, '5f', 'i']
		vl_headers=[]
		for i in range(nfields):
			my_dict=self.read_list_from_file(items, forms, '>')
			vl_headers.append(my_dict)
		self.vlevel_headers=vl_headers
	
	def get_chunk_headers(self, nchunks):
		items=["record_len1","struct_id","chunk_id","chunk_data_offset","size","unused_si32","info","record_len2"]
		forms=['i']*5 + ['2i', '%(n)dc' %{'n':self.MDV_CHUNK_INFO_LEN}, 'i']
		ch_headers=[]
		for i in range(nchunks):
			my_dict=self.read_list_from_file(items, forms, '>')
			ch_headers.append(my_dict)
		self.chunk_headers=ch_headers
		
	def get_compression_info(self):
		#we assume the file pointer is at the right spot!
		items=['magic_cookie','nbytes_uncompressed','nbytes_compressed','nbytes_coded','spare']
		forms=['I']*4+['2I']
		my_dict=self.read_list_from_file(items, forms, '>')
		self.current_compression_info=my_dict

	def read_a_field(self, fnum, **kwargs):
		debug=kwargs.get('debug', False)
		try:
			sweep_list= getattr(self, self.field_headers[fnum]['field_name'])
			if debug: print 'Getting data from object'
		except AttributeError:
			if debug: print 'No data found in object, populating'
			self.fileptr.seek(self.field_headers[fnum]['field_data_offset'])
			items=['vlevel_offsets', 'vlevel_nbytes']
			forms=['%(n)dI' %{'n':self.master_header['nsweeps']}]*2
			sweep_info_dict=self.read_list_from_file(items, forms, '>')
			sweep_list=zeros([self.master_header['nsweeps'], self.master_header['nrays'],self.master_header['ngates']], dtype=float32)
			for sw in range(self.master_header['nsweeps']):
				if debug: print "doint sweep ", sw
				self.get_compression_info()
				#print self.current_compression_info
				compressed_data=self.fileptr.read(self.current_compression_info['nbytes_coded'])
				#pp=open('/tmp/foo2.gz', 'wb')
				#pp.write(compressed_data)
				#pp.close()
				cd_fobj=StringIO.StringIO(compressed_data)
				gzip_file_handle = gzip.GzipFile(fileobj=cd_fobj)
				n_pts=self.master_header['ngates']*self.master_header['nrays']
				if self.field_headers[fnum]['encoding_type'] == self.ENCODING_INT8:
					form='>%(n)dB' %{'n':n_pts}
				elif self.field_headers[fnum]['encoding_type'] == self.ENCODING_INT16:
					form='>%(n)dH' %{'n':n_pts}
				elif self.field_headers[fnum]['encoding_type'] == self.ENCODING_FLOAT32:
					form='>%(n)df' %{'n':n_pts}
				#print "using: ", form
				decompressed_data = gzip_file_handle.read(struct.calcsize(form))
				gzip_file_handle.close()
				my_array=array(struct.unpack(form, decompressed_data)).reshape((self.master_header['nrays'],self.master_header['ngates'])).astype(float32)
				putmask(my_array, my_array==self.field_headers[fnum]['bad_data_value'], [NaN])
				ret_array=my_array*self.field_headers[fnum]['scale']+ self.field_headers[fnum]['bias']
				sweep_list[sw, :, :]=ret_array
			setattr(self, self.field_headers[fnum]['field_name'], sweep_list)
		return sweep_list
	
	def get_all_fields(self):
		my_list=[]
		for i in range(self.master_header['nfields']):
			print "reading field"
			my_list.append(self.read_a_field(i))
	
	def read_radar_info(self, nbytes):
		weareat=self.fileptr.tell()
		items=["radar_id","radar_type","nfields","ngates","samples_per_beam","scan_type","scan_mode","nfields_current","field_flag","polarization","follow_mode","prf_mode","spare_ints","radar_constant","altitude_km","latitude_deg","longitude_deg","gate_spacing_km","start_range_km","horiz_beam_width_deg","vert_beam_width_deg","pulse_width_us","prf_hz","wavelength_cm","xmit_peak_pwr_watts","receiver_mds_dbm","receiver_gain_db","antenna_gain_db","system_gain_db","unambig_vel_mps","unambig_range_km","measXmitPowerDbmH_dbm","measXmitPowerDbmV_dbm","prt_s","prt2_s","spare_floats","radar_name","scan_type_name"]
		forms=['i']*12 + ['2i'] +['f']*22 + ['4f', '%(n)dc' %{'n':self.DS_LABEL_LEN}, '%(n)dc' %{'n':self.DS_LABEL_LEN}]
		#for i in range(len(forms)):
		#	print items[i], ':', forms[i]
		radar_info_dict=self.read_list_from_file(items, forms, '>')
		self.radar_info=radar_info_dict
		chunk_read_bytes=self.fileptr.tell()-weareat
		#print self.fileptr.tell()-weareat
		#print nbytes-chunk_read_bytes
		self.fileptr.seek(nbytes-chunk_read_bytes,1)
	
	def get_elevs(self,nbytes):
		SIZE_FLOAT=4.0
		nelevations=floor(nbytes/SIZE_FLOAT)
		#print nelevations
		fid_start=self.fileptr.tell()
		form='%(n)df' %{'n':nelevations}
		eledata=self.fileptr.read(struct.calcsize(form))
		elevations=array(struct.unpack(form, eledata))
		self.elevations=elevations
		chunk_read_bytes=self.fileptr.tell()-fid_start
		#print chunk_read_bytes
		self.fileptr.seek(nbytes-chunk_read_bytes,1)
		#print nbytes-chunk_read_bytes
	
	def get_calib(self, nbytes):
		fid_start=self.fileptr.tell()
		items=["radar_name","year","month","day","hour","minute","second","wavelength_cm","beamwidth_h_deg","beamwidth_v_deg","antenna_gain_h_db","antenna_gain_v_db","pulse_width_us","xmit_power_h_dbm","xmit_power_v_dbm","twoway_waveguide_loss_h_db","twoway_waveguide_loss_v_db","twoway_radome_loss_h_db","twoway_radome_loss_v_db","filter_loss_db","radar_constant_h_db","radar_constant_v_db","noise_h_co_dbm","noise_h_cx_dbm","noise_v_co_dbm","noise_v_cx_dbm","rx_gain_h_co_dbm","rx_gain_h_cx_dbm","rx_gain_v_co_dbm","rx_gain_v_cx_dbm","zh1km_co_dbz","zh1km_cx_dbz","zv1km_co_dbz","zv1km_cx_dbz","sun_h_co_dbm","sun_h_cx_dbm","sun_v_co_dbm","sun_v_cx_dbm","noise_source_h_dbm","noise_source_v_dbm","power_meas_loss_h_db","power_meas_loss_v_db","coupler_fwd_loss_h_db","coupler_fwd_loss_v_db","zdr_bias_db","ldr_h_bias_db","ldr_v_bias_db","system_phidp_deg","test_pulse_h_dbm","test_pulse_v_dbm","rx_slope_h_co_db","rx_slope_h_cx_db","rx_slope_v_co_db","rx_slope_v_cx_db","I0_h_co_dbm","I0_h_cx_dbm","I0_v_co_dbm","I0_v_cx_dbm","spare"]
		forms=[ '%(n)dc' %{'n':self.DS_RADAR_CALIB_NAME_LEN}]+ ['i']*6 + ['f']*51 + ['14f']
		radar_calib_dict=self.read_list_from_file(items, forms, '>')
		self.calib_info=radar_calib_dict
		chunk_read_bytes=self.fileptr.tell()-fid_start
		self.fileptr.seek(nbytes-chunk_read_bytes,1)
		
	def read_chunks(self, **kwargs):
		debug=kwargs.get('debug', False)
		cal_info=[]
		self.fileptr.seek(self.chunk_headers[0]['chunk_data_offset'])
		for i in range(self.master_header['nchunks']):
			self.fileptr.seek(self.chunk_headers[i]['chunk_data_offset'])
			#print self.CHUNK_DSRADAR_ELEVATIONS, '==', self.chunk_headers[i]['chunk_id']
			if self.chunk_headers[i]['chunk_id']==self.CHUNK_DSRADAR_PARAMS:
				if debug:print 'Getting radar info'
				self.read_radar_info(self.chunk_headers[i]['size'])
			elif self.chunk_headers[i]['chunk_id']==self.CHUNK_DSRADAR_CALIB:#self.CHUNK_DSRADAR_ELEVATIONS:
				if debug:print 'getting elevations'
				self.get_elevs(self.chunk_headers[i]['size'])
			elif self.chunk_headers[i]['chunk_id']==self.CHUNK_DSRADAR_ELEVATIONS:#self.CHUNK_DSRADAR_CALIB:
				if debug:print 'getting cal'
				self.get_calib(self.chunk_headers[i]['size'])
	
	def calc_geometry(self):
		range_km=self.field_headers[1]['grid_minx']+arange(self.master_header['ngates'])*self.field_headers[1]['grid_dx']
		if self.field_headers[1]['proj_type']==self.PROJ_RHI_RADAR:
			self.scan_type='rhi'
			el_deg=self.field_headers[1]['grid_miny']+arange(self.master_header['nrays'])*self.field_headers[1]['grid_dy']
			az_deg=self.vlevel_headers[1]['level'][0:self.master_header['nsweeps']]
		if self.field_headers[1]['proj_type']==self.PROJ_POLAR_RADAR:
			self.scan_type='ppi'
			az_deg=self.field_headers[1]['grid_miny']+arange(self.master_header['nrays'])*self.field_headers[1]['grid_dy']
			el_deg=self.vlevel_headers[1]['level'][0:self.master_header['nsweeps']]
		self.az_deg=az_deg
		self.range_km=range_km
		self.el_deg=el_deg
	
	def append_time(self):
		time_objects=['time_begin', 'time_end', 'time_centroid']
		time_dict={}
		base_time=datetime.datetime(1970, 1, 1, 00, 00)
		for obj in time_objects:
			time_dict.update({obj: base_time+datetime.timedelta(seconds=self.master_header[obj])})
		self.times=time_dict
	
	def append_carts(self):
		#simple calculation involving 4/3 earth radius
		xx=zeros([self.master_header['nsweeps'], self.master_header['nrays'],self.master_header['ngates']], dtype=float32)
		yy=zeros([self.master_header['nsweeps'], self.master_header['nrays'],self.master_header['ngates']], dtype=float32)
		zz=zeros([self.master_header['nsweeps'], self.master_header['nrays'],self.master_header['ngates']], dtype=float32)
		if self.scan_type=='rhi':
			rg,ele=meshgrid(self.range_km, self.el_deg)
			for aznum in range(self.master_header['nsweeps']):
				azg=ones(rg.shape, dtype=float32)*self.az_deg[aznum]
				x,y,z=radar_coords_to_cart(rg, azg, ele)
				zz[aznum, :, :]=z
				xx[aznum, :, :]=x
				yy[aznum, :, :]=y
		if self.scan_type=='ppi':
			rg,azg=meshgrid(self.range_km, self.az_deg)
			for elnum in range(self.master_header['nsweeps']):
				ele=ones(rg.shape, dtype=float32)*self.el_deg[elnum]
				x,y,z=radar_coords_to_cart(rg, azg, ele)
				zz[elnum, :, :]=z
				xx[elnum, :, :]=x
				yy[elnum, :, :]=y
		self.carts={'x':xx, 'y':yy, 'z':zz}
	
	def index_fields(self):
		flds=[self.field_headers[i]['field_name'] for i in range(len(self.field_headers))]
		self.fields=flds

