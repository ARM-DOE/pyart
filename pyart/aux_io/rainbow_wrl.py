"""
pyart.aux_io.rainbow
====================

Routines for reading RAINBOW files (Used by SELEX) using the wradlib library

This routine has been tested to read rainbow5 files version 5.22.3, 5.34.16 and 5.35.1.
Since the rainbow file format is evolving constanly there is no guaranty that it can work with other versions.
If necessary, the user should adapt to code according to its own version.

Data types read by this routine:
Reflectivity: dBZ, dBuZ, dBZv, dBuZv
Velocity: V, Vu, Vv, Vvu
Spectrum width: W, Wu, Wv, Wvu
Differential reflectivity: ZDR, ZDRu
Co-polar correlation coefficient: RhoHV, RhoHVu
Co-polar differential phase: PhiDP, uPhiDP, uPhiDPu
Specific differential phase: KDP, uKDP, uKDPu
Signal quality parameters: SQI, SQIu, SQIv, SQIvu
Temperature: TEMP
Position of the range bin respect to the ISO0: ISO0

History

    V0.1 20160629 Jordi Figueras i Ventura, MeteoSwiss, first prototype
        
    read_rainbow_wrl

"""

# specific modules for this function
import os
import wradlib as wrl

import datetime

import numpy as np

from ..config import FileMetadata, get_fillvalue
from ..io.common import make_time_unit_str, _test_arguments
from ..core.radar import Radar
from ..exceptions import MissingOptionalDependency

RAINBOW_FIELD_NAMES = {
    'W': 'spectrum_width',
    'Wv': 'spectrum_width_vv', # non standard name
    'Wu': 'unfiltered_spectrum_width', # non standard name
    'Wvu': 'unfiltered_spectrum_width_vv', # non standard name
    'V': 'velocity',
    'Vv': 'velocity_vv', # non standard name
    'Vu': 'unfiltered_velocity', # non standard name
    'Vvu': 'unfiltered_velocity_vv', # non standard name
    'dBZ': 'reflectivity',
    'dBZv': 'reflectivity_vv',      # non standard name
    'dBuZ': 'unfiltered_reflectivity', # non standard name
    'dBuZv': 'unfiltered_reflectivity_vv', # non standard name
    'ZDR': 'differential_reflectivity',
    'ZDRu': 'unfiltered_differential_reflectivity', # non standard name
    'RhoHV': 'cross_correlation_ratio',
    'RhoHVu': 'unfiltered_cross_correlation_ratio', # non standard name
    'PhiDP': 'differential_phase',
    'uPhiDP': 'uncorrected_differential_phase', # non standard name
    'uPhiDPu': 'uncorrected_unfiltered_differential_phase', # non standard name 
    'KDP': 'specific_differential_phase',
    'uKDP': 'uncorrected_specific_differential_phase', # non standard name
    'uKDPu': 'uncorrected_unfiltered_specific_differential_phase', # non standard name
    'SQI': 'signal_quality_index', # non standard name
    'SQIv': 'signal_quality_index_vv', # non standard name
    'SQIu': 'unfiltered_signal_quality_index', # non standard name
    'SQIvu': 'unfiltered_signal_quality_index_vv', # non standard name
    'TEMP': 'temperature', # non standard name
    'ISO0': 'iso0', # non standard name
}


def read_rainbow_wrl(filename, field_names=None, additional_metadata=None, file_field_names=False, exclude_fields=None, **kwargs):
    """
    Read a RAINBOW file.
    
    Parameters
    ----------
    filename : str
        Name of the RAINBOW file to read.
    field_names : dict, optional
        Dictionary mapping RAINBOW field names to radar field names. If a
        data type found in the file does not appear in this dictionary or has
        a value of None it will not be placed in the radar.fields dictionary.
        A value of None, the default, will use the mapping defined in the
        Py-ART configuration file.
    additional_metadata : dict of dicts, optional
        Dictionary of dictionaries to retrieve metadata during this read.
        This metadata is not used during any successive file reads unless
        explicitly included.  A value of None, the default, will not
        introduct any addition metadata and the file specific or default
        metadata as specified by the Py-ART configuration file will be used.
    file_field_names : bool, optional
        True to use the MDV data type names for the field names. If this
        case the field_names parameter is ignored. The field dictionary will
        likely only have a 'data' key, unless the fields are defined in
        `additional_metadata`.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters.
    
    
    Returns
    -------
    radar : Radar
        Radar object containing data from RAINBOW file.
    
    """
    
    # check if it is the right file. Open it and read it
    bfile = os.path.basename(filename)
    if bfile.endswith('.vol') or bfile.endswith('.azi') or bfile.endswith('.ele'):
    
        # create metadata retrieval object
        if field_names is None:
            field_names = RAINBOW_FIELD_NAMES
        filemetadata = FileMetadata('RAINBOW', field_names, additional_metadata,
                                file_field_names, exclude_fields)
                                
        rbf=wrl.io.read_Rainbow(filename, loaddata=True)
            
        # check if file contains a known field
        # all slices should have the same data type
        nslices=int(rbf['volume']['scan']['pargroup']['numele'])
        if nslices > 1:
            datatype=rbf['volume']['scan']['slice'][0]['slicedata']['rawdata']['@type']
        else:
            datatype=rbf['volume']['scan']['slice']['slicedata']['rawdata']['@type']
            
        field_name = filemetadata.get_field_name(datatype)
        if field_name is not None:
        
            # get definitions from filemetadata class
            latitude = filemetadata('latitude')
            longitude = filemetadata('longitude')
            altitude = filemetadata('altitude')     
            metadata = filemetadata('metadata')
            sweep_start_ray_index = filemetadata('sweep_start_ray_index')
            sweep_end_ray_index = filemetadata('sweep_end_ray_index')
            sweep_number = filemetadata('sweep_number')
            sweep_mode = filemetadata('sweep_mode')
            fixed_angle = filemetadata('fixed_angle')
            elevation = filemetadata('elevation')
            _range = filemetadata('range')
            azimuth = filemetadata('azimuth')
            _time = filemetadata('time')
            field_dic = filemetadata(field_name)
            
            # other metadata
            scan_rate=filemetadata('scan_rate')
            frequency=filemetadata('frequency')
            
            # get general file information          
                
            # position
            if 'sensorinfo' in rbf['volume'].keys() :
                latitude['data']=np.array([rbf['volume']['sensorinfo']['lat']], dtype='float64')
                longitude['data']=np.array([rbf['volume']['sensorinfo']['lon']], dtype='float64')
                altitude['data']=np.array([rbf['volume']['sensorinfo']['alt']], dtype='float64')
                frequency['data']=3e8/ float(rbf['volume']['sensorinfo']['wavelen'])
            elif 'radarinfo' in rbf['volume'].keys() :
                latitude['data']=np.array([rbf['volume']['radarinfo']['@lat']], dtype='float64')
                longitude['data']=np.array([rbf['volume']['radarinfo']['@lon']], dtype='float64')
                altitude['data']=np.array([rbf['volume']['radarinfo']['@alt']], dtype='float64')
                frequency['data']=3e8/ float(rbf['volume']['sensorinfo']['wavelen'])
            
            # get number of rays and number of range bins per sweep         
            if nslices > 1:
                # sweep_number (is the sweep index)         
                sweep_number['data'] = np.arange(nslices, dtype='int32')
                
                rays_per_sweep=np.empty(nslices, dtype='int32')  # number of angles per ray         
                nbins=np.empty(nslices, dtype='int32')  # number of range bins per ray in sweep 
                for i in range(nslices):
                    # number of rays per sweep              
                    rays_per_sweep[i]=int(rbf['volume']['scan']['slice'][i]['slicedata']['rawdata']['@rays'])
                
                    # number of range bins per ray in sweep
                    nbins[i]=int(rbf['volume']['scan']['slice'][i]['slicedata']['rawdata']['@bins'])
            
                # all sweeps have to have the same number of range bins
                if any(nbins != nbins[0]):
                    raise ValueError('number of range bins changes between sweeps')
                nbins=nbins[0]
                
                # total number of rays
                total_rays = sum(rays_per_sweep)
                
                # sweep start ray index, sweep end ray index
                ssri = np.cumsum(np.append([0], rays_per_sweep[:-1])).astype('int32')
                seri = np.cumsum(rays_per_sweep).astype('int32') - 1
                sweep_start_ray_index['data'] = ssri
                sweep_end_ray_index['data'] = seri
                
                # angle step
                angle_step=float(rbf['volume']['scan']['slice'][0]['anglestep'])
                
                # range
                r_res=float(rbf['volume']['scan']['slice'][0]['rangestep'])*1000. # range resolution in m
                if 'start_range' in rbf['volume']['scan']['slice'][0].keys() :
                    start_range=float(rbf['volume']['scan']['slice'][0]['start_range'])*1000. # range start in m
                else:
                    start_range=0.
                _range['data'] = np.linspace(start_range+r_res/2., float(nbins-1.)*r_res+r_res/2., nbins, dtype='float32')
                
                # antenna speed
                if 'antspeed' in rbf['volume']['scan']['slice'][0]:
                    ant_speed=float(rbf['volume']['scan']['slice'][0]['antspeed'])
                else:
                    ant_speed=10.
                    print('WARNING: Unable to read antenna speed. Default value of '+str(ant_speed)+' deg/s will be used')
              
            else:
                # sweep_number (is the sweep index)         
                sweep_number['data'] = np.array([0], dtype='int32')
                
                rays_per_sweep=np.array([rbf['volume']['scan']['slice']['slicedata']['rawdata']['@rays']], dtype='int32')
                nbins=int(rbf['volume']['scan']['slice']['slicedata']['rawdata']['@bins'])
                total_rays=rays_per_sweep[0]
                
                ssri=0
                seri=total_rays-1
                sweep_start_ray_index['data']=[ssri]
                sweep_end_ray_index['data']=[seri]
                
                # angle step
                angle_step=float(rbf['volume']['scan']['slice']['anglestep'])
                
                # range
                r_res=float(rbf['volume']['scan']['slice']['rangestep'])*1000. # range resolution in m
                if 'start_range' in rbf['volume']['scan']['slice'].keys() :
                    start_range=float(rbf['volume']['scan']['slice']['start_range'])*1000. # range start in m
                else:
                    start_range=0.
                _range['data'] = np.linspace(start_range+r_res/2., float(nbins-1.)*r_res+r_res/2., nbins, dtype='float32')
                
                # antenna speed
                if 'antspeed' in rbf['volume']['scan']['slice']:
                    ant_speed=float(rbf['volume']['scan']['slice']['antspeed'])
                else:
                    ant_speed=10.
                    print('WARNING: Unable to read antenna speed. Default value of '+str(ant_speed)+' deg/s will be used')
                
            if bfile.endswith('.vol') or bfile.endswith('.azi') :               
                # sweep_mode            
                sweep_mode['data'] = np.array(nslices * ['azimuth_surveillance'])
                
                # scan_type    
                scan_type = 'ppi'
                
                # read data from file
                if nslices> 1:
                    # containers for data
                    t_fixed_angle=np.empty(nslices, dtype='float64')
                    angle_data=np.empty(total_rays, dtype='float64')
                    el_data=np.empty(total_rays, dtype='float64')
                    time_data=np.empty(total_rays, dtype='float64')
                    fdata = np.ma.zeros((total_rays, nbins), dtype='float32')
                    for i in range(nslices):                
                        # fixed angle                   
                        t_fixed_angle[i]=float(rbf['volume']['scan']['slice'][i]['posangle'])
    
                        # elevation
                        el_data[ssri[i]:seri[i]+1]=t_fixed_angle[i]
                        
                        # azimuth                       
                        if (len(rbf['volume']['scan']['slice'][i]['slicedata']['rayinfo']) == 2):
                            angle_start_bin=rbf['volume']['scan']['slice'][i]['slicedata']['rayinfo'][0]['data']
                            angle_start=np.array(angle_start_bin/65536.*360., dtype='float64')
                            
                            angle_stop_bin=rbf['volume']['scan']['slice'][i]['slicedata']['rayinfo'][1]['data']                                     
                            angle_stop=np.array(angle_stop_bin/65536.*360., dtype='float64')
                        else:
                            angle_start_bin=rbf['volume']['scan']['slice'][i]['slicedata']['rayinfo']['data']
                            angle_start=np.array(angle_start_bin/65536.*360., dtype='float64')
                            
                            angle_stop=angle_start+angle_step
                            
                        angle_data[ssri[i]:seri[i]+1]=np.angle( (np.exp(1.j*np.deg2rad(angle_start)) + np.exp(1.j*np.deg2rad(angle_stop)) ) / 2., deg=True)
                        
                        # time                  
                        date_sweep=rbf['volume']['scan']['slice'][i]['slicedata']['@date']
                        time_sweep=rbf['volume']['scan']['slice'][i]['slicedata']['@time']
                        datetime_sweep = datetime.datetime.strptime(date_sweep+' '+time_sweep, '%Y-%m-%d %H:%M:%S')
                        sweep_start_epoch = (datetime_sweep - datetime.datetime(1970, 1, 1)).total_seconds()
                                            
                        if (angle_stop[-1] > angle_start[0]) and ((angle_stop[-1]-angle_start[0])/rays_per_sweep[i] > angle_step):
                            sweep_duration=(angle_stop[-1]-angle_start[0])/ant_speed
                        else :
                            sweep_duration=(angle_stop[-1]+360.-angle_start[0])/ant_speed                   
                        time_angle=sweep_duration/rays_per_sweep[i]
                        
                        sweep_end_epoch=sweep_start_epoch+sweep_duration
                        if i == 0:  
                            volume_start_epoch=sweep_start_epoch+0.
                            start_time = datetime.datetime.utcfromtimestamp(volume_start_epoch)
                        time_data[ssri[i]:seri[i]+1]=np.linspace(sweep_start_epoch+time_angle/2., sweep_end_epoch-time_angle/2., num=rays_per_sweep[i])
                                            
                        # data                  
                        sweep_data_bin=rbf['volume']['scan']['slice'][i]['slicedata']['rawdata']['data']
                        datamin=float(rbf['volume']['scan']['slice'][i]['slicedata']['rawdata']['@min'])
                        datamax=float(rbf['volume']['scan']['slice'][i]['slicedata']['rawdata']['@max'])
                        datadepth=float(rbf['volume']['scan']['slice'][i]['slicedata']['rawdata']['@depth'])
                        sweep_data=np.array(datamin+sweep_data_bin*(datamax-datamin)/2 ** datadepth, dtype='float32')
                        
                        # fill invalid data with fill value
                        is_invalid= sweep_data_bin == 0
                        sweep_data[is_invalid.nonzero()]= get_fillvalue()
                        
                        # put phidp data in the range [-180, 180]
                        if (datatype == 'PhiDP') or (datatype == 'uPhiDP') or (datatype == 'uPhiDPu'):
                            is_above_180= sweep_data > 180.
                            sweep_data[is_above_180.nonzero()]-=360.
                        
                        sweep_data=np.reshape(sweep_data, [rays_per_sweep[i], nbins])
                        fdata[ssri[i]:seri[i]+1, :] = sweep_data
                else:
                    # fixed angle
                    t_fixed_angle=float(rbf['volume']['scan']['slice']['posangle'])
                    fixed_angle['data']=np.array([t_fixed_angle], dtype='float64')
                                    
                    # elevation
                    el_data=np.empty(total_rays, dtype='float64')
                    el_data[:]=t_fixed_angle                    
                    
                    # azimuth                   
                    if (len(rbf['volume']['scan']['slice']['slicedata']['rayinfo']) == 2):
                        angle_start_bin=rbf['volume']['scan']['slice']['slicedata']['rayinfo'][0]['data']
                        angle_start=np.array(angle_start_bin/65536.*360., dtype='float64')
                        
                        angle_stop_bin=rbf['volume']['scan']['slice']['slicedata']['rayinfo'][1]['data']                                
                        angle_stop=np.array(angle_stop_bin/65536.*360., dtype='float64')
                    else:
                        angle_start_bin=rbf['volume']['scan']['slice']['slicedata']['rayinfo']['data']
                        angle_start=np.array(angle_start_bin/65536.*360., dtype='float64')
                        
                        angle_stop=angle_start+angle_step
                        
                    angle_data=np.angle( (np.exp(1.j*np.deg2rad(angle_start)) + np.exp(1.j*np.deg2rad(angle_stop)) ) / 2., deg=True)
                    
                    # time
                    date_sweep=rbf['volume']['scan']['slice']['slicedata']['@date']
                    time_sweep=rbf['volume']['scan']['slice']['slicedata']['@time']
                    datetime_sweep = datetime.datetime.strptime(date_sweep+' '+time_sweep, '%Y-%m-%d %H:%M:%S')
                    sweep_start_epoch = (datetime_sweep - datetime.datetime(1970, 1, 1)).total_seconds()
                                                            
                    if angle_stop[-1] > angle_start[0] and ((angle_stop[-1]-angle_start[0])/rays_per_sweep[0] > angle_step):
                        sweep_duration=(angle_stop[-1]-angle_start[0])/ant_speed
                    else :
                        sweep_duration=(angle_stop[-1]+360.-angle_start[0])/ant_speed
                    time_angle=sweep_duration/rays_per_sweep[0]
                    
                    sweep_end_epoch=sweep_start_epoch+sweep_duration
                    volume_start_epoch=sweep_start_epoch+0.
                    start_time = datetime.datetime.utcfromtimestamp(volume_start_epoch)
                    time_data=np.array(np.linspace(sweep_start_epoch+time_angle/2., sweep_end_epoch-time_angle/2., num=rays_per_sweep[0]), dtype='float64')
                        
                    # data
                    sweep_data_bin=rbf['volume']['scan']['slice']['slicedata']['rawdata']['data']
                    datamin=float(rbf['volume']['scan']['slice']['slicedata']['rawdata']['@min'])
                    datamax=float(rbf['volume']['scan']['slice']['slicedata']['rawdata']['@max'])
                    datadepth=float(rbf['volume']['scan']['slice']['slicedata']['rawdata']['@depth'])
                    fdata=np.array(datamin+sweep_data_bin*(datamax-datamin)/2 ** datadepth, dtype='float32')
                    
                    # fill invalid data with fill value
                    is_invalid= sweep_data_bin == 0
                    fdata[is_invalid.nonzero()]= get_fillvalue()
                    
                    # put phidp data in the range [-180, 180]
                    if (datatype == 'PhiDP') or (datatype == 'uPhiDP') or (datatype == 'uPhiDPu'):
                        is_above_180= fdata > 180.
                        fdata[is_above_180.nonzero()]-=360.
                    
                    fdata=np.reshape(fdata, [rays_per_sweep[0], nbins])
                    
                fixed_angle['data']=np.array(t_fixed_angle, dtype='float64')
                azimuth['data'] = angle_data
                elevation['data'] = el_data
                
                _time['data']=time_data-volume_start_epoch
                _time['units'] = make_time_unit_str(start_time)
                
                # fields
                fields = {}                
                # create field dictionary                                           
                field_dic['_FillValue'] = get_fillvalue()
                field_dic['data'] = fdata
                fields[field_name] = field_dic          
            else :              
                # sweep_mode            
                sweep_mode['data'] = np.array(['elevation_surveillance'])
                
                # scan_type    
                scan_type = 'rhi'
                
                # fixed angle
                t_fixed_angle=float(rbf['volume']['scan']['slice']['posangle'])
                fixed_angle['data']=np.array([t_fixed_angle], dtype='float64')
                                
                # azimuth
                az_data=np.empty(total_rays, dtype='float64')
                az_data[:]=t_fixed_angle
                azimuth['data'] = az_data
                
                # elevation
                if (len(rbf['volume']['scan']['slice']['slicedata']['rayinfo']) == 2):
                    angle_start_bin=rbf['volume']['scan']['slice']['slicedata']['rayinfo'][0]['data']
                    angle_start=np.array(angle_start_bin/65536.*360., dtype='float64')
                    ind=(angle_start>225.).nonzero()
                    angle_start[ind]-=360.
                
                    angle_stop_bin=rbf['volume']['scan']['slice']['slicedata']['rayinfo'][1]['data']
                    angle_stop=np.array(angle_stop_bin/65536.*360., dtype='float64')
                    ind=(angle_stop>225.).nonzero()
                    angle_stop[ind]-=360.
                else:
                    angle_start_bin=rbf['volume']['scan']['slice']['slicedata']['rayinfo']['data']
                    angle_start=np.array(angle_start_bin/65536.*360., dtype='float64')
                    ind=(angle_start>225.).nonzero()
                    angle_start[ind]-=360.
                    
                    angle_stop=angle_start+angle_step
                    
                elevation['data'] = np.angle( (np.exp(1.j*np.deg2rad(angle_start)) + np.exp(1.j*np.deg2rad(angle_stop)) ) / 2., deg=True)
                
                # time
                date_sweep=rbf['volume']['scan']['slice']['slicedata']['@date']
                time_sweep=rbf['volume']['scan']['slice']['slicedata']['@time']
                datetime_sweep = datetime.datetime.strptime(date_sweep+' '+time_sweep, '%Y-%m-%d %H:%M:%S')
                sweep_start_epoch = (datetime_sweep - datetime.datetime(1970, 1, 1)).total_seconds()
                                    
                if angle_stop[-1] > angle_start[0]:
                    sweep_duration=(angle_stop[-1]-angle_start[0])/ant_speed
                else :
                    sweep_duration=(angle_start[0]-angle_stop[-1])/ant_speed
                time_angle=sweep_duration/rays_per_sweep[0]
                
                sweep_end_epoch=sweep_start_epoch+sweep_duration
                volume_start_epoch=sweep_start_epoch+0.
                start_time = datetime.datetime.utcfromtimestamp(volume_start_epoch)
                time_data=np.array(np.linspace(sweep_start_epoch+time_angle/2., sweep_end_epoch-time_angle/2., num=rays_per_sweep[0]), dtype='float64')
                _time['data']=time_data-volume_start_epoch
                _time['units'] = make_time_unit_str(start_time) 
                    
                # data
                sweep_data_bin=rbf['volume']['scan']['slice']['slicedata']['rawdata']['data']
                datamin=float(rbf['volume']['scan']['slice']['slicedata']['rawdata']['@min'])
                datamax=float(rbf['volume']['scan']['slice']['slicedata']['rawdata']['@max'])
                datadepth=float(rbf['volume']['scan']['slice']['slicedata']['rawdata']['@depth'])
                sweep_data=np.array(datamin+sweep_data_bin*(datamax-datamin)/2 ** datadepth, dtype='float32')
                
                # fill invalid data with fill value
                is_invalid= sweep_data_bin == 0
                sweep_data[is_invalid.nonzero()]= get_fillvalue()
                
                # put phidp data in the range [-180, 180]
                if (datatype == 'PhiDP') or (datatype == 'uPhiDP') or (datatype == 'uPhiDPu'):
                    is_above_180= sweep_data > 180.
                    sweep_data[is_above_180.nonzero()]-=360.
                
                sweep_data=np.reshape(sweep_data, [rays_per_sweep[0], nbins])
                
                # fields
                fields = {}
                # create field dictionary                                           
                field_dic['_FillValue'] = get_fillvalue()
                field_dic['data'] = sweep_data
                fields[field_name] = field_dic
        
#           # metadata
#           metadata['instrument_name']=radar_id

                

            # instrument_parameters
            instrument_parameters=dict()
            instrument_parameters.update({'frequency' : frequency})
                
            return Radar(_time, _range, fields, metadata, scan_type, latitude, longitude, altitude, sweep_number, 
            sweep_mode, fixed_angle, sweep_start_ray_index, sweep_end_ray_index, azimuth, elevation, instrument_parameters=instrument_parameters)
        else:
            raise ValueError('Field Name Unknown')
    else:
        raise ValueError('Only data files with extension .vol, .azi or .ele are supported')






