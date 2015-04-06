"""
ADD FILE DESCRIPTION
"""

import datetime
from netCDF4 import num2date, date2num
import numpy as np

import .mdv as MDV

def write_grid_mdv(filename, grid):
    """
    ADD FUNCTION DESCRIPTION HERE
    
    
    obs: xy grid must be regular
        data must be pre-encoded
        max nz =122
    """
    #XXX proper encolding (uint8, uint16 or float32) should be done before calling this function
    #first of all firm field list
    fields=grid.fields.keys()
    
    #mount empty mdv file
    grid_shape = grid.fields[fields[0]]['data'].shape
    nz, ny, nx = grid_shape
    nfields = len(fields)
    if nz > MDV.MDV_MAX_VLEVELS:
        import warnings
        warnings.warn('%i vlevels exceed MDV_MAX_VLEVELS = %i. Extra levels will be ignored'%(nz,MDV.MDV_MAX_VLEVELS))
        nz = MDV.MDV_MAX_VLEVELS
    mdv=MDV.MdvFile(None)
    mdv.field_headers = mdv._get_field_headers(nfields)
    mdv.vlevel_headers = msc._get_vlevel_headers(nfields)
    mdv.fields_data = [None]*nfields
    mdv.projection = 'flat'
    
    #fill headers
    d = mdv.master_header
    if "time_start" in grid.axes.keys():
        d["time_begin"] = time_dict_to_unixtime(grid.axes["time_start"])
    if "time_end" in grid.axes.keys():
        d["time_end"] = time_dict_to_unixtime(grid.axes["time_end"])
    if "time" in grid.axes.keys():
        d["time_centroid"] = time_dict_to_unixtime(grid.axes["time"])
    else:
        d["time_centroid"] = d["time_begin"]
    d["data_dimension"] = 3 #XXX are grid's always 3d?
    d["data_collection_type"] = 3  # DATA_SYNTHESIS, I don't realy know, so miscellaneous!
    if grid.axes['z_disp']['unit']=='m' or grid.axes['z_disp']['unit']=='meters':
        d["native_vlevel_type"] = 4
    elif grid.axes['z_disp']['unit']=='°' or grid.axes['z_disp']['unit']=='degree':
        d["native_vlevel_type"] = 9
    d["vlevel_type"] = d["native_vlevel_type"]
    d["nfields"] = nfields
    d["max_nx"] = nx
    d["max_ny"] = ny
    d["max_nz"] = nz
    d["time_written"] = (datetime.utcnow()-datetime.datetime(1970, 1, 1, 00, 00)).total_seconds()
    #try metadata, if not use axes
    if "radar_0_lon" in grid.metadata.keys():
        d["sensor_lon"] = grid.metadata["radar_0_lon"]
    elif 'lat' in grid.axes.keys:
        d["sensor_lon"] = grid.axes['lon']['data'][0]
    if "radar_0_lat" in grid.metadata.keys():
        d["sensor_lat"] = grid.metadata["radar_0_lat"]
    elif 'lon' in grid.axes.keys:
        d["sensor_lat"] = grid.axes['lat']['data'][0]
    if "radar_0_alt" in grid.metadata.keys():
        d["sensor_alt"] = grid.metadata["radar_0_alt"]
    elif 'alt' in grid.axes.keys:
        d["sensor_alt"] = grid.axes['alt']['data'][0]

    # there is no mandatory metadata to use in the following
    #d["data_set_info"] = 
    #d["data_set_name"] = 
    #d["data_set_source"] = 

    for ifield,field in enumerate(fields):
        d = mdv.field_headers[ifield]
        l = mdv.vlevels_headers[ifield]
        
        #fill fields_header
        d["nx"] = nx                 # 17i: 10
        d["ny"] = ny                  # 17i: 11
        d["nz"] = nz                # 17i: 12
        d["proj_type"] = MDV.PROJ_FLAT
        dtype = grid.fields[field]['data'].dtype
        if dtype == np.uint8:
            d["encoding_type"] = MDV.ENCODING_INT8
             d["data_element_nbytes"] = 1
        elif dtype == np.uint16:
            d["encoding_type"] = MDV.ENCODING_INT16
             d["data_element_nbytes"] = 2
        elif dtype == np.float32:
            d["encoding_type"] = MDV.ENCODING_FLOAT32
             d["data_element_nbytes"] = 4
        else:
            raise("Unsuported encoding, please encode data as uint8, uint16 or float32 before calling this function")
        d["compression_type"] = 3 #zlib

        d["scaling_type"] = 4 #SCALING_SPECIFIED (by the user)
        d["native_vlevel_type"] = mdv.master_header["vlevel_type"]     # 9i: 4
        d["vlevel_type"] = mdv.master_header["vlevel_type"]            # 9i: 5
        d["data_dimension"] = 3 #XXX are grid's always 3d?
        d["proj_origin_lat"] = grid.axes['lat']['data'][0]       # f
        d["proj_origin_lon"] = grid.axes['lon']['data'][0]        # f
        d["grid_dx"] = (grid.axes["x_disp"]['data'][1]-grid.axes["x_disp"]['data'][0])/1000.               # 12f: 2
        d["grid_dy"] = (grid.axes["y_disp"]['data'][1]-grid.axes["y_disp"]['data'][0])/1000.                # 12f: 3
        d["grid_minx"] = (grid.axes["x_disp"]['data'][0])/1000.             # 12f: 5
        d["grid_miny"] = (grid.axes["y_disp"]['data'][0])/1000.              # 12f: 6
        if "scale_factor" in grid.fields[field].keys():
            d["scale"] = grid.fields[field]["scale_factor"]                    # 12f: 8
        if "add_offset" in grid.fields[field].keys():
            d["bias"] = grid.fields[field]["add_offset"]  
        
        # bad_data prioritise _FillValue
        if "_FillValue" in grid.fields[field].keys():
            d["bad_data_value"] = grid.fields[field]["_FillValue"]        # 12f: 10
        elif "missing_value" in grid.fields[field].keys():
            d["bad_data_value"] = grid.fields[field]["missing_value"]        # 12f: 10
         # missing_data prioritise missing_value
        if "missing_value" in grid.fields[field].keys():
            d["missing_data_value"] = grid.fields[field]["missing_value"]        # 12f: 10
        elif "_FillValue" in grid.fields[field].keys():
            d["missing_data_value"] = grid.fields[field]["_FillValue"]        # 12f: 10
        
        d["min_value"] = np.amax(grid.fields[field]['data'])              # 5f: 1
        d["max_value"] = np.amin(grid.fields[field]['data'])              # 5f: 2
        if "standard_name" in grid.fields[field].keys():
            d["field_name_long"] = grid.fields[field]["standard_name"]
        elif "long_name" in  grid.fields[field].keys():
            d["field_name_long"] = grid.fields[field]["long_name"]
        d["field_name"] = field         # 16c
        if "units" in  grid.fields[field].keys():
            d["units"] = grid.fields[field]["units"]
        d["transform"] = "none" #XXX not implemented
        
        #fill vlevels_header
        typ=[0]*122
        level=[0]*122
        for iz in range(nz):
            typ[iz] = d["vlevel_type"]
            typ[iz] = grid.axes["z_disp"]["data"][iz]
        d["type"] = typ
        d["level"] = level
        
        #put data to field
        mdv.fields_data[ifield] = grid.fields[field]["data"]
    
    #write the file
    mdv.write(filename)
    
#XXX move to common
def time_dict_to_unixtime(d):
    """ convert a dict containing NetCDF style time information to unixtime, i.e second since 1st Jan 1970 00:00 """
    if 'calendar' in d
        calendar = d['calendar']
    else:
        calendar = 'standard'
    
    date = num2date(d['data'],d['units'],calendar)
    epoch= datetime.datetime(1970, 1, 1, 00, 00)
    return (date-epoch).total_seconds()
    
 
