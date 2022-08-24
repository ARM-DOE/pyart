import numpy as np
from datatree import DataTree
import xarray as xr

def create_dataset_from_sweep(radar, sweep=0, field=None):
    """
    Creates an xarray.Dataset from sweeps
    """
    
    # Grab the range
    range = radar.range['data']
    
    # Grab the elevation
    elevation = radar.get_elevation(sweep)
    
    # Grab the azimuth
    azimuth = radar.get_azimuth(sweep)
    
    # Grab the x, y, z values
    x, y, z = radar.get_gate_x_y_z(sweep)
    
    # Grab the lat, lon, and elevation
    lat, lon, alt = radar.get_gate_lat_lon_alt(sweep)
    
    # Add the fields
    field_dict = {}
    for field in list(radar.fields):
        field_dict[field]=(["azimuth", "range"], radar.get_field(sweep, field))
    
    ds = xr.Dataset(
        data_vars=field_dict,
         coords=dict(
             x=(["azimuth", "range"], x),
             y=(["azimuth", "range"], y),
             z=(["azimuth", "range"], z),
             lat=(["azimuth", "range"], lat),
             lon=(["azimuth", "range"], lon),
             alt=(["azimuth", "range"], alt),
             azimuth=azimuth,
             elevation=elevation,
             sweep=np.array([sweep]),
             range=range),
    )
    
    return ds


def convert_to_xradar(radar):
    """
    Converts from radar to xradar
    """
    
    ds_list = []
    sweeps = radar.sweep_number['data']
    for sweep in sweeps:
        ds_list.append(create_dataset_from_sweep(radar, sweep))
    
    # Convert the numpy array to a list of strings
    dict_keys = [x for x in list(sweeps.astype(str))]
    
    # Zip the keys and datasets and turn into a dictionary
    dict_of_dsets = dict(zip(dict_keys, ds_list))
    
    return DataTree.from_dict(dict_of_dsets)