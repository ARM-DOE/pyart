"""
Calculate the composite reflectivity 

"""

import copy

import numpy as np
from netCDF4 import num2date
from pandas import to_datetime


def composite_reflectivity(radar, field='reflectivity', gatefilter=None):

    """
    Composite Reflectivity 
    
    Often a legacy product, composite reflectivity is: 
    "A display or mapping of the maximum radar reflectivity factor at any
    altitude as a function of position on the ground." - AMS Glossary
    This is more useful for the dry regions of the world, where maximum
    reflectivity values are found aloft, as opposed to the lowest scan.
    Alternatively this is useful for comparing to NWP since composite Z
    is a standard output of NWP. 

    Why this is not as easy as one would think: Turns out the data are
    not natively stored with index 0 being azimuth 0. Likely due to the
    physical spinning of the radar antenna. 

    Author: Randy J. Chase (@dopplerchase) 

    Parameters
    ----------
    radar : Radar
        Radar object used.
    gatefilter : GateFilter
        GateFilter instance. None will result in no gatefilter mask being
        applied to data. 

    Returns
    -------
    radar : Radar
        The radar object containing the radar dimensions, metadata and
        composite field.

    """
    # loop over all measured sweeps 
    for sweep in radar.sweep_number['data']:

        # get start and stop index numbers
        s_idx = radar.sweep_start_ray_index['data'][sweep]
        e_idx = radar.sweep_end_ray_index['data'][sweep] + 1

        # grab radar data 
        z = radar.get_field(sweep, field)

        # Use gatefilter
        if gatefilter is not None:
            mask_sweep = gatefilter.gate_excluded[s_idx:e_idx, :]
            z = np.ma.masked_array(z, mask_sweep)

        # extract lat lons 
        lon = radar.gate_longitude['data'][s_idx:e_idx, :]
        lat = radar.gate_latitude['data'][s_idx:e_idx, :]

        # get azimuth
        az = radar.azimuth['data'][s_idx:e_idx]
        # get order of azimuths 
        az_ids = np.argsort(az)

        # reorder azs so they are in order 
        z = z[az_ids]
        lon = lon[az_ids]
        lat = lat[az_ids]
 
        # if the first sweep, store re-ordered lons/lats
        if sweep == 0:
            lon_0 = copy.deepcopy(lon)
            lon_0[-1, :] = lon_0[0, :]
            lat_0 = copy.deepcopy(lat)
            lat_0[-1, :] = lat_0[0, :]

        # if 360 scan, upsample to super res
        if lon.shape[0] < 720:
            # upsample via kron
            z = np.kron(z, np.ones((2, 1)))

        # if first sweep, create new dim, otherwise concat them up 
        if sweep == 0:
            z_stack = copy.deepcopy(z[np.newaxis, :, :])
        else:
            z_stack = np.concatenate([z_stack,z[np.newaxis, :, :]])

    # now that the stack is made, take max across vertical 
    compz = z_stack.max(axis=0)

    # since we are using the whole volume scan, report mean time 
    dtime = to_datetime(num2date(radar.time['data'],
                        radar.time['units']).astype(str))
    dtime = dtime.mean()

    # return dict, because this is was pyart does with lots of things 
    comp_dict = {}
    comp_dict['longitude'] = {'data': lon_0, 'units': 'degrees',
                              'info': 'reordered longitude grid, [az,range]'}
    comp_dict['latitude'] = {'data': lat_0, 'units': 'degrees',
                             'info': 'reordered latitude grid, [az,range]'}
    comp_dict['composite_reflectivity'] = {
        'data': compz, 'units': 'dBZ',
        'info': 'composite refelctivity computed from calculating the max radar value in each radar gate vertically after reordering'}
    comp_dict['time'] = {'data': dtime, 'units': 'timestamp',
                         'info': 'mean time of all scans'}
    return comp_dict
