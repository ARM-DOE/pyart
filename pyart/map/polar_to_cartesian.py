__author__ = 'wolfensb'

import numpy as np
import matplotlib.pyplot as plt
from scipy import spatial

KE = 4/3. # Constant in the 4/3 earth radius model

def polar_to_cartesian(radar, sweep, field_name, cart_res = 75,
                        max_range = None, mapping = None):  

    """
    Interpolates a PPI or RHI scan in polar coordinates to a regular cartesian
    grid of South-North and West-East coordinates (for PPI) or distance at
    ground and altitude coordinates (for RHI)

    Parameters
    ----------
    radar : Radar
        Radar containing differential phase field.
    sweep : int
        Sweep number to project to cartesian coordinates.
    field_name : str
        Name of the radar field to be interpolated
    cart_res : int, optional
        Resolution (in m.) of the cartesian grid to which polar data is 
        interpolated
    max_range : int, optional
        Maximal allowed range (in m.) from radar for gates to be interpolated 
    mapping : dict, optional
        Dictionnary of mapping indexes (from polar to cartesian), gets returned
        by the function (see below). Can be used as input when interpolating
        sequentially several variables for the same scan, to save significant
        time
        
    Returns
    -------
    coords : tuple of 2 arrays
        2D coordinates of the cartesian grid
    cart_data : 2D array
        Interpolated radar measurements (on the cartesian grid)
    mapping,: dict
        Dictionnary of mapping indexes (from polar to cartesian),which contains
        the indexes mapping the polar grid to the cartesian grid as well as some
        metadata.
    """
    
    # Get data to be interpolated
    pol_data = radar.get_field(sweep,field_name)
            
    is_ppi = radar.sweep_mode['data'][0] == 'ppi'

    if mapping:
        # Check if mapping is usable:
        if is_ppi != mapping['is_ppi']:
            print('Input mapping does not correspond to given scan type, ignoring it')
            mapping = None
        elif mapping['dim_pol'] != pol_data.shape:
            print('Input mapping does not correspond to dimensions of given field' \
                ', ignoring it')
            mapping = None
        else: 
            cart_res = mapping['res']
            max_range = mapping['max_range']        
    # Get distances of radar data
    r = radar.range['data']

    if not max_range:
        max_range = np.max(r)
        
    # Cut data at max_range
    pol_data_cut = pol_data.copy()
    pol_data_cut = pol_data_cut[:,r<max_range]

    # Unmask array
    pol_data_cut = np.ma.masked_values(pol_data, np.nan)

    # One specificity of using the kd-tree is that we need to pad the array
    # with nans at large ranges and angles smaller and larger
    pol_data_cut = np.pad(pol_data_cut,pad_width=((1,1),(0,1)),
                      mode='constant',constant_values=np.nan)
    # Get dim of data
    [N,M] = pol_data_cut.shape

    # Get angles of radar data
    if is_ppi:
        theta = radar.azimuth['data']
    else:
        theta = radar.elevation['data']
    
    # We need to pad theta and r as well    
    theta = np.hstack([theta[0]-0.1,theta,theta[-1]+0.1])

    r = np.hstack([r,r[-1]+0.1])

    r_grid_p, theta_grid_p = np.meshgrid(r,theta)
    
    # Generate regular cartesian grid
    if is_ppi:
        x_vec = np.arange(-max_range-cart_res,max_range+cart_res,cart_res)
        y_vec = np.arange(-max_range-cart_res,max_range+cart_res,cart_res)
    else:
        x_vec = np.arange((max_range-cart_res)*np.cos(np.radians(np.max(theta)))
                          ,max_range+cart_res,cart_res)
        y_vec = np.arange(0,max_range+cart_res,cart_res)
        
    x_grid_c, y_grid_c = np.meshgrid(x_vec, y_vec)

    if is_ppi:
        theta_grid_c= np.degrees(np.arctan2(-x_grid_c, -y_grid_c)+np.pi)
        r_grid_c = (np.sqrt(x_grid_c**2 + y_grid_c**2))
    else:
        theta_grid_c = np.degrees(-(np.arctan2(x_grid_c, y_grid_c)-np.pi/2))
        E = get_earth_radius(radar.latitude['data'])
        r_grid_c = np.sqrt((E*KE*np.sin(np.radians(theta_grid_c)))**2 + \
        2*E*KE*y_grid_c + y_grid_c**2) - E *KE * np.sin(np.radians(theta_grid_c))

    if not mapping:
        # Kd-tree construction and query
        kdtree = spatial.cKDTree(np.vstack((r_grid_p.ravel(),
                                            theta_grid_p.ravel())).T)
        _,mapping_idx = kdtree.query(np.vstack((r_grid_c.ravel(),
                                            theta_grid_c.ravel())).T, k=1)
        
        mapping = {'idx':mapping_idx,'max_range':max_range,'res':cart_res,
                   'is_ppi':is_ppi, 'dim_pol':pol_data.shape}
                   
    cart_data = pol_data_cut.ravel()[mapping['idx']]
    cart_data = np.reshape(cart_data,x_grid_c.shape)
                
    return (x_vec,y_vec),cart_data, mapping

def cartesian_to_polar(radar, x,y,data):
    is_ppi = radar.sweep_mode['data'][0] == 'ppi'
    
    x_grid_c, y_grid_c = np.meshgrid(x,y)
    
    if is_ppi:
        theta_grid_c= np.degrees(np.arctan2(-x_grid_c, -y_grid_c)+np.pi)
        r_grid_c = (np.sqrt(x_grid_c**2 + y_grid_c**2))
    else:
        theta_grid_c = np.degrees(-(np.arctan2(x_grid_c, y_grid_c)-np.pi/2))
        E = _get_earth_radius(radar.latitude['data'])
        r_grid_c = np.sqrt((E*KE*np.sin(np.radians(theta_grid_c)))**2 + \
        2*E*KE*y_grid_c + y_grid_c**2) - E *KE * np.sin(np.radians(theta_grid_c))
    
    # Get distances of radar data
    r = radar.range['data']
    # Get angles of radar data
    if is_ppi:
        theta = radar.azimuth['data']
    else:
        theta = radar.elevation['data']
    
    r_grid_p, theta_grid_p = np.meshgrid(r,theta)
    
    # Kd-tree construction and query
    print('dd')
    kdtree = spatial.cKDTree(np.vstack((r_grid_c.ravel(),
                                        theta_grid_c.ravel())).T)
    print('fu')                          
    _,mapping_idx = kdtree.query(np.vstack((r_grid_p.ravel(),
                                        theta_grid_p.ravel())).T, k=1)
                                    
    cart_data = data.ravel()[mapping_idx]
    cart_data = np.reshape(cart_data,(len(r),len(theta)))
    
    return r,theta, cart_data
            
def get_earth_radius(latitude):
    # Two extreme earth radius
    a=6378.1370*1000
    b=6356.7523*1000

    earth_radius=np.sqrt(((a**2*np.cos(latitude))**2+
        (b**2*np.sin(latitude))**2)/((a*np.cos(latitude))**2+
        (b*np.sin(latitude))**2))
        
    return earth_radius


   
if __name__ == '__main__':
    file_rad ='/ltedata/Payerne_2014/Radar/Proc_data/2014/03/22/MXPol-polar-20140322-112758-RHI-330_0.nc'
    from cosmo_pol.radar import pyart_wrapper
    part = pyart_wrapper.PyradMXPOL(file_rad)

    loc,c1, m = polar_to_cartesian(part,0,'Zh')
    plt.figure()
    plt.contourf(loc[0],loc[1],c1,levels=np.arange(-10,60,1))
    
    loc,c2, m = polar_to_cartesian(part,0,'Zdr')
    plt.figure()
    plt.contourf(loc[0],loc[1],c2,levels=np.arange(-1,3,0.1))

    
    
#    b = np.zeros((114,410)).ravel()
#    n = np.zeros((114,410)).ravel()
#    c = c1.ravel()
#    idx = np.arange(0,len(m['idx']))
#    
#    for i in idx:
#        try:
#            b[m['idx'][i]] += c[i]
#            n[m['idx'][i]] += 1
#            print('d')
#        except:
#            pass
#    b = b/n
#    b = np.reshape(b,((114,410)))