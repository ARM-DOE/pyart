"""
Function for extracting the radar column above a target
given position in latitude, longitude

"""

import numpy as np
import time

from ..core.transforms import antenna_vectors_to_cartesian

# for testing. 
t0 = time.time()

def sphere_distance(radar_latitude, target_latitude, radar_longitude, 
                    target_longitude):
    """
    Calculated of the great circle distance between radar and target

    Assumptions
    -----------
    Radius of the Earth = 6371 km / 6371000 meters
    Distance is calculated for a smooth sphere
    Radar and Target are at the same altitude (need to check)

    Parameter
    ---------
    radar_latitude : float, [degrees]
        latitude of the radar in degrees
    target_latitude : float, [degrees]
        latitude of the target in degrees
    radar_longitude : float, [degrees]
        longitude of the radar in degrees
    target_longitude : float, [degrees]
        longitude of the target in degress

    Returns
    -------
    distance : float, [meters]
        Great-Circle Distance between radar and target in meters
    """
    # check if latitude, longitudes are valid
    check_latitude(radar_latitude)
    check_latitude(target_latitude)
    check_longitude(radar_longitude)
    check_longitude(target_longitude)
    
    # convert latitudeitude/longitudegitudes to radians
    radar_latitude = radar_latitude * (np.pi/180.)
    target_latitude = target_latitude * (np.pi/180.)
    radar_longitude = radar_longitude * (np.pi/180.)
    target_longitude = target_longitude * (np.pi/180.)
    
    # difference in latitudeitude
    d_latitude = (target_latitude - radar_latitude)
    # difference in longitudegitude
    d_longitude = (target_longitude - radar_longitude)
    
    # Haversine formula
    numerator = ((np.sin(d_latitude/2.0)**2.0)
                 + np.cos(radar_latitude) * np.cos(target_latitude)
                 * (np.sin(d_longitude/2.0)**2.0))
    distance = 2 * 6371000 * np.arcsin(np.sqrt(numerator))

    # return the output
    return distance


def for_azimuth(radar_latitude, target_latitude, radar_longitude,
                target_longitude):
    """
    Calculation of inital bearing alongitudeg a great-circle arc
    Known as Forward Azimuth Angle.

    Assumptions
    -----------
    Radius of the Earth = 6371 km / 6371000 meters
    Distance is calculatitudeed for a smooth sphere
    Radar and Target are at the same altitude (need to check)

    Parameters
    ----------
    radar_latitude : float, [degrees]
        latitude of the radar in degrees
    target_latitude : float, [degrees]
        latitude of the target in degrees
    radar_longitude : float, [degrees]
        longitude of the radar in degrees
    target_longitude : float, [degrees]
        longitude of the target in degress

    Returns
    -------
    azimuth : float, [degrees]
        azimuth angle from the radar where
        target is located within the scan.
        output is in degrees.
    """
    # check the input latitude/longitude
    check_latitude(radar_latitude)
    check_latitude(target_latitude)
    check_longitude(radar_longitude)
    check_longitude(target_longitude)
    
    # convert latitudeitude/longitudegitudes to radians
    radar_latitude = radar_latitude * (np.pi/180.)
    target_latitude = target_latitude * (np.pi/180.)
    radar_longitude = radar_longitude * (np.pi/180.)
    target_longitude = target_longitude * (np.pi/180.)
    
    # Differnce in longitudegitudes
    d_longitude = target_longitude - radar_longitude
    
    # Determine x,y coordinates for arc tangent function
    corr_y = np.sin(d_longitude) * np.cos(target_latitude)
    corr_x = ((np.cos(radar_latitude) * np.sin(target_latitude))
              - (np.sin(radar_latitude) * (np.cos(target_latitude) 
              * np.cos(d_longitude))))
    
    # Determine forward azimuth angle
    azimuth = np.arctan2(corr_y, corr_x) * (180./np.pi)

    # Return the output as a function of 0-360 degrees
    if azimuth < 0:
        azimuth += 360.

    return azimuth


def get_field_location(radar, latitude, longitude):

    """
    Given the location (in latitude, longitude) of a target, extract the
    radar column above that point for further analysis.

    Parameters
    ----------
    radar : pyart.core.Radar Object
        Py-ART Radar Object from which distance to the target, along
        with gates above the target, will be calculated.
    latitude : float, [degrees]
        Latitude, in degrees North, of the target.
    longitude : float, [degrees]
        Longitude, in degrees East, of the target.

    Function Calls
    --------------
    sphere_distance
    for_azimuth
    get_column_rays

    Returns
    -------
    column : dict
        Dictionary containing the radar column above the target for
        the various fields within the radar object, as well as,
        the height gates above that point.

    Dictionary keys:
    ----------------
    date : string
        date of the radar scan in YYYY-MM-DD format
    time : string
        beginning time of each scan in HH:MM:SSZ format
    xgate : list
        East-west distance to target in m
    ygate : list
        North-South distance to target in m
    zgate : list
        Height above target in m
    tardis : float
        Distance from radar to target in m
    tarazi : float
        Forward azimuth angle from radar to target in degrees
    moments : list
        Radar fields for the column above the target.
        Dependent on which radar fields are within the
        input radar object.

    """
    # Make sure latitude, longitudes are valid
    check_latitude(latitude)
    check_longitude(longitude)
    
    # initiate a dictionary to hold the striped out meat data from each scan
    meta = {'date': [], 'time': [], 'xgate': [], 'ygate': [],
            'zgate': [], 'distance': [], 'azimuth': []
            }
    
    # initiate a diciontary to hold the moment data.
    moment = {key: [] for key in radar.fields.keys()}
    
    # call the sphere_distance function
    dis = sphere_distance(radar.latitude['data'][0],
                          latitude,
                          radar.longitude['data'][0],
                          longitude)
    
    # call the for_azimuth function
    azim = for_azimuth(radar.latitude['data'][0],
                       latitude,
                       radar.longitude['data'][0],
                       longitude)
    
    # call the get_column_ray function
    ray = get_column_rays(radar, azim)
    
    # add the parameters to the meta dictionary
    meta['date'].append(radar.time['units'].split(' ')[2].split('T')[0])
    meta['time'].append(radar.time['units'].split(' ')[2].split('T')[-1])
    meta['distance'].append(dis / 1000.)
    meta['azimuth'].append(azim)
    
    # Determine the gates for the rays
    (rhi_x,
     rhi_y,
     rhi_z) = antenna_vectors_to_cartesian(radar.range['data'],
                                           radar.azimuth['data'][ray],
                                           radar.elevation['data'][ray],
                                           edges=True)
    # Calculatitudee distance from the x,y coordinates to target
    rhidis = np.sqrt((rhi_x**2) + (rhi_y**2)) * np.sign(rhi_z)
    t1 = time.time()
    for i in range(len(ray)):
        tar_gate = np.argmin(abs(rhidis[i, 1:] - (dis)))
        for key in moment:
            moment[key].append(radar.fields[key]['data'][ray[i], tar_gate])
        # fill the meta dictionary for further analysis
        meta['xgate'].append(rhi_x[i, tar_gate])
        meta['ygate'].append(rhi_y[i, tar_gate])
        meta['zgate'].append(rhi_z[i, tar_gate])
    t2 = time.time()
    print("LOOP TIME: ", t2-t1)
    # Combine the meta and moment data
    column = dict(meta)
    column.update(moment)

    return column


def get_column_rays(radar, azimuth):

    """
    Given the location (in latitude,longitude) of a target, return the rays that
    correspond to radar column above the target.

    Parameters
    ----------
    radar    : Radar Object
               Py-ART Radar Object from which distance to the target, along
               with gates above the target, will be calculated.
    azimuth  : float,int
               forward azimuth angle from radar to target in degrees.

    Returns
    -------
    nrays : List
            radar ray indices that correspond to the column above a
            target location
    """
    if isinstance(azimuth, int) or isinstance(azimuth, float) is True:
        if (azimuth <= 0) or (azimuth >= 360):
            raise ValueError("azimuth not valid (not between 0-360 degrees)")
    else:
        raise TypeError("radar longitudegitude type not valid." 
                        + " Please convert input to be an int or float")
    # define a list to hold the valid rays
    rays = []
    # check to see which radar scan
    if radar.scan_type == 'rhi':
        for i in range(radar.sweep_number['data'].shape[0]):
            nstart = radar.sweep_start_ray_index['data'][i]
            nstop = radar.sweep_end_ray_index['data'][i]
            counter = 0
            for j in range(nstart, nstop):
                if (abs(radar.azimuth['data'][nstart + counter] -
                        azimuth) < 1):
                    rays.append(nstart+counter)
                counter += 1
    else:
        # taken from pyart.graph.RadarDisplay.get_azimuth_rhi_data_x_y_z
        for sweep in radar.iter_slice():
            sweep_azi = radar.azimuth['data'][sweep]
            nray = np.argmin(np.abs(sweep_azi - azimuth))
            rays.append(nray + sweep.start)
    # make sure rays were found
    if len(rays) == 0:
        raise ValueError("No rays were found between azimuth and target")

    return rays


def check_latitude(latitude):

    """
    Function to check if input latitude is valid for type and value
    """
    if (isinstance(latitude, int) or isinstance(latitude, float) or
            isinstance(latitude, np.floating)) is True:
        if (latitude <= -90) or (latitude >= 90):
            raise ValueError("Latitude not between -90 and 90 degrees, need to"
                              + "convert to values between -90 and 90")
    else:
        raise TypeError("Latitude type not valid, need to convert input to be" 
                         + " an int or float")


def check_longitude(longitude):

    """
    Function to check if input latitude is valid for type and value
    """
    if (isinstance(longitude, int) or isinstance(longitude, float) or
            isinstance(longitude, np.floating)) is True:
        if (longitude <= -180) or (longitude >= 180):
            raise ValueError("Longitude not valid between -180 and 180"
                             + " degrees, need to convert to values between"
                             + "  -180 and 180")
    else:
        raise TypeError("Longitude type not valid, need to convert input to"
                         + " be an int or float")

t3 = time.time()
print("TOTAL TIME: ", t3-t0)