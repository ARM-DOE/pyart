"""
Function for extracting the radar column above a target
given position in lat, lon

"""

import numpy as np

from ..core.transforms import antenna_vectors_to_cartesian


def sphere_distance(rad_lat, tar_lat, rad_lon, tar_lon):

    """
    Calculation of the great circle distance between radar and target

    Assumptions
    -----------
    Radius of the Earth = 6371 km / 6371000 meters
    Distance is calculated for a smooth sphere
    Radar and Target are at the same altitude (need to check)

    Parameter
    ---------
    rad_lat : float, [degrees]
        latitude of the radar in degrees
    tar_lat : float, [degrees]
        latidude of the target in degrees
    rad_lon : float, [degrees]
        longitude of the radar in degrees
    tar_lon : float, [degrees]
        longitude of the target in degress

    Returns
    -------
    distance : float, [meters]
        Great-Circle Distance between radar and target in meters
    """
    # check if lat, lons are valid
    check_lat(rad_lat)
    check_lat(tar_lat)
    check_lon(rad_lon)
    check_lon(tar_lon)
    # convert latitude/longitudes to radians
    rad_lat = rad_lat * (np.pi/180.)
    tar_lat = tar_lat * (np.pi/180.)
    rad_lon = rad_lon * (np.pi/180.)
    tar_lon = tar_lon * (np.pi/180.)
    # difference in latitude
    d_lat = (tar_lat - rad_lat)
    # difference in longitude
    d_lon = (tar_lon - rad_lon)
    # Haversine formula
    numerator = ((np.sin(d_lat/2.0)**2.0)
                 + np.cos(rad_lat) * np.cos(tar_lat)
                 * (np.sin(d_lon/2.0)**2.0))
    distance = 2 * 6371000 * np.arcsin(np.sqrt(numerator))

    # return the output
    return distance


def for_azimuth(rad_lat, tar_lat, rad_lon, tar_lon):

    """
    Calculation of inital bearing along a great-circle arc
    Known as Forward Azimuth Angle.

    Assumptions
    -----------
    Radius of the Earth = 6371 km / 6371000 meters
    Distance is calculated for a smooth sphere
    Radar and Target are at the same altitude (need to check)

    Parameters
    ----------
    rad_lat : float, [degrees]
        latitude of the radar in degrees
    tar_lat : float, [degrees]
        latidude of the target in degrees
    rad_lon : float, [degrees]
        longitude of the radar in degrees
    tar_lon : float, [degrees]
        longitude of the target in degress

    Returns
    -------
    azimuth : float, [degrees]
        azimuth angle from the radar where
        target is located within the scan.
        output is in degrees.
    """
    # check the input lat/lon
    check_lat(rad_lat)
    check_lat(tar_lat)
    check_lon(rad_lon)
    check_lon(tar_lon)
    # convert latitude/longitudes to radians
    rad_lat = rad_lat * (np.pi/180.)
    tar_lat = tar_lat * (np.pi/180.)
    rad_lon = rad_lon * (np.pi/180.)
    tar_lon = tar_lon * (np.pi/180.)
    # Differnce in longitudes
    d_lon = tar_lon - rad_lon
    # Determine x,y coordinates for arc tangent function
    corr_y = np.sin(d_lon) * np.cos(tar_lat)
    corr_x = ((np.cos(rad_lat) * np.sin(tar_lat))
              - (np.sin(rad_lat) * np.cos(tar_lat) * np.cos(d_lon)))
    # Determine forward azimuth angle
    azimuth = np.arctan2(corr_y, corr_x) * (180./np.pi)

    # Return the output as a function of 0-360 degrees
    if azimuth < 0:
        azimuth += 360.

    return azimuth


def get_fields_latlon(radar, lat, lon):

    """
    Given the location (in lat,lon) of a target, extract the radar column
    above that point for further analysis.

    Parameters
    ----------
    radar : Radar Object
        Py-ART Radar Object from which distance to the target, along
        with gates above the target, will be calculated.
    lat : float, [degrees]
        Latitude, in degrees North, of the target.
    lon : float, [degrees]
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
    # Make sure lat, lons are valid
    check_lat(lat)
    check_lon(lon)
    # initiate a dictionary to hold the striped out meat data from each scan
    meta = {'date': [], 'time': [], 'xgate': [], 'ygate': [],
            'zgate': [], 'distance': [], 'azimuth': []
            }
    # initiate a diciontary to hold the moment data.
    moment = {key: [] for key in radar.fields.keys()}
    # call the sphere_distance function
    dis = sphere_distance(radar.latitude['data'][0],
                          lat,
                          radar.longitude['data'][0],
                          lon)
    # call the for_azimuth function
    azim = for_azimuth(radar.latitude['data'][0],
                       lat,
                       radar.longitude['data'][0],
                       lon)
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
    # Calculate distance from the x,y coordinates to target
    rhidis = np.sqrt((rhi_x**2) + (rhi_y**2)) * np.sign(rhi_z)
    for i in range(len(ray)):
        tar_gate = np.argmin(abs(rhidis[i, 1:] - (dis)))
        for key in moment:
            moment[key].append(radar.fields[key]['data'][ray[i], tar_gate])
        # fill the meta dictionary for further analysis
        meta['xgate'].append(rhi_x[i, tar_gate])
        meta['ygate'].append(rhi_y[i, tar_gate])
        meta['zgate'].append(rhi_z[i, tar_gate])
    # Combine the meta and moment data
    column = dict(meta)
    column.update(moment)

    return column


def get_column_rays(radar, azimuth):

    """
    Given the location (in lat,lon) of a target, return the rays that
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
        raise TypeError("radar longitude type not valid")
    # define a list to hold the valid rays
    rays = []
    # check to see which radar scan
    if radar.scan_type == 'rhi':
        for i in range(radar.sweep_number['data'].shape[0]):
            nstart = radar.sweep_start_ray_index['data'][i]
            nstop = radar.sweep_end_ray_index['data'][i]
            counter = 0
            for j in range(nstart, nstop):
                if (abs(radar.azimuth['data'][nstart+counter] -
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


def check_lat(lat):

    """
    Function to check if input latitude is valid for type and value
    """
    if (isinstance(lat, int) or isinstance(lat, float) or
            isinstance(lat, np.floating)) is True:
        if (lat <= -90) or (lat >= 90):
            raise ValueError("latitude not valid")
    else:
        raise TypeError("latitude type not valid")


def check_lon(lon):

    """
    Function to check if input latitude is valid for type and value
    """
    if (isinstance(lon, int) or isinstance(lon, float) or
            isinstance(lon, np.floating)) is True:
        if (lon <= -180) or (lon >= 180):
            raise ValueError("longitude not valid")
    else:
        raise TypeError("longitude type not valid")
