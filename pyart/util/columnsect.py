"""
Function for extracting the radar column above a target
given position in latitude, longitude

"""

import numpy as np
import pandas as pd
import xarray as xr

from ..core.transforms import antenna_vectors_to_cartesian
from .datetime_utils import datetime_from_radar


def column_vertical_profile(
    radar, latitude, longitude, azimuth_spread=3, spatial_spread=3
):
    """
    Given the location (in latitude, longitude) of a target, return the rays
    that correspond to radar column above the target, allowing for user
    defined range of azimuths and range gates to be included within this
    extraction.

    Parameters
    ----------
    radar : pyart.core.Radar Object
        Py-ART Radar Object from which distance to the target, along
        with gates above the target, will be calculated.
    latitude : float, [degrees]
        Latitude, in degrees North, of the target.
    longitude : float, [degrees]
        Longitude, in degrees East, of the target.
    azimuth_spread : int
        Number of azimuth angles to include within extraction list
    spatial_range : int
        Number of range gates to include within the extraction

    Function Calls
    --------------
    sphere_distance
    for_azimuth
    get_sweep_rays
    subset_fields
    assemble_column

    Returns
    -------
    column : xarray
        Xarray Dataset containing the radar column above the target for
        the various fields within the radar object.

    References
    ----------
    Murphy, A. M., A. Ryzhkov, and P. Zhang, 2020: Columnar Vertical
    Profile (CVP) Methodology for Validating Polarimetric Radar Retrievals
    in Ice Using In Situ Aircraft Measurements. J. Atmos. Oceanic Technol.,
    37, 1623–1642, https://doi.org/10.1175/JTECH-D-20-0011.1.

    Bukovčić, P., A. Ryzhkov, and D. Zrnić, 2020: Polarimetric Relations for
    Snow Estimation—Radar Verification. J. Appl. Meteor. Climatol.,
    59, 991–1009, https://doi.org/10.1175/JAMC-D-19-0140.1.
    """

    # Define the spatial range to use within extraction
    spatial_range = radar.range["meters_between_gates"] * spatial_spread

    # Define a dictionary structure to contain the extracted features
    total_moment = {key: [] for key in radar.fields.keys()}
    total_moment.update({"height": [], "time_offset": []})

    # Define the start of the radar volume
    base_time = np.datetime64(datetime_from_radar(radar).isoformat(), "ns")

    # call the sphere_distance function
    dis = sphere_distance(
        radar.latitude["data"][0], latitude, radar.longitude["data"][0], longitude
    )
    # calculate forward azimuth angle
    forazi = for_azimuth(
        radar.latitude["data"][0], latitude, radar.longitude["data"][0], longitude
    )

    # Iterate through radar sweeps, extract desired section
    for sweep in radar.iter_slice():
        moment = {key: [] for key in radar.fields.keys()}
        zgates = []
        gate_time = []

        # call the new sweep rays
        center, spread = get_sweep_rays(
            radar.azimuth["data"][sweep], forazi, azimuth_spread=azimuth_spread
        )

        # add the start indice of each ray
        center = [x + sweep.start for x in center]
        spread = [x + sweep.start for x in spread]
        # Correct the spread indices to remove the centerline
        spread = [x for x in spread if x not in center]

        # For the ray(s) directly over the target, extract and average fields
        for ray in center:
            # Convert gates from antenna or cartesian coordinates
            (rhi_x, rhi_y, rhi_z) = antenna_vectors_to_cartesian(
                radar.range["data"],
                radar.azimuth["data"][ray],
                radar.elevation["data"][ray],
                edges=False,
            )
            # Calculate distance to target
            rhidis = np.sqrt((rhi_x**2) + (rhi_y**2)) * np.sign(rhi_z)
            # Calculate target gate
            tar_gate = np.nonzero(np.abs(rhidis[0, :] - dis) < spatial_range)[
                0
            ].tolist()
            # Subset the radar fields for the target locations
            subset = subset_fields(radar, ray, tar_gate)
            # Add back to the total dictionary
            moment = {key: moment[key] + subset[key] for key in moment}
            # Add radar elevation to height gates
            # to define height as center of each gate above sea level
            zgates.append(np.ma.mean(rhi_z[0, tar_gate] + radar.altitude["data"][0]))
            # Determine the time for the individual gates
            gate_time.append(radar.time["data"][ray])

        # Convert to Cartesian Coordinates
        # Determine the center of each gate for the subsetted rays.
        for ray in spread:
            (rhi_x, rhi_y, rhi_z) = antenna_vectors_to_cartesian(
                radar.range["data"],
                radar.azimuth["data"][ray],
                radar.elevation["data"][ray],
                edges=False,
            )
            # Calculate distance to target
            rhidis = np.sqrt((rhi_x**2) + (rhi_y**2)) * np.sign(rhi_z)
            # Calculate target gate
            tar_gate = np.nonzero(np.abs(rhidis[0, :] - dis) < spatial_range)[
                0
            ].tolist()
            # Subset the radar fields for the target locations
            subset = subset_fields(radar, ray, tar_gate)
            # Add back to the sweep dictionary
            moment = {key: moment[key] + subset[key] for key in moment}
            # Add radar elevation to height gates
            # to define height as center of each gate above sea level
            zgates.append(np.ma.mean(rhi_z[0, tar_gate] + radar.altitude["data"][0]))
            # Determine the time for the individual gates
            gate_time.append(radar.time["data"][ray])

        # Average all azimuth moments into a single value for the sweep
        for key in total_moment:
            if key == "height":
                total_moment[key].append(np.ma.mean(np.ma.masked_invalid(zgates)))
            elif key == "time_offset":
                total_moment[key].append(np.round(np.ma.mean(np.array(gate_time)), 4))
            else:
                total_moment[key].append(
                    np.round(np.ma.mean(np.ma.masked_invalid(moment[key])), 4)
                )

    # Add the base time for the radar
    total_moment.update({"base_time": base_time})

    # Convert to xarray
    return assemble_column(radar, total_moment, forazi, dis, latitude, longitude)


def sphere_distance(radar_latitude, target_latitude, radar_longitude, target_longitude):
    """
    Calculated of the great circle distance between radar and target

    Assumptions
    -----------
    Radius of the Earth = 6371 km / 6371000 meters
    Distance is calculated for a smooth sphere
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
    distance : float, [meters]
        Great-Circle Distance between radar and target in meters
    """
    # check if latitude, longitudes are valid
    check_latitude(radar_latitude)
    check_latitude(target_latitude)
    check_longitude(radar_longitude)
    check_longitude(target_longitude)

    # convert latitudeitude/longitudegitudes to radians
    radar_latitude = radar_latitude * (np.pi / 180.0)
    target_latitude = target_latitude * (np.pi / 180.0)
    radar_longitude = radar_longitude * (np.pi / 180.0)
    target_longitude = target_longitude * (np.pi / 180.0)

    # difference in latitudeitude
    d_latitude = target_latitude - radar_latitude
    # difference in longitudegitude
    d_longitude = target_longitude - radar_longitude

    # Haversine formula
    numerator = (np.sin(d_latitude / 2.0) ** 2.0) + np.cos(radar_latitude) * np.cos(
        target_latitude
    ) * (np.sin(d_longitude / 2.0) ** 2.0)
    distance = 2 * 6371000 * np.arcsin(np.sqrt(numerator))

    # return the output
    return distance


def for_azimuth(radar_latitude, target_latitude, radar_longitude, target_longitude):
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
    radar_latitude = radar_latitude * (np.pi / 180.0)
    target_latitude = target_latitude * (np.pi / 180.0)
    radar_longitude = radar_longitude * (np.pi / 180.0)
    target_longitude = target_longitude * (np.pi / 180.0)

    # Differnce in longitudegitudes
    d_longitude = target_longitude - radar_longitude

    # Determine x,y coordinates for arc tangent function
    corr_y = np.sin(d_longitude) * np.cos(target_latitude)
    corr_x = (np.cos(radar_latitude) * np.sin(target_latitude)) - (
        np.sin(radar_latitude) * (np.cos(target_latitude) * np.cos(d_longitude))
    )

    # Determine forward azimuth angle
    azimuth = np.arctan2(corr_y, corr_x) * (180.0 / np.pi)

    # Return the output as a function of 0-360 degrees
    if azimuth < 0:
        azimuth += 360.0

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
    column : xarray DataSet
        Xarray Dataset containing the radar column above the target for
        the various fields within the radar object.
    """

    # Make sure latitude, longitudes are valid
    check_latitude(latitude)
    check_longitude(longitude)

    # initiate a dictionary to hold the gates above the radar.
    zgate = []

    # initiate a diciontary to hold the moment data.
    moment = {key: [] for key in radar.fields.keys()}

    # call the sphere_distance function
    dis = sphere_distance(
        radar.latitude["data"][0], latitude, radar.longitude["data"][0], longitude
    )

    # call the for_azimuth function
    azim = for_azimuth(
        radar.latitude["data"][0], latitude, radar.longitude["data"][0], longitude
    )

    # call the get_column_ray function
    ray = get_column_rays(radar, azim)

    # Determine the center of each gate for the subsetted rays.
    (rhi_x, rhi_y, rhi_z) = antenna_vectors_to_cartesian(
        radar.range["data"],
        radar.azimuth["data"][ray],
        radar.elevation["data"][ray],
        edges=False,
    )
    # Calculate distance from the x,y coordinates to target
    rhidis = np.sqrt((rhi_x**2) + (rhi_y**2)) * np.sign(rhi_z)
    for i in range(len(ray)):
        tar_gate = np.argmin(abs(rhidis[i, 1:] - (dis)))
        for key in moment:
            if radar.fields[key]["data"][ray[i], tar_gate] is np.ma.masked:
                moment[key].append(np.nan)
            else:
                moment[key].append(radar.fields[key]["data"][ray[i], tar_gate])
        # Add radar elevation to height gates
        # to define height as center of each gate above sea level
        zgate.append(rhi_z[i, tar_gate] + radar.altitude["data"][0])

    # Determine the time at the center of each ray within the column
    # Define the start of the radar volume as a numpy datetime object for xr
    base_time = np.datetime64(datetime_from_radar(radar).isoformat(), "ns")
    # Convert Py-ART radar object time (time since volume start) to time delta
    # Add to base time to have sequential time within the xr Dataset
    # for easier future merging/work
    combined_time = []
    for i in range(len(ray)):
        delta = pd.to_timedelta(radar.time["data"][ray[i]], unit="s")
        total_time = base_time + delta
        combined_time.append(total_time.to_numpy())

    # Create a blank list to hold the xarray DataArrays
    ds_container = []
    da_meta = [
        "units",
        "standard_name",
        "long_name",
        "valid_max",
        "valid_min",
        "coordinates",
    ]
    # Convert the moment dictionary to xarray DataArray.
    # Apply radar object meta data to DataArray attribute
    for key in moment:
        if key != "height":
            da = xr.DataArray(
                moment[key], coords=dict(height=zgate), name=key, dims=["height"]
            )
            for tag in da_meta:
                if tag in radar.fields[key]:
                    da.attrs[tag] = radar.fields[key][tag]
            # Append to ds container
            ds_container.append(da.to_dataset(name=key))

    # Add additional DataArrays 'base_time' and 'time_offset'
    # if not present within the radar object.
    da_base = xr.DataArray(base_time, name="base_time")
    da_offset = xr.DataArray(
        combined_time, coords=dict(height=zgate), name="time_offset", dims=["height"]
    )
    ds_container.append(da_base.to_dataset(name="base_time"))
    ds_container.append(da_offset.to_dataset(name="time_offset"))

    # Create a xarray DataSet from the DataArrays
    column = xr.merge(ds_container)

    # Assign Attributes for the Height and Times
    height_des = (
        "Height Above Sea Level [in meters] for the Center of Each"
        + " Radar Gate Above the Target Location"
    )
    column.height.attrs.update(
        long_name="Height of Radar Beam",
        units="m",
        standard_name="height",
        description=height_des,
    )

    column.base_time.attrs.update(long_name="UTC Reference Time", units="seconds")

    time_long = "Time in Seconds Since Volume Start"
    time_des = (
        "Time in Seconds Since Volume Start that Cooresponds"
        + " to the Center of Each Height Gate"
        + " Above the Target Location"
    )
    column.time_offset.attrs.update(
        long_name=time_long, units="seconds", description=time_des
    )

    # Assign Global Attributes to the DataSet
    column.attrs["distance_from_radar"] = str(np.around(dis / 1000.0, 3)) + " km"
    column.attrs["azimuth"] = str(np.around(azim, 3)) + " degrees"
    column.attrs["latitude_of_location"] = str(latitude) + " degrees"
    column.attrs["longitude_of_location"] = str(longitude) + " degrees"
    return column


def get_column_rays(radar, azimuth):
    """
    Given the location (in latitude,longitude) of a target, return the rays
    that correspond to radar column above the target.

    Parameters
    ----------
    radar : Radar Object
        Py-ART Radar Object from which distance to the target, along
        with gates above the target, will be calculated.
    azimuth : float,int
        forward azimuth angle from radar to target in degrees.

    Returns
    -------
    nrays : List
        radar ray indices that correspond to the column above a
        target location.
    """
    if isinstance(azimuth, int) or isinstance(azimuth, float) is True:
        if (azimuth <= 0) or (azimuth >= 360):
            raise ValueError("azimuth not valid (not between 0-360 degrees)")
    else:
        raise TypeError(
            "radar azimuth type not valid."
            " Please convert input to be an int or float."
        )
    # define a list to hold the valid rays
    rays = []
    # check to see which radar scan
    if radar.scan_type == "rhi":
        for i in range(radar.sweep_number["data"].shape[0]):
            nstart = radar.sweep_start_ray_index["data"][i]
            nstop = radar.sweep_end_ray_index["data"][i]
            counter = 0
            for j in range(nstart, nstop):
                if abs(radar.azimuth["data"][nstart + counter] - azimuth) < 1:
                    rays.append(nstart + counter)
                counter += 1
    else:
        # taken from pyart.graph.RadarDisplay.get_azimuth_rhi_data_x_y_z
        for sweep in radar.iter_slice():
            sweep_azi = radar.azimuth["data"][sweep]
            nray = np.argmin(np.abs(sweep_azi - azimuth))
            rays.append(nray + sweep.start)
    # make sure rays were found
    if len(rays) == 0:
        raise ValueError("No rays were found between azimuth and target")

    return rays


def get_sweep_rays(sweep_azi, azimuth, azimuth_spread=0):
    """
    Extract the specific rays for a given azimuth from a radar sweep

    Azimuth spread determines the +/- degrees azimuth to include within
    the extraction by multipling the azimuth resolution by input value.

    Parameters
    ----------
    radar_sweep : pyart.core.radar object
        Radar Sweep from which the rays are extracted from
    azimuth : float [degrees]
        Forward Azimuth Angle from Radar to Target in Degreees
    Azimuth_Spread : int
        Number of azimuth angles to include within extraction list

    Returns
    -------
    center_rays : list [integers]
        List of integers cooresponding to ray indices within the azimuth
        directly over the target
    spread_rays : list [integers]
        List of integers cooresponding to ray indices within the spread
        of azimuths emcompassing the target

    """
    # determine resolution of azimuth angles
    resolution = np.round((sweep_azi[1] - sweep_azi[0]), 3)

    centerline = np.nonzero(np.abs(sweep_azi - azimuth) < 0.5)[0].tolist()
    spread = np.nonzero(np.abs(sweep_azi - azimuth) < (resolution * azimuth_spread))[
        0
    ].tolist()

    return centerline, spread


def subset_fields(radar, ray, target_gates):
    """
    Parameter
    ---------
    radar : pyart.core.radar object
        Radar Sweep from which fields are extracted from the target locations
    target_gates : list
        List containing indices for the gates of interest

    Returns
    -------
    fields : dict
        dictionary containing averaged subset fields for target location

    """
    # initiate a diciontary to hold the moment data.
    moment = {key: [] for key in radar.fields.keys()}

    # Iterate over input rays and average according to input method
    # future - allow users to input weights for spatial averaging
    for key in moment:
        if key != "height":
            if np.ma.all(radar.fields[key]["data"][ray, target_gates]) is np.ma.masked:
                moment[key].append(np.nan)
            else:
                moment[key].append(
                    np.ma.mean(radar.fields[key]["data"][ray, target_gates])
                )

    return moment


def assemble_column(radar, total_moment, azimuth, distance, latitude, longitude):
    """
    With a dictionary containing the extracted fields from a radar sweep,
    assemble individual gates and fields into an xarray DataSet

    Parameters
    ----------
    total_moment : dict
        Dictionary containing the extracted fields from the radar object.
        File requires at least the height of the individual gates and
        the start time of the volumetric scan.
    azimuth : float, [degrees]
        azimuth angle from the radar where
        target is located within the scan.
        output is in degrees.
    distance : float, [meters]
        Great-Circle Distance between radar and target in meters
    latitude : float, [degrees]
        Latitude of the target in degrees
    longitude : float, [degrees]
        Longitude of the target in degrees

    Returns
    -------
    column : xarray DataSet
        Xarray Dataset containing the radar column above the target for
        the various fields within the radar object.

    """
    # Create a blank list to hold the xarray DataArrays
    ds_container = []
    da_meta = [
        "units",
        "standard_name",
        "long_name",
        "valid_max",
        "valid_min",
        "coordinates",
    ]
    # Skip these fields and apply meta data after creation
    # of the Xarray DataSet
    skip = ["height", "base_time"]
    # Convert the moment dictionary to xarray DataArray.
    # Apply radar object meta data to DataArray attribute
    for key in total_moment:
        if key not in skip:
            # Convert Masked Array elements to NaNs for Xarray
            total_moment[key] = [
                np.nan if x is np.ma.masked else x for x in total_moment[key]
            ]
            # Convert to Xarray DataArray, set derived height as
            # dimension/coordinates
            da = xr.DataArray(
                total_moment[key],
                coords=dict(height=total_moment["height"]),
                name=key,
                dims=["height"],
            )
            # Add meta data for the radar fields
            if key != "time_offset":
                for tag in da_meta:
                    if tag in radar.fields[key]:
                        da.attrs[tag] = radar.fields[key][tag]
            # Append to ds container
            ds_container.append(da.to_dataset(name=key))

    # Create a xarray DataSet from the DataArrays
    column = xr.merge(ds_container)

    # Add the scan times back into the merged column
    column["base_time"] = total_moment["base_time"]
    column.base_time.attrs.update(
        long_name=(
            "Start time of individual radar scan volumes "
            + " from which column are extracted "
        ),
        units="UTC Time",
    )

    # Assign Attributes for the Height and Times
    height_des = (
        "Height Above Sea Level [in meters] for the Center of Each"
        + " Radar Gate Above the Target Location"
    )
    column.height.attrs.update(
        long_name="Height of Radar Beam",
        units="m",
        standard_name="height",
        description=height_des,
    )

    time_long = "Time in Seconds Since Volume Start to the Center of Each Gate"
    time_des = (
        "Time in Seconds Since Volume Start (i.e. base_time) that Cooresponds"
        + " to the Center of Each Height Gate"
        + " Above the Target Location"
    )
    column.time_offset.attrs.update(
        long_name=time_long, units="seconds", description=time_des
    )

    # Add latitude, longitude as variable in extracted column
    column["latitude"] = latitude
    column.latitude.attrs.update(
        long_name="Latitude of Location Column is Extracted Above", units="deg"
    )
    column["longitude"] = longitude
    column.longitude.attrs.update(
        long_name="Longitude of Location Column is Extracted Above", units="deg"
    )

    # Assign Global Attributes to the DataSet
    column.attrs["distance_from_radar"] = str(np.around(distance / 1000.0, 3)) + " km"
    column.attrs["azimuth"] = str(np.around(azimuth, 3)) + " degrees"
    column.attrs["latitude_of_location"] = str(latitude) + " degrees"
    column.attrs["longitude_of_location"] = str(longitude) + " degrees"

    # Drop duplicated heights, keep latest value
    column = column.drop_duplicates(dim="height", keep="last")

    return column


def check_latitude(latitude):
    """
    Function to check if input latitude is valid for type and value.

    Parameters
    ----------
    latitude : int, float
        Latitude of a location that should be between 90S and 90N

    """
    if (
        isinstance(latitude, int)
        or isinstance(latitude, float)
        or isinstance(latitude, np.floating)
    ) is True:
        if (latitude <= -90) or (latitude >= 90):
            raise ValueError(
                "Latitude not between -90 and 90 degrees, need to "
                "convert to values between -90 and 90"
            )
    else:
        raise TypeError(
            "Latitude type not valid, need to convert input to be an int or float"
        )


def check_longitude(longitude):
    """
    Function to check if input latitude is valid for type and value.

    Parameters
    ----------
    longitude : int, float
        Longitude of a location taht should be between 180W and 180E
    """
    if (
        isinstance(longitude, int)
        or isinstance(longitude, float)
        or isinstance(longitude, np.floating)
    ) is True:
        if (longitude <= -180) or (longitude >= 180):
            raise ValueError(
                "Longitude not valid between -180 and 180"
                " degrees, need to convert to values between"
                "  -180 and 180"
            )
    else:
        raise TypeError(
            "Longitude type not valid, need to convert input to be an int or float"
        )
