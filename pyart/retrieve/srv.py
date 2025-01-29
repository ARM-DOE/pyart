"""
Calculation of storm-relative velocity from a radar object. Code written by
Edward C. Wolff. Modifications for single-sweep files suggested by Leanne
Blind.

"""

import math

import numpy as np

from ..config import get_field_name


def storm_relative_velocity(
    radar, direction=None, speed=None, field=None, u=None, v=None
):
    """
    This function calculates storm-relative Doppler velocities.

    Parameters
    ----------
    radar: Radar
        Radar object used.
    direction: float or string
        Direction of the storm motion vector (where north equals 0 degrees).
        Accepts a float or a string with the abbreviation of a cardinal or
        ordinal/intercardinal direction (for example: N, SE, etc.). If both
        speed/direction and u/v are specified, speed/direction will be used.
    speed: string
        Speed of the storm motion vector.
        Units should be identical to those in the provided radar
        object. If both speed/direction and u/v are specified, speed/direction
        will be used.
    field: string, optional
        Velocity field to use for storm-relative calculation. A value of None
        will use the default field name as defined in the Py-ART configuration
        file.
    u: float, optional
        U-component of the storm motion
    v: float, optional
        V-component of the storm motion

    Returns
    -------
    sr_data : dict
        Field dictionary containing storm-relative Doppler velocities in the
        same units as original velocities and the specified storm speed.
        Array is stored under the 'data' key.

    """
    # Parse the field parameter
    if field is None:
        field = get_field_name("velocity")

    # Obtain velocity data and copy the array
    sr_data = radar.fields[field]["data"].copy()

    # Specify cardinal directions that can be interpreted
    direction_dict = {
        "N": 0,
        "NE": 45,
        "E": 90,
        "SE": 135,
        "S": 180,
        "SW": 225,
        "W": 270,
        "NW": 315,
    }

    # Set the direction of the storm motion vector
    # When speed and direction are specified
    if direction is not None and speed is not None:
        if isinstance(direction, int) or isinstance(direction, float):
            alpha = direction
        elif isinstance(direction, str):
            if direction in direction_dict.keys():
                alpha = direction_dict[direction]
            else:
                raise ValueError("Direction string must be cardinal/ordinal direction")
        else:
            raise ValueError("Direction must be an integer, float, or string")
    # When u and v are specified
    elif u is not None:
        if v is not None:
            speed = np.sqrt((u**2) + (v**2))
            direction = 90 - np.rad2deg(math.atan2(v / speed, u / speed))
            if direction < 0:
                direction = direction + 360
        else:
            raise ValueError("Must specify both u and v components")
    else:
        raise ValueError("Must specify either speed and direction or u and v")

    # Calculates the storm relative velocities
    # If the radar file contains only one sweep (e.g. some research radars)
    if len(radar.sweep_number["data"]) == 1:
        sweep = 0
        start, end = radar.get_start_end(sweep)
        angle_array = radar.get_azimuth(sweep=sweep)
        ray_array = np.arange(start, end, 1)
        for count, ray in enumerate(ray_array):
            correction = speed * np.cos(np.deg2rad(alpha - angle_array[count]))
            sr_data[ray] = radar.fields[field]["data"][ray] - correction
    # If the radar file contains several sweeps, one volume scan (e.g. NEXRAD)
    else:
        for sweep in radar.sweep_number["data"]:
            start, end = radar.get_start_end(sweep)
            angle_array = radar.get_azimuth(sweep=sweep)
            ray_array = np.arange(start, end + 1, 1)
            for count, ray in enumerate(ray_array):
                correction = speed * np.cos(np.deg2rad(alpha - angle_array[count]))
                sr_data[ray] = radar.fields[field]["data"][ray] - correction

    return sr_data
