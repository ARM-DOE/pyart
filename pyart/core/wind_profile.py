"""
pyart.core.wind_profile
=======================

Storage of wind profiles.

.. autosummary::
    :toctree: generated/

    HorizontalWindProfile

"""


import numpy as np


class HorizontalWindProfile(object):
    """
    Horizontal wind profile.

    Parameters
    ----------
    height : array-like, 1D
        Heights in meters above sea level at which horizontal winds were
        sampled.
    speed : array-like, 1D
        Horizontal wind speed in meters per second at each height sampled.
    direction : array-like, 1D
        Horizontal wind direction in degrees at each height sampled.

    Other Parameters
    ----------------
    latitude : array-like, 1D, optional
        Latitude in degrees north at each height sampled.
    longitude : array-like, 1D, optional
        Longitude in degrees east at each height sampled.

    Attributes
    ----------
    height : array, 1D
        Heights in meters above sea level at which horizontal winds were
        sampled.
    speed : array, 1D
        Horizontal wind speed in meters per second at each height.
    direction : array, 1D
        Horizontal wind direction in degrees at each height.
    u_wind : array, 1D
        U component of horizontal winds in meters per second at each height.
    v_wind : array, 1D
        V component of horizontal winds in meters per second at each height.

    """

    def __init__(self, height, speed, direction, latitude=None,
                 longitude=None):
        """ initialize """
        if len(height) != len(speed) or len(height) != len(direction):
            raise ValueError("Wind parameters must have the same length.")
        self.height = np.asanyarray(height)
        self.speed = np.asanyarray(speed)
        self.direction = np.asanyarray(direction)
        self._parse_location_data(latitude, longitude)

    @classmethod
    def from_u_and_v(cls, height, u_wind, v_wind):
        """
        Create a HorizontalWindProfile instance from U and V components.

        Parameters
        ----------
        height : array-like, 1D
            Heights in meters above sea level at which horizontal winds were
            sampled.
        u_wind : array-like, 1D
            U component of horizontal wind speed in meters per second.
        v_wind : array-like, 1D
            V component of horizontal wind speed in meters per second.

        """
        u_wind = np.asanyarray(u_wind)
        v_wind = np.asanyarray(v_wind)
        speed = np.sqrt(u_wind*u_wind + v_wind*v_wind)
        direction = np.rad2deg(np.arctan2(-u_wind, -v_wind))
        direction[direction < 0] += 360
        return cls(height, speed, direction)

    @property
    def u_wind(self):
        """ U component of horizontal wind in meters per second. """
        # U = -sin(direction) * speed
        u_wind = -np.sin(np.deg2rad(self.direction)) * self.speed
        return u_wind

    @property
    def v_wind(self):
        """ V component of horizontal wind in meters per second. """
        # V = -cos(direction) * speed
        v_wind = -np.cos(np.deg2rad(self.direction)) * self.speed
        return v_wind

    def _parse_location_data(self, latitude, longitude):
        """ Parse profile location data. """
        if latitude is not None:
            if len(self.height) != len(latitude):
                raise ValueError('Latitude data must have same length.')
            self.latitude = np.asanyarray(latitude)
        if longitude is not None:
            if len(self.height) != len(longitude):
                raise ValueError('Longitude data must have same length.')
            self.longitude = np.asanyarray(longitude)
