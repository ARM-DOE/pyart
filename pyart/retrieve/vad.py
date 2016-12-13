"""
pyart.retrieve.velocity_azimuth_display
=======================================

Retrieval of VADs from a radar object.

Code adapted by Scott Collis,

"Below is some code from Susan Rennie (based off Michelson et al 2000)
from the Centre for Australian Weather and Climate Research (CAWCR)
for calculating VADs. I used her one tilt version and wrote my own
adaptation into Py-ART PyRadar object.. Note I convert to U and V before
averaging.. my concern is that if you have θ rapidly varying between
360° and 0° this will average to a nonsense number" (Collis).

.. autosummary::
    :toctreeL generated/
    :template: dev_template.rst

    Velocity Azimuth Display

"""

import numpy as np
from pyart.core import HorizontalWindProfile


def velocity_azimuth_display(radar, velocity=None,
                             z_start=500, z_end=10500,
                             z_count=101):
    """
    Velocity azimuth display.

    Parameters
    ----------
    radar : Radar
        Radar object used.
    velocity : string
        Velocity field to use for VAD calculation.
        If None, the default velocity field will be used.

    Other Parameters
    ----------------
    z_start : int
        Z location to begin VAD calculation in meters,
        default=0.
    z_end : int
        Z location to end VAD calculation in meters,
        default=10500.
    z_count : int
        Amount of data points used between z_start and z_end.
        Data points are evenly spread out using np.linspace,
        default=101.

    Returns
    -------
    height : array
        Heights in meters above sea level at which horizontal winds were
        sampled.
    speed : array
        Horizontal wind speed in meters per second at each height.
    direction : array
        Horizontal wind direction in degrees at each height.
    u_wind : array
        U-wind mean in meters per second.
    v_wind : array
        V-wind mean in meters per second.

    Reference
    ----------
    Michelson, D. B., Andersson, T., Koistinen, J., Collier, C. G., Riedl, J.,
    Szturc, J., Gjertsen, U., Nielsen, A. and Overgaard, S. (2000) BALTEX Radar
    Data Centre Products and their Methodologies. In SMHI Reports. Meteorology
    and Climatology.Swedish Meteorological and Hydrological Institute, Norrkoping.

    """

    def interval_mean(data, current_z, wanted_z):
        """ Find the mean of data indexed by current_z
            at wanted_z on intervals wanted_z+/- delta
            wanted_z. """
        delta = wanted_z[1] - wanted_z[0]
        pos_lower = [np.argsort((current_z - (
            wanted_z[i] - delta / 2.0))**2)[0]
                     for i in range(len(wanted_z))]
        pos_upper = [np.argsort((current_z - (
            wanted_z[i] + delta / 2.0))**2)[0]
                     for i in range(len(wanted_z))]
        new_values = np.array([data[pos_lower[i]:pos_upper[i]].mean()
                               for i in range(len(pos_upper))])
        return new_values

    def sd_to_uv(speed, direction):
        """ Takes speed and direction to create u_mean and v_mean."""
        return (speed * np.sin(direction), speed * np.cos(direction))

    def vad_algorithm(velocity_field, azimuth, elevation):
        """ Calculates VAD for a scan, returns speed and angle
        outdic = vad_algorithm(velocity_field, azimuth, elevation)
        velocity_field is a 2D array, azimuth is a 1D array,
        elevation is a number.
        All in degrees, m outdic contains speed, angle, variance.
        """
        # Possible name change nbins to ngates.
        nrays, ngates = velocity_field.shape
        nrays2 = nrays / 2
        velocity_count = np.empty((nrays2, ngates, 2))
        velocity_count[:, :, 0] = velocity_field[0:nrays2, :]
        velocity_count[:, :, 1] = velocity_field[nrays2:, :]
        sinaz = np.sin(np.deg2rad(azimuth))
        cosaz = np.cos(np.deg2rad(azimuth))
        sumv = np.ma.sum(velocity_count, 2)
        vals = np.isnan(sumv)
        vals2 = np.vstack((vals, vals))
        # Jonathan, `== False` to `is False` raises an error
        # At the moment `== False` is `now == 0` on the line following
        # this comment. 
        count = np.sum(np.isnan(sumv) == 0, 0)
        aa = count < 8
        vals[:, aa] = 0
        vals2[:, aa] = 0
        count = np.float_(count)
        count[aa] = np.nan
        u_m = np.array([np.nansum(sumv, 0) / (2 * count)])
        count[aa] = 0

        cminusu_mcos = np.zeros((nrays, ngates))
        cminusu_msin = np.zeros((nrays, ngates))
        sincos = np.zeros((nrays, ngates))
        sin2 = np.zeros((nrays, ngates))
        cos2 = np.zeros((nrays, ngates))

        for i in range(ngates):
            cminusu_mcos[:, i] = cosaz * (velocity_field[:, i] - u_m[:, i])
            cminusu_msin[:, i] = sinaz * (velocity_field[:, i] - u_m[:, i])
            sincos[:, i] = sinaz * cosaz
            sin2[:, i] = sinaz**2
            cos2[:, i] = cosaz**2

        cminusu_mcos[vals2] = np.nan
        cminusu_msin[vals2] = np.nan
        sincos[vals2] = np.nan
        sin2[vals2] = np.nan
        cos2[vals2] = np.nan
        sumcminu_mcos = np.nansum(cminusu_mcos, 0)
        sumcminu_msin = np.nansum(cminusu_msin, 0)
        sumsincos = np.nansum(sincos, 0)
        sumsin2 = np.nansum(sin2, 0)
        sumcos2 = np.nansum(cos2, 0)
        a_value = (sumcminu_mcos - (sumsincos*sumcminu_msin / sumsin2)) / (
            sumcos2 - (sumsincos**2) / sumsin2)
        b_value = (sumcminu_msin - a_value*sumsincos) / sumsin2
        speed = np.sqrt(a_value**2 + b_value**2) / np.cos(
            np.deg2rad(elevation))
        angle = np.arctan2(a_value, b_value)

        # Took out code producing variance at the moment, not sure if needed.
        # crv = np.empty((nrays, nbins))
        # for i in range(nbins):
            # crv[:, i] = np.sin(np.radians(azimuth) + angle[i])*speed[i]
        # Need to change Vn name.
        # Vn = velocity_field.copy()
        # Vn[vals2 == True] = np.nan
        # var = np.nansum((crv - Vn)**2, 0) / (sum(np.isnan(Vn) == False) - 2)
        # return {'speed': speed, 'angle': angle, 'variance': var}

        return {'speed': speed, 'angle': angle}

    def vad_calculation(radar):
        """ Appends speed, height, and angle, then uses sd_to_uv
        function and pyart's HorizontalWindProfile function to
        obtain u_wind and v_wind. """
        speed = []
        angle = []
        height = []
        z_gate_data = radar.gate_z['data']
        z_want = np.linspace(z_start, z_end, z_count)

        for i in range(len(radar.sweep_start_ray_index['data'])):
            index_start = radar.sweep_start_ray_index['data'][i]
            index_end = radar.sweep_end_ray_index['data'][i]
            if (index_end - index_start) % 2 == 0:
                print("even, all good")
            else:
                index_end = index_end - 1

            if velocity is None:
                velocity_used = 'velocity'
            else:
                velocity_used = velocity
            velocity_field = radar.fields[
                velocity_used]['data'][index_start:index_end, :]
            azimuth = radar.azimuth['data'][index_start:index_end]
            elevation = radar.fixed_angle['data'][i]

            one_level = vad_algorithm(velocity_field,
                                      azimuth, elevation)
            not_garbage = np.isfinite(one_level['speed'])
            print('max height', z_gate_data[index_start, :][
                np.where(not_garbage)].max(), ' meters')
            speed.append(one_level['speed'][np.where(not_garbage)])
            angle.append(one_level['angle'][np.where(not_garbage)])
            height.append(z_gate_data[index_start, :][np.where(not_garbage)])

        speed_array = np.concatenate(speed)
        angle_array = np.concatenate(angle)
        height_array = np.concatenate(height)
        arg_order = height_array.argsort()
        speed_ordered = speed_array[arg_order]
        height_ordered = height_array[arg_order]
        angle_ordered = angle_array[arg_order]

        u_ordered, v_ordered = sd_to_uv(speed_ordered, angle_ordered)
        u_mean = interval_mean(u_ordered, height_ordered, z_want)
        v_mean = interval_mean(v_ordered, height_ordered, z_want)

        vad = HorizontalWindProfile.from_u_and_v(z_want, u_mean,
                                                 v_mean)
        return vad
    vad = vad_calculation(radar)
    return vad
