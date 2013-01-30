"""
Bindings to the University of Washington 4DD code through the NASA TRMM
RSL code
"""

import datetime

import netCDF4
from numpy import where, isnan, ma

from pyart.io import radar
from pyart.io import py4dd, rsl_utils


class dealiaser:
    """
    Create a first guess field from a sounding and use this to attempt to
    dealias and return dopler velocities.

    Usage:
        radar_field_object = get_dealias_field_sounding(
            radar_obj, sounding_heights, sounding_wind_speeds,
            sounding_wind_direction)

    Parameters
    ----------
    radar_obj: a Py-RART radar object, needs to have nyquist defined in
               the inst_params and have at least reflectivity_horizontal and
               mean_doppler_velocity sounding_heights, sounding_wind_speeds,
    sounding_wind_direction : Numpy arrays in meters, degrees and m/s
    datetime_sounding : A datetime object for the mean time of the profiler

    Note
    ----
    Due to limitations in the C code do not call with numpy arrays over 900
    elements long

    """

    def __init__(self, myradar, sonding_heights, sounding_wind_speeds,
                 sounding_wind_direction, datetime_sounding):
        self.radar = myradar
        self.rsl_radar = rsl_utils.radar_to_rsl(myradar, {
            'reflectivity_horizontal': 'DZ', 'mean_doppler_velocity': 'VR'})

        dayofyear = (datetime.datetime(
            datetime_sounding.year, datetime_sounding.month,
            datetime_sounding.day) -
            datetime.datetime(datetime_sounding.year, 01, 01)).days + 1
        juldate = (datetime_sounding.year -
                   int(("%(d)d" % {'d': datetime_sounding.year})[0:2])
                   * 100) * 1000 + dayofyear
        self.fulljuldate = (juldate * 10000 + datetime_sounding.hour * 100 +
                            datetime_sounding.minute)
        self.sounding_heights = sounding_heights
        self.sounding_wind_speeds = sounding_wind_speeds
        self.sounding_wind_direction = sounding_wind_direction

    def __call__(self, prep=1, low_dbz=-10, filt=1, rsl_badval=131072,
                 fill_value=-9999.0):
        self.my_new_volume, self.sonde_volume = py4dd.dealias_radar_array(
            self.rsl_radar,
            None,
            self.sounding_heights,
            self.sounding_wind_speeds,
            self.sounding_wind_direction,
            self.fulljuldate,
            prep=prep,
            LOWDBZ=low_dbz,
            filt=filt)

        dealiased_data = radar.create_cube_array_lim(
            self.my_new_volume[0],
            self.my_new_volume.contents.h.nsweeps,
            self.my_new_volume.contents.sweeps[0].h.nrays)

        dealiased_data[where(isnan(dealiased_data))] = fill_value
        dealiased_data[where(dealiased_data == rsl_badval)] = fill_value
        # fetch metadata
        meta = self.radar.get_mdv_meta(self.rsl_radar, 'VEL_COR')
        dealiased_fielddict = {'data': ma.masked_equal(
            dealiased_data, -9999.0).reshape(
                dealiased_data.shape[0] * dealiased_data.shape[1],
                dealiased_data.shape[2])}
        dealiased_fielddict.update(meta)
        return dealiased_fielddict


def find_time_in_interp_sonde(interp_sonde, radar_datetime):
    """
    Take a netcdf4 object pointing to an ARM Interpsonde file and get the
    correct time
    """
    sonde_datetimes = netCDF4.num2date(interp_sonde.variables['time'][:],
                                       interp_sonde.variables['time'].units)
    selected = sorted(sonde_datetimes,
                      key=lambda x: abs(x - radar_datetime))[0]
    my_index = list(sonde_datetimes).index(selected)
    return (interp_sonde.variables['height'][:],
            interp_sonde.variables['wspd'][my_index, :],
            interp_sonde.variables['wdir'][my_index, :])
