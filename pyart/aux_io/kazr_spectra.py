"""
Utilities for reading ARM KAZR spectra files.

"""

import datetime
import re

import numpy as np
try:
    import xarray as xr
    _XARRAY_AVAILABLE = True
except ImportError:
    _XARRAY_AVAILABLE = False

from ..config import FileMetadata
from ..core.radar_spectra import RadarSpectra
from ..exceptions import MissingOptionalDependency
from ..io.common import _test_arguments


KAZR_SPECTRA_FIELD_NAMES = {'spectra': 'spectra'}

def read_kazr_spectra(filename, field_names=None, additional_metadata=None,
                      file_field_names=False, exclude_fields=None,
                      include_fields=None, **kwargs):
    """
    Read a ARM KAZR spectra netCDF file.

    Parameters
    ----------
    filename : str
        Name of KAZR spectra netCDF file to read data from.
    field_names : dict, optional
        Dictionary mapping field names in the file names to radar_spectra field names.
        Unlike other read functions, fields not in this dictionary or having a
        value of None are still included in the radar_spectra.fields dictionary, to
        exclude them use the `exclude_fields` parameter. Fields which are
        mapped by this dictionary will be renamed from key to value.
    additional_metadata : dict of dicts, optional
        This parameter is not used, it is included for uniformity.
    file_field_names : bool, optional
        True to force the use of the field names from the file in which
        case the `field_names` parameter is ignored. False will use to
        `field_names` parameter to rename fields.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar spectra object. This is applied
        after the `file_field_names` and `field_names` parameters. Set
        to None to include all fields specified by include_fields.
    include_fields : list or None, optional
        List of fields to include from the radar spectra object. This is applied
        after the `file_field_names` and `field_names` parameters. Set
        to None to include all fields not specified by exclude_fields.
    delay_field_loading : bool
        True to delay loading of field data from the file until the 'data'
        key in a particular field dictionary is accessed. In this case
        the field attribute of the returned RadarSpectra object will contain
        LazyLoadDict objects not dict objects. Delayed field loading will not
        provide any speedup in file where the number of gates vary between
        rays (ngates_vary=True) and is not recommended.

    Returns
    -------
    radar_spectra : RadarSpectra
        RadarSpectra object.

    """
    if not _XARRAY_AVAILABLE:
        raise MissingOptionalDependency(
            "Xarray is required to read KAZR spectra files "
            "but is not installed!")
    # test for non empty kwargs
    _test_arguments(kwargs)

     # create metadata retrieval object
    if field_names is None:
        field_names = KAZR_SPECTRA_FIELD_NAMES
    filemetadata = FileMetadata('kazr_spectra', field_names,
                                additional_metadata,
                                file_field_names, exclude_fields)
    # read the data
    xrobj = xr.open_dataset(filename)
    # retrieve the correct time format
    # time format changed in later KAZR files
    try:
        times = [(datetime.timedelta(seconds=sec) +
                  datetime.datetime.utcfromtimestamp(
                      xrobj.base_time.values.astype("O")/1e9))
                 for sec in xrobj.time_offset.values]
    except TypeError:
        ts = [(i - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
              for i in xrobj.time_offset.values]
        bs = (
            xrobj.base_time.values - np.datetime64(
                '1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
        times = [(datetime.timedelta(seconds=sec) +
                  datetime.datetime.utcfromtimestamp(bs.astype("O")/1e9))
                 for sec in ts]
    times = np.array(times)
    time_dict = filemetadata('time')
    base_time = ('seconds since ' + datetime.datetime.strftime(
        times[0], '%Y-%m-%DT%H:%M:%SZ'))
    time_array = [(x-times[0]).seconds for x in times]
    time_dict['units'] = base_time
    time = xr.DataArray(np.array(time_array), attrs=time_dict, dims='time')

    rng_array = np.squeeze(xrobj.range.values)
    rng_dict = filemetadata('range')
    _range = xr.DataArray(rng_array, attrs=rng_dict, dims='range')

    #metadata = filemetadata('metadata')
    # Retrieve the ARM site info
    if 'nsa' in filename:
        xrobj.attrs['site_id'] = 'nsa'
    elif 'ena' in filename:
        xrobj.attrs['site_id'] = 'ena'
    else:
        xrobj.attrs['site_id'] = 'sgp'

    lat_dict = filemetadata('latitude')
    lon_dict = filemetadata('longitude')
    alt_dict = filemetadata('altitude')
    lat, lon, alt = get_kazr_location(xrobj.attrs['site_id'])
    latitude = xr.DataArray(np.array(lat, dtype='float32'), attrs=lat_dict)
    longitude = xr.DataArray(np.array(lon, dtype='float32'), attrs=lon_dict)
    altitude = xr.DataArray(np.array(alt, dtype='float32'), attrs=alt_dict)

    elevation_dict = filemetadata('elevation')
    elevation = xr.DataArray(np.array(90.0*np.ones(len(time_array)),
                                      dtype='float32'),
                             attrs=elevation_dict, dims='time')

    azimuth_dict = filemetadata('azimuth')
    azimuth = xr.DataArray(np.array(360.0*np.ones(len(time_array)),
                                    dtype='float32'),
                           attrs=azimuth_dict, dims='time')

    fix_agl_dict = filemetadata('fixed_angle')
    fixed_angle = xr.DataArray(np.array(90.0, dtype='float32'),
                               attrs=fix_agl_dict)

    sweep_number_dict = filemetadata('sweep_number')
    sweep_number = xr.DataArray(np.array([1], dtype='int32'),
                                attrs=sweep_number_dict)

    sweep_mode_dict = filemetadata('sweep_mode')
    sweep_mode = xr.DataArray(np.array('spectra'),
                              attrs=sweep_mode_dict)

    sweep_start_ray_index_dict = filemetadata('sweep_start_ray_index')
    sweep_start_ray_index = xr.DataArray(np.array([0], dtype='int32'),
                                         attrs=sweep_start_ray_index_dict)
    sweep_end_ray_index_dict = filemetadata('sweep_end_ray_index')
    sweep_end_ray_index = xr.DataArray(np.array([len(times)-1], dtype='int32'),
                                       attrs=sweep_end_ray_index_dict)
    spectra_array = [_get_spectra(xrobj, x) for x in xrobj.time.values]
    spectras = np.stack(spectra_array)

    c = 299792458

    npulses_max = xrobj.speclength.values
    wavelength = c/float(re.findall(r"[-+]?\d*\.\d+|\d+",
                                    xrobj.radar_operating_frequency)[0])*1e-9
    scan_type = 'vpt'

    metadata = xrobj.attrs
    metadata['wavelength'] = wavelength
    velocity_dict = xrobj.velocity_bins.attrs
    bins = xrobj.velocity_bins.values
    velocity_bins = xr.DataArray(bins, attrs=velocity_dict,
                                 dims='npulses_max')
    xrobj.close()

    return RadarSpectra(time, _range, spectras, metadata, scan_type,
                        latitude, longitude, altitude,
                        sweep_number, sweep_mode, fixed_angle,
                        sweep_start_ray_index, sweep_end_ray_index,
                        azimuth, elevation, npulses_max, velocity_bins)


def _get_spectra(xrobj, time):
    """ Retrieves the spectra values using time and locator mask. """
    the_spectra_locs = xrobj.locator_mask.sel(time=time).values
    the_spectra_loc = np.zeros((len(the_spectra_locs),
                                len(xrobj.speclength.values)))
    the_spectra_locs[np.isnan(the_spectra_locs)] = -9999.0
    for i, locs in enumerate(the_spectra_locs):
        if locs != -9999.0:
            the_spectra_loc[i, :] = xrobj.spectra.values[int(locs), :]
        else:
            the_spectra_loc[i, :] = np.nan*np.ones(len(
                xrobj.speclength.values))
    return the_spectra_loc


def get_kazr_location(arm_site):
    """
    Return the latitude, longitude and altitude of a ARM KAZR Radar.

    Parameters
    ----------
    arm_site : str
        Three letter ARM site acroynm.

    Returns
    -------
    lat, lon, alt : float
        Latitude (in degrees), longitude (in degrees), and altitude
        (in meters above mean sea level) of the KAZR radar.

    """
    loc = ARM_SITES[arm_site.upper()]
    return loc['lat'], loc['lon'], loc['alt']

ARM_SITES = {
    "SGP": {'lat': 36.605, 'lon': -97.485, 'alt': 318},
    "ENA": {'lat': 39.0916, 'lon': -28.0257, 'alt': 30},
    "NSA": {'lat': 71.323, 'lon': -156.609, 'alt': 8}, }
