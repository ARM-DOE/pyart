"""
Routines for reading IPHEx NOXP files.

"""

import datetime

import numpy as np
import netCDF4

from ..config import FileMetadata
from ..io.common import make_time_unit_str  # , _test_arguments
from ..core.radar import Radar

# Only providing updated names for the most common fields.
# Many more fields in the files, and these are read in but not renamed.
NOXP_FIELD_NAMES = {
    # uncorrected reflectivity, horizontal
    'ZH': 'reflectivity',
    # uncorrected differential reflectivity
    'ZDR': 'differential_reflectivity',
    'RHOHV': 'cross_correlation_ratio',
    'KDP': 'differential_phase',
    'KDPM': 'specific_differential_phase',
    'VR': 'velocity',
    'SVR': 'spectrum_width',
}


def read_noxp_iphex_nc(filename, field_names=None, additional_metadata=None,
                       file_field_names=False, exclude_fields=None, **kwargs):
    """
    Read a NOXP IPHEX netCDF file.

    Parameters
    ----------
    filename : str
        Name of the netCDF file to read.
    field_names : dict, optional
        Dictionary mapping netCDF field names to radar field names. If a
        data type found in the file does not appear in this dictionary or has
        a value of None it will not be placed in the radar.fields dictionary.
        A value of None, the default, will use the mapping defined in the
        Py-ART configuration file.
    additional_metadata : dict of dicts, optional
        Dictionary of dictionaries to retrieve metadata from during this read.
        This metadata is not used during any successive file reads unless
        explicitly included.  A value of None, the default, will not
        introduct any addition metadata and the file specific or default
        metadata as specified by the Py-ART configuration file will be used.
    file_field_names : bool, optional
        True to use the netCDF data type names for the field names. If this
        case the field_names parameter is ignored. The field dictionary will
        likely only have a 'data' key, unless the fields are defined in
        `additional_metadata`.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters. Set
        to None to include all fields specified by include_fields.
    include_fields : list or None, optional
        List of fields to include from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters. Set
        to None to include all fields not specified by exclude_fields.

    Returns
    -------
    radar : Radar
        Radar object containing data from netCDF file.

    """
    # TODO before moving to pyart.io
    # * unit test
    # * add default field mapping, etc to default config
    # * auto-detect file type with pyart.io.read function
    # * instrument parameters
    # * add additional checks for HOW attributes
    # * support for other objects (SCAN, XSEC)

    # test for non empty kwargs
    # Getting import error on this function, skipping for now.
    # Only issues a warning anyway.
#    _test_arguments(kwargs)

    # create metadata retrieval object
    if field_names is None:
        field_names = NOXP_FIELD_NAMES
    filemetadata = FileMetadata('cfradial', field_names, additional_metadata,
                                file_field_names, exclude_fields,
                                include_fields)

    # read the data
    ncobj = netCDF4.Dataset(filename)
    ncvars = ncobj.variables

    # One sweep per file
    nsweeps = 1

    # latitude, longitude and altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')
    latitude['data'] = np.array(ncvars['Latitude'])
    longitude['data'] = np.array(ncvars['Longitude'])
    altitude['data'] = np.array(ncvars['Altitude'])

    # metadata
    metadata = filemetadata('metadata')
    metadata['source'] = 'NOAA NSSL'
    metadata['original_container'] = 'noxp_iphex_nc'
    metadata['system'] = ncobj.Radar

    # sweep_start_ray_index, sweep_end_ray_index
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')
    rays_per_sweep = np.shape(ncvars['AZ'][:])
    ssri = np.cumsum(np.append([0], rays_per_sweep[:-1])).astype('int32')
    seri = np.cumsum(rays_per_sweep).astype('int32') - 1
    sweep_start_ray_index['data'] = ssri
    sweep_end_ray_index['data'] = seri

    # sweep_number
    sweep_number = filemetadata('sweep_number')
    sweep_number['data'] = np.arange(nsweeps, dtype='int32')

    # sweep_mode
    sweep_mode = filemetadata('sweep_mode')
    scan_type = getattr(ncobj, 'Scan type').lower()
    if scan_type in ['ppi', 'rhi']:
        sweep_mode['data'] = np.array(nsweeps * ['manual_' + scan_type])
    else:  # Guessing that if not RHI or PPI, then a pointing scan
        sweep_mode['data'] = np.array(nsweeps * ['pointing'])

    # fixed_angle
    fixed_angle = filemetadata('fixed_angle')
    if scan_type == 'rhi':
        sweep_el = ncvars['AZ'][0]
    else:
        sweep_el = ncvars['EL'][0]
    fixed_angle['data'] = np.array([sweep_el], dtype='float32')

    # elevation
    elevation = filemetadata('elevation')
    elevation['data'] = np.array(ncvars['EL'])

    # range
    _range = filemetadata('range')
    rng = 1000.0 * ncvars['Range'][:].T
    _range['data'] = rng[0]
    _range['meters_to_center_of_first_gate'] = np.mean(rng[0][0:2])
    # Gate spacing mostly but not completely constant
    _range['meters_between_gates'] = np.median(np.diff(_range['data']))
    _range['spacing_is_constant'] = 0

    # azimuth
    azimuth = filemetadata('azimuth')
    azimuth['data'] = np.array(ncvars['AZ'])

    # time
    # Basing everything around 6/3/2014 because better than 1/0/0000
    # for datetime module.
    _time = filemetadata('time')
    frac_since_basetime = ncvars['time'][0] - 735762.0
    start_time = \
        datetime.datetime(2014, 6, 3) + datetime.timedelta(frac_since_basetime)
    _time['units'] = make_time_unit_str(start_time)
    _time['data'] = 3600.0 * 24.0 * \
        (ncvars['time'][:]-ncvars['time'][0]).astype('float32')

    # fields
    # nearly all variables w/ dimensions of 'Gate', 'Time' are fields
    keys = [k for k, v in ncvars.items()
            if v.dimensions == ('Gate', 'Time') and
            k not in ['Range', 'Distance', 'Height']]
    fields = {}
    for key in keys:
        field_name = filemetadata.get_field_name(key)
        if field_name is None:
            if exclude_fields is not None and key in exclude_fields:
                continue
            if include_fields is not None and not key in include_fields:
                continue
            field_name = key
        fields[field_name] = _ncvar_to_dict(ncvars[key])

    # instrument_parameters
    instrument_parameters = {}
    prt = filemetadata('prt')
    prt['data'] = 1.0 / float(ncobj.PRF[0:5])
    instrument_parameters['prt'] = prt
    frequency = filemetadata('frequency')
    frequency['frequency'] = float(ncobj.Frequency[0:4]) * 10e9
    instrument_parameters['frequency'] = frequency
    hits = filemetadata('n_samples')
    hits['data'] = ncobj.Hits
    instrument_parameters['hits'] = hits

    # Go get the radar object!
    return Radar(
        _time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)


def _ncvar_to_dict(ncvar):
    """ Convert a NetCDF Dataset variable to a dictionary. """
    # copy all attribute except for scaling parameters
    d = dict((k, getattr(ncvar, k)) for k in ncvar.ncattrs()
             if k not in ['scale_factor', 'add_offset'])
    d['data'] = ncvar[:].T  # Data originally stored as (gate, azimuth)
    if np.isscalar(d['data']):
        # netCDF4 1.1.0+ returns a scalar for 0-dim array, we always want
        # 1-dim+ arrays with a valid shape.
        d['data'] = np.array(d['data'])
        d['data'].shape = (1, )
    return d
