"""
Functions for accessing Common Data Model (CDM) NEXRAD Level 2 files.

"""

from datetime import datetime, timedelta
import os

import netCDF4
import numpy as np

from .nexrad_common import get_nexrad_location
from ..config import FileMetadata, get_fillvalue
from ..core.radar import Radar
from .common import make_time_unit_str, _test_arguments


def read_nexrad_cdm(filename, field_names=None, additional_metadata=None,
                    file_field_names=False, exclude_fields=None,
                    include_fields=None, station=None, **kwargs):
    """
    Read a Common Data Model (CDM) NEXRAD Level 2 file.

    Parameters
    ----------
    filename : str
        File name or URL of a Common Data Model (CDM) NEXRAD Level 2 file.
        File of in this format can be created using the NetCDF Java Library
        tools [1]_. A URL of a OPeNDAP file on the UCAR THREDDS Data
        Server [2]_ is also accepted the netCDF4 library has been compiled
        with OPeNDAP support.
    field_names : dict, optional
        Dictionary mapping NEXRAD moments to radar field names. If a
        data type found in the file does not appear in this dictionary or has
        a value of None it will not be placed in the radar.fields dictionary.
        A value of None, the default, will use the mapping defined in the
        metadata configuration file.
    additional_metadata : dict of dicts, optional
        Dictionary of dictionaries to retrieve metadata from during this read.
        This metadata is not used during any successive file reads unless
        explicitly included. A value of None, the default, will not
        introduct any addition metadata and the file specific or default
        metadata as specified by the metadata configuration file will be used.
    file_field_names : bool, optional
        True to use the NEXRAD field names for the field names. If this
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
    station : str
        Four letter ICAO name of the NEXRAD station used to determine the
        location in the returned radar object. This parameter is only
        used when the location is not contained in the file, which occur
        in older NEXRAD files. If the location is not provided in the file
        and this parameter is set to None the station name will be determined
        from the filename.

    Returns
    -------
    radar : Radar
        Radar object containing all moments and sweeps/cuts in the volume.
        Gates not collected are masked in the field data.

    References
    ----------
    .. [1] http://www.unidata.ucar.edu/software/netcdf-java/documentation.htm
    .. [2] http://thredds.ucar.edu/thredds/catalog.html

    """
    # test for non empty kwargs
    _test_arguments(kwargs)

    # create metadata retrieval object
    filemetadata = FileMetadata('nexrad_cdm', field_names,
                                additional_metadata, file_field_names,
                                exclude_fields, include_fields)

    # open the file
    dataset = netCDF4.Dataset(filename)
    dattrs = dataset.ncattrs()
    dvars = dataset.variables
    if 'cdm_data_type' not in dattrs or dataset.cdm_data_type != 'RADIAL':
        raise IOError('%s is not a valid CDM NetCDF file' % (filename))

    # determine the scan information
    scan_info = _scan_info(dvars)
    radials_per_scan = [max(s['nradials']) for s in scan_info]
    ngates_per_scan = [max(s['ngates']) for s in scan_info]
    ngates = max(ngates_per_scan)
    nrays = sum(radials_per_scan)
    nsweeps = len(scan_info)

    # extract data which changes depending on scan,
    # specifically time, azimuth, elevation and fixed angle data, as well as
    # the moment data.
    time_data = np.empty((nrays, ), dtype='float64')
    azim_data = np.empty((nrays, ), dtype='float32')
    elev_data = np.empty((nrays, ), dtype='float32')
    fixed_agl_data = np.empty((nsweeps, ), dtype='float32')
    fdata = {
        'Reflectivity':
        np.ma.masked_equal(np.ones((nrays, ngates), dtype='float32'), 1),
        'RadialVelocity':
        np.ma.masked_equal(np.ones((nrays, ngates), dtype='float32'), 1),
        'SpectrumWidth':
        np.ma.masked_equal(np.ones((nrays, ngates), dtype='float32'), 1),
        'DifferentialReflectivity':
        np.ma.masked_equal(np.ones((nrays, ngates), dtype='float32'), 1),
        'CorrelationCoefficient':
        np.ma.masked_equal(np.ones((nrays, ngates), dtype='float32'), 1),
        'DifferentialPhase':
        np.ma.masked_equal(np.ones((nrays, ngates), dtype='float32'), 1),
    }

    ray_i = 0
    for scan_index, scan_dic in enumerate(scan_info):

        var_index = scan_dic['index'][0]
        nradials = scan_dic['nradials'][0]

        time_var = scan_dic['time_vars'][0]
        azimuth_var = scan_dic['azimuth_vars'][0]
        elevation_var = scan_dic['elevation_vars'][0]
        nradials = scan_dic['nradials'][0]

        end = ray_i + nradials
        time_data[ray_i:end] = dvars[time_var][var_index][:nradials]
        azim_data[ray_i:end] = dvars[azimuth_var][var_index][:nradials]
        elev_data[ray_i:end] = dvars[elevation_var][var_index][:nradials]
        fixed_agl_data[scan_index] = np.mean(
            dvars[elevation_var][var_index][:nradials])

        for i, moment in enumerate(scan_dic['moments']):

            moment_index = scan_dic['index'][i]
            m_ngates = scan_dic['ngates'][i]
            m_nradials = scan_dic['nradials'][i]

            if moment.endswith('_HI'):
                fdata_name = moment[:-3]
            else:
                fdata_name = moment

            sweep = _get_moment_data(dvars[moment], moment_index, m_ngates)
            fdata[fdata_name][ray_i:ray_i + m_nradials, :m_ngates] = (
                sweep[:m_nradials, :m_ngates])

        ray_i += nradials

    # time
    time = filemetadata('time')
    first_time_var = scan_info[0]['time_vars'][0]
    time_start = datetime.strptime(dvars[first_time_var].units[-20:],
                                   "%Y-%m-%dT%H:%M:%SZ")
    time_start = time_start + timedelta(seconds=int(time_data[0]/1000))
    time['data'] = time_data/1000. - int(time_data[0]/1000)
    time['units'] = make_time_unit_str(time_start)

    # range
    _range = filemetadata('range')
    max_ngates_scan_index = ngates_per_scan.index(ngates)
    scan_dic = scan_info[max_ngates_scan_index]
    max_ngates_moment_index = scan_dic['ngates'].index(ngates)
    distance_var = scan_dic['distance_vars'][max_ngates_moment_index]
    _range['data'] = dvars[distance_var][:]
    _range['meters_to_center_of_first_gate'] = _range['data'][0]
    _range['meters_between_gates'] = _range['data'][1] - _range['data'][0]

    # fields
    fields = {}
    for moment_name, moment_data in fdata.items():
        field_name = filemetadata.get_field_name(moment_name)
        field_dic = filemetadata(field_name)
        field_dic['_FillValue'] = get_fillvalue()
        field_dic['data'] = moment_data
        fields[field_name] = field_dic

    # metadata
    metadata = filemetadata('metadata')
    metadata['original_container'] = 'NEXRAD Level II'

    # scan_type
    scan_type = 'ppi'

    # latitude, longitude, altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')

    # use the locations in the NetCDF file is available
    if ((hasattr(dataset, 'StationLatitude') and
         hasattr(dataset, 'StationLongitude') and
         hasattr(dataset, 'StationElevationInMeters'))):
        lat = dataset.StationLatitude
        lon = dataset.StationLongitude
        alt = dataset.StationElevationInMeters
    else:
        # if no locations in the file look them up from station name.
        if station is None:
            # determine the station name from the filename
            # this will fail in some cases, in which case station
            # should be implicitly provided in the function call.
            station = os.path.basename(filename)[:4].upper()
        lat, lon, alt = get_nexrad_location(station)
    latitude['data'] = np.array([lat], dtype='float64')
    longitude['data'] = np.array([lon], dtype='float64')
    altitude['data'] = np.array([alt], dtype='float64')

    # sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index
    # sweep_end_ray_index
    sweep_number = filemetadata('sweep_number')
    sweep_mode = filemetadata('sweep_mode')
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')

    sweep_number['data'] = np.arange(nsweeps, dtype='int32')
    sweep_mode['data'] = np.array(
        nsweeps * ['azimuth_surveillance'], dtype='S')
    rays_per_scan = list(radials_per_scan)
    sweep_end_ray_index['data'] = np.cumsum(rays_per_scan, dtype='int32') - 1

    rays_per_scan.insert(0, 0)
    sweep_start_ray_index['data'] = np.cumsum(rays_per_scan[:-1],
                                              dtype='int32')

    # azimuth, elevation, fixed_angle
    azimuth = filemetadata('azimuth')
    elevation = filemetadata('elevation')
    fixed_angle = filemetadata('fixed_angle')
    azimuth['data'] = azim_data
    elevation['data'] = elev_data
    fixed_angle['data'] = fixed_agl_data

    dataset.close()
    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=None)


def _scan_info(dvars):
    """ Return a list of information on the scans in the volume. """
    # determine the time of the sweep start
    time_variables = [k for k in dvars.keys() if k.startswith('time')]

    scan_start_times = set([])
    for var in time_variables:
        for time in dvars[var][:, 0]:
            scan_start_times.add(time)
    scan_start_times = list(scan_start_times)
    scan_start_times.sort()

    # build the scan_info list
    time_var_to_moment = {          # time variable to moment conversion
        'timeR': 'Reflectivity',
        'timeV': 'RadialVelocity',
        'timeD': 'DifferentialReflectivity',
        'timeC': 'CorrelationCoefficient',
        'timeP': 'DifferentialPhase',
        'timeR_HI': 'Reflectivity_HI',
        'timeV_HI': 'RadialVelocity_HI',
        'timeD_HI': 'DifferentialReflectivity_HI',
        'timeC_HI': 'CorrelationCoefficient_HI',
        'timeP_HI': 'DifferentialPhase_HI',
    }

    scan_info = [{'start_time': t, 'time_vars': [], 'moments': [],
                  'nradials': [], 'ngates': [], 'elevation_vars': [],
                  'azimuth_vars': [], 'distance_vars': [], 'index': []}
                 for t in scan_start_times]

    for time_var in time_variables:
        for time_var_i, time in enumerate(dvars[time_var][:, 0]):

            scan_index = scan_start_times.index(time)
            scan_dic = scan_info[scan_index]
            moment = time_var_to_moment[time_var]
            _populate_scan_dic(scan_dic, time_var, time_var_i, moment, dvars)

            # corner cases, timeV is a dimension for RadialVelocity AND
            # SpectrumWidth
            if time_var == 'timeV':
                _populate_scan_dic(scan_dic, time_var, time_var_i,
                                   'SpectrumWidth', dvars)
            if time_var == 'timeV_HI':
                _populate_scan_dic(scan_dic, time_var, time_var_i,
                                   'SpectrumWidth_HI', dvars)

    return scan_info


def _populate_scan_dic(scan_dic, time_var, time_var_i, moment, dvars):
    """ Populate a dictionary in the scan_info list. """
    if time_var.endswith('HI'):
        var_suffix = time_var[-4:]
    else:
        var_suffix = time_var[-1:]

    ngates = dvars['numGates' + var_suffix][time_var_i]
    nradials = dvars['numRadials' + var_suffix][time_var_i]

    scan_dic['time_vars'].append(time_var)
    scan_dic['index'].append(time_var_i)
    scan_dic['moments'].append(moment)
    scan_dic['elevation_vars'].append('elevation' + var_suffix)
    scan_dic['azimuth_vars'].append('azimuth' + var_suffix)
    scan_dic['distance_vars'].append('distance' + var_suffix)
    scan_dic['ngates'].append(ngates)
    scan_dic['nradials'].append(nradials)
    return


def _get_moment_data(moment_var, index, ngates):
    """ Retieve moment data for a given scan. """

    # mask, scale and offset
    moment_var.set_auto_maskandscale(False)
    raw_moment_data = moment_var[index][:, :ngates]
    if '_Unsigned' in moment_var.ncattrs():
        if raw_moment_data.dtype == np.int8:
            raw_moment_data = raw_moment_data.view('uint8')
        if raw_moment_data.dtype == np.int16:
            raw_moment_data = raw_moment_data.view('uint16')

    raw_moment_data = np.ma.masked_less_equal(raw_moment_data, 1)

    if 'scale_factor' in moment_var.ncattrs():
        scale = moment_var.scale_factor
    else:
        scale = 1.0

    if 'add_offset' in moment_var.ncattrs():
        add_offset = moment_var.add_offset
    else:
        add_offset = 0.0

    return raw_moment_data * scale + add_offset
