"""
pyart.aux_io.rainbow
====================

Routines for reading RAINBOW files (Used by SELEX) using the wradlib library

.. autosummary::
    :toctree: generated/

    read_rainbow_wrl
    _get_angle
    _get_data
    _get_time

"""

# specific modules for this function
import os

try:
    import wradlib
    _WRADLIB_AVAILABLE = True
    # `read_rainbow` as of wradlib version 1.0.0
    try:
        from wradlib.io import read_Rainbow as read_rainbow
    except ImportError:
        from wradlib.io import read_rainbow
except:
    _WRADLIB_AVAILABLE = False


import datetime

import numpy as np

from ..config import FileMetadata, get_fillvalue
from ..io.common import make_time_unit_str, _test_arguments
from ..core.radar import Radar
from ..exceptions import MissingOptionalDependency

RAINBOW_FIELD_NAMES = {
    'W': 'spectrum_width',
    'Wv': 'spectrum_width_vv',  # non standard name
    'Wu': 'unfiltered_spectrum_width',  # non standard name
    'Wvu': 'unfiltered_spectrum_width_vv',  # non standard name
    'V': 'velocity',
    'Vv': 'velocity_vv',  # non standard name
    'Vu': 'unfiltered_velocity',  # non standard name
    'Vvu': 'unfiltered_velocity_vv',  # non standard name
    'dBZ': 'reflectivity',
    'dBZv': 'reflectivity_vv',       # non standard name
    'dBuZ': 'unfiltered_reflectivity',  # non standard name
    'dBuZv': 'unfiltered_reflectivity_vv',  # non standard name
    'ZDR': 'differential_reflectivity',
    'ZDRu': 'unfiltered_differential_reflectivity',  # non standard name
    'RhoHV': 'cross_correlation_ratio',
    'RhoHVu': 'unfiltered_cross_correlation_ratio',  # non standard name
    'PhiDP': 'differential_phase',
    'uPhiDP': 'uncorrected_differential_phase',  # non standard name
    'uPhiDPu':
        'uncorrected_unfiltered_differential_phase',  # non standard name
    'KDP': 'specific_differential_phase',
    'uKDP': 'uncorrected_specific_differential_phase',  # non standard name
    'uKDPu':                                            # non standard name
        'uncorrected_unfiltered_specific_differential_phase',
    'SQI': 'signal_quality_index',  # non standard name
    'SQIv': 'signal_quality_index_vv',  # non standard name
    'SQIu': 'unfiltered_signal_quality_index',  # non standard name
    'SQIvu': 'unfiltered_signal_quality_index_vv',  # non standard name
    'TEMP': 'temperature',  # non standard name
    'ISO0': 'iso0',  # non standard name
}


def read_rainbow_wrl(filename, field_names=None, additional_metadata=None,
                     file_field_names=False, exclude_fields=None,
                     include_fields=None, **kwargs):
    """
    Read a RAINBOW file.
    This routine has been tested to read rainbow5 files version 5.22.3,
    5.34.16 and 5.35.1.
    Since the rainbow file format is evolving constantly there is no guaranty
    that it can work with other versions.
    If necessary, the user should adapt to code according to its own
    file version and raise an issue upstream.

    Data types read by this routine:
    Reflectivity: dBZ, dBuZ, dBZv, dBuZv
    Velocity: V, Vu, Vv, Vvu
    Spectrum width: W, Wu, Wv, Wvu
    Differential reflectivity: ZDR, ZDRu
    Co-polar correlation coefficient: RhoHV, RhoHVu
    Co-polar differential phase: PhiDP, uPhiDP, uPhiDPu
    Specific differential phase: KDP, uKDP, uKDPu
    Signal quality parameters: SQI, SQIu, SQIv, SQIvu
    Temperature: TEMP
    Position of the range bin respect to the ISO0: ISO0

    Parameters
    ----------
    filename : str
        Name of the RAINBOW file to read.
    field_names : dict, optional
        Dictionary mapping RAINBOW field names to radar field names. If a
        data type found in the file does not appear in this dictionary or has
        a value of None it will not be placed in the radar.fields dictionary.
        A value of None, the default, will use the mapping defined in the
        Py-ART configuration file.
    additional_metadata : dict of dicts, optional
        Dictionary of dictionaries to retrieve metadata during this read.
        This metadata is not used during any successive file reads unless
        explicitly included.  A value of None, the default, will not
        introduct any addition metadata and the file specific or default
        metadata as specified by the Py-ART configuration file will be used.
    file_field_names : bool, optional
        True to use the MDV data type names for the field names. If this
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
        Radar object containing data from RAINBOW file.

    """

    # check that wradlib is available
    if not _WRADLIB_AVAILABLE:
        raise MissingOptionalDependency(
            "wradlib is required to use read_rainbow_wrl but is not installed")

    # test for non empty kwargs
    _test_arguments(kwargs)

    # check if it is the right file. Open it and read it
    bfile = os.path.basename(filename)
    supported_file = (bfile.endswith('.vol') or bfile.endswith('.azi') or
                      bfile.endswith('.ele'))
    if not supported_file:
        raise ValueError(
            'Only data files with extension .vol, .azi or .ele are supported')

    # create metadata retrieval object
    if field_names is None:
        field_names = RAINBOW_FIELD_NAMES
    filemetadata = FileMetadata('RAINBOW', field_names, additional_metadata,
                                file_field_names, exclude_fields,
                                include_fields)

    rbf = read_rainbow(filename, loaddata=True)

    # check the number of slices
    nslices = int(rbf['volume']['scan']['pargroup']['numele'])
    if nslices > 1:
        single_slice = False
        common_slice_info = rbf['volume']['scan']['slice'][0]
    else:
        single_slice = True
        common_slice_info = rbf['volume']['scan']['slice']

    # check the data type
    # all slices should have the same data type
    datatype = common_slice_info['slicedata']['rawdata']['@type']
    field_name = filemetadata.get_field_name(datatype)
    if field_name is None:
        raise ValueError('Field Name Unknown')

    # get definitions from filemetadata class
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')
    metadata = filemetadata('metadata')
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')
    sweep_number = filemetadata('sweep_number')
    sweep_mode = filemetadata('sweep_mode')
    fixed_angle = filemetadata('fixed_angle')
    elevation = filemetadata('elevation')
    _range = filemetadata('range')
    azimuth = filemetadata('azimuth')
    _time = filemetadata('time')
    field_dic = filemetadata(field_name)

    # other metadata
    scan_rate = filemetadata('scan_rate')
    frequency = filemetadata('frequency')

    # get general file information

    # position and radar frequency
    if 'sensorinfo' in rbf['volume'].keys():
        latitude['data'] = np.array(
            [rbf['volume']['sensorinfo']['lat']], dtype='float64')
        longitude['data'] = np.array(
            [rbf['volume']['sensorinfo']['lon']], dtype='float64')
        altitude['data'] = np.array(
            [rbf['volume']['sensorinfo']['alt']], dtype='float64')
        frequency['data'] = np.array(
            [3e8 / float(rbf['volume']['sensorinfo']['wavelen'])],
            dtype='float64')
    elif 'radarinfo' in rbf['volume'].keys():
        latitude['data'] = np.array(
            [rbf['volume']['radarinfo']['@lat']], dtype='float64')
        longitude['data'] = np.array(
            [rbf['volume']['radarinfo']['@lon']], dtype='float64')
        altitude['data'] = np.array(
            [rbf['volume']['radarinfo']['@alt']], dtype='float64')
        frequency['data'] = np.array(
            [3e8 / float(rbf['volume']['radarinfo']['wavelen'])],
            dtype='float64')

    # antenna speed
    if 'antspeed' in common_slice_info:
        ant_speed = float(common_slice_info['antspeed'])
    else:
        ant_speed = 10.
        print('WARNING: Unable to read antenna speed. Default value of ' +
              str(ant_speed) + ' deg/s will be used')

    # angle step
    angle_step = float(common_slice_info['anglestep'])

    # sweep_number (is the sweep index)
    sweep_number['data'] = np.arange(nslices, dtype='int32')

    # get number of rays and number of range bins per sweep
    rays_per_sweep = np.empty(nslices, dtype='int32')

    if single_slice:
        rays_per_sweep[0] = int(
            common_slice_info['slicedata']['rawdata']['@rays'])
        nbins = int(common_slice_info['slicedata']['rawdata']['@bins'])
        ssri = np.array([0], dtype='int32')
        seri = np.array([rays_per_sweep[0] - 1], dtype='int32')
    else:
        # number of range bins per ray in sweep
        nbins_sweep = np.empty(nslices, dtype='int32')
        for i in range(nslices):
            slice_info = rbf['volume']['scan']['slice'][i]
            # number of rays per sweep
            rays_per_sweep[i] = int(
                slice_info['slicedata']['rawdata']['@rays'])

            # number of range bins per ray in sweep
            nbins_sweep[i] = int(
                slice_info['slicedata']['rawdata']['@bins'])

        # all sweeps have to have the same number of range bins
        if any(nbins_sweep != nbins_sweep[0]):
            raise ValueError('number of range bins changes between sweeps')
        nbins = nbins_sweep[0]
        ssri = np.cumsum(np.append([0], rays_per_sweep[:-1])).astype('int32')
        seri = np.cumsum(rays_per_sweep).astype('int32') - 1

    # total number of rays and sweep start ray index and end
    total_rays = sum(rays_per_sweep)
    sweep_start_ray_index['data'] = ssri
    sweep_end_ray_index['data'] = seri

    # range
    r_res = float(common_slice_info['rangestep']) * 1000.
    if 'start_range' in common_slice_info.keys():
        start_range = float(common_slice_info['start_range']) * 1000.
    else:
        start_range = 0.
    _range['data'] = np.linspace(
        start_range+r_res / 2., float(nbins - 1.) * r_res+r_res / 2.,
        nbins).astype('float32')

    # containers for data
    t_fixed_angle = np.empty(nslices, dtype='float64')
    moving_angle = np.empty(total_rays, dtype='float64')
    static_angle = np.empty(total_rays, dtype='float64')
    time_data = np.empty(total_rays, dtype='float64')
    fdata = np.ma.zeros((total_rays, nbins), dtype='float32',
                        fill_value=get_fillvalue())

    # read data from file
    if bfile.endswith('.vol') or bfile.endswith('.azi'):
        scan_type = 'ppi'
        sweep_mode['data'] = np.array(nslices * ['azimuth_surveillance'])
    else:
        scan_type = 'rhi'
        sweep_mode['data'] = np.array(['elevation_surveillance'])

    # read data from file:
    for i in range(nslices):
        if single_slice:
            slice_info = common_slice_info
        else:
            slice_info = rbf['volume']['scan']['slice'][i]

        # fixed angle
        t_fixed_angle[i] = float(slice_info['posangle'])

        # fixed angle (repeated for each ray)
        static_angle[ssri[i]: seri[i]+1] = t_fixed_angle[i]

        # moving angle
        moving_angle[ssri[i]: seri[i]+1], angle_start, angle_stop = (
            _get_angle(slice_info['slicedata']['rayinfo'],
                       angle_step=angle_step, scan_type=scan_type))

        # time
        time_data[ssri[i]:seri[i]+1], sweep_start_epoch = (
                _get_time(slice_info['slicedata']['@date'],
                          slice_info['slicedata']['@time'],
                          angle_start[0], angle_stop[-1], angle_step,
                          rays_per_sweep[i], ant_speed, scan_type=scan_type))

        if i == 0:
            volume_start_epoch = sweep_start_epoch + 0.
            start_time = (
                datetime.datetime.utcfromtimestamp(volume_start_epoch))

        # data
        fdata[ssri[i]:seri[i]+1, :] = _get_data(
            slice_info['slicedata']['rawdata'],
            rays_per_sweep[i], nbins)

    if bfile.endswith('.vol') or bfile.endswith('.azi'):
        azimuth['data'] = moving_angle
        elevation['data'] = static_angle
    else:
        azimuth['data'] = static_angle
        elevation['data'] = moving_angle

    fixed_angle['data'] = t_fixed_angle

    _time['data'] = time_data-volume_start_epoch
    _time['units'] = make_time_unit_str(start_time)

    # fields
    fields = {}
    # create field dictionary
    field_dic['_FillValue'] = get_fillvalue()
    field_dic['data'] = fdata
    fields[field_name] = field_dic

    # metadata
    # metadata['instrument_name'] = radar_id

    # instrument_parameters
    instrument_parameters = dict()
    instrument_parameters.update({'frequency': frequency})

    return Radar(_time, _range, fields, metadata, scan_type, latitude,
                 longitude, altitude, sweep_number, sweep_mode, fixed_angle,
                 sweep_start_ray_index, sweep_end_ray_index, azimuth,
                 elevation, instrument_parameters=instrument_parameters)


def _get_angle(ray_info, angle_step=None, scan_type='ppi'):
    """
    obtains the ray angle start, stop and center

    Parameters
    ----------
    ray_info : dictionary of dictionaries
        contains the ray info
    angle_step : float
        Optional. The angle step. Used in case there is no information of
        angle stop. Otherwise ignored.
    scan_type : str
        Default ppi. scan_type. Either ppi or rhi.

    Returns
    -------
    moving_angle : numpy array
        the central point of the angle [Deg]
    angle_start :
        the starting point of the angle [Deg]
    angle_stop :
        the end point of the angle [Deg]

    """
    bin_to_deg = 360./65536.

    def _extract_angles(data):
        angle = np.array(data * bin_to_deg, dtype='float64')
        if scan_type == 'rhi':
            ind = (angle > 225.).nonzero()
            angle[ind] -= 360.
        return angle

    try:
        angle_start = _extract_angles(ray_info['data'])
        if angle_step is None:
            raise ValueError('Unknown angle step')
        angle_stop = angle_start + angle_step
    except TypeError:
        angle_start = _extract_angles(ray_info[0]['data'])
        angle_stop = _extract_angles(ray_info[1]['data'])

    moving_angle = np.angle((np.exp(1.j * np.deg2rad(angle_start)) +
                            np.exp(1.j * np.deg2rad(angle_stop))) / 2.,
                            deg=True)
    moving_angle[moving_angle < 0.] += 360.  # [0, 360]

    return moving_angle, angle_start, angle_stop


def _get_data(rawdata, nrays, nbins):
    """
    Obtains the raw data

    Parameters
    ----------
    rawdata : dictionary of dictionaries
        contains the raw data information
    nrays : int
        Number of rays in sweep
    nbins : int
        Number of bins in ray

    Returns
    -------
    data : numpy array
        the data

    """
    databin = rawdata['data']
    datamin = float(rawdata['@min'])
    datamax = float(rawdata['@max'])
    datadepth = float(rawdata['@depth'])
    datatype = rawdata['@type']

    data = np.array(datamin+databin*(datamax-datamin) / 2 ** datadepth,
                    dtype='float32')

    # fill invalid data with fill value
    mask = databin == 0
    data[mask.nonzero()] = get_fillvalue()

    # put phidp data in the range [-180, 180]
    if (datatype == 'PhiDP') or (datatype == 'uPhiDP') or (
            datatype == 'uPhiDPu'):
        is_above_180 = data > 180.
        data[is_above_180.nonzero()] -= 360.

    data = np.reshape(data, [nrays, nbins])
    mask = np.reshape(mask, [nrays, nbins])

    masked_data = np.ma.array(data, mask=mask, fill_value=get_fillvalue())

    return masked_data


def _get_time(date_sweep, time_sweep, first_angle_start, last_angle_stop,
              angle_step, nrays, ant_speed, scan_type='ppi'):
    """
    Computes the time at the center of each ray

    Parameters
    ----------
    date_sweep, time_sweep : str
        the date and time of the sweep
    first_angle_start : float
        The starting point of the first angle in the sweep
    last_angle_stop : float
        The end point of the last angle in the sweep
    nrays : int
        Number of rays in sweep
    ant_speed : float
        antenna speed [deg/s]
    scan_type : str
        Default ppi. scan_type. Either ppi or rhi.

    Returns
    -------
    time_data : numpy array
        the time of each ray
    sweep_start_epoch : float
        sweep start time in seconds since 1.1.1970

    """
    datetime_sweep = datetime.datetime.strptime(
        date_sweep+' '+time_sweep, '%Y-%m-%d %H:%M:%S')
    sweep_start_epoch = (
        datetime_sweep - datetime.datetime(1970, 1, 1)).total_seconds()
    if scan_type == 'ppi':
        if (last_angle_stop > first_angle_start) and (
                (last_angle_stop-first_angle_start) /
                nrays > angle_step):
            sweep_duration = (last_angle_stop - first_angle_start) / ant_speed
        else:
            sweep_duration = (
                last_angle_stop + 360. - first_angle_start) / ant_speed
    else:
        if last_angle_stop > first_angle_start:
            sweep_duration = (last_angle_stop - first_angle_start) / ant_speed
        else:
            sweep_duration = (first_angle_start - last_angle_stop) / ant_speed

    time_angle = sweep_duration/nrays

    sweep_end_epoch = sweep_start_epoch+sweep_duration

    time_data = np.linspace(
        sweep_start_epoch+time_angle / 2.,
        sweep_end_epoch-time_angle / 2., num=nrays)

    return time_data, sweep_start_epoch
