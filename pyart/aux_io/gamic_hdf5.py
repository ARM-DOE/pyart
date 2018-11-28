"""
pyart.aux_io.read_gamic
=======================

Utilities for reading gamic hdf5 files.

.. autosummary::
    :toctree: generated/

    read_gamic
    _get_instrument_params
    _avg_radial_angles
    _prt_mode_from_unfolding

"""

# TODO to move out of aux_io namespace:
# * unit tests - need small sample files.
# * auto-detect file type with pyart.io.read function
# * move to pyart.io namespace

import datetime
import warnings

import numpy as np

from ..config import FileMetadata, get_fillvalue
from ..io.common import make_time_unit_str, _test_arguments
from ..core.radar import Radar
try:
    from .gamicfile import GAMICFile
    _H5PY_AVAILABLE = True
except ImportError:
    _H5PY_AVAILABLE = False
from ..exceptions import MissingOptionalDependency


LIGHT_SPEED = 2.99792458e8  # speed of light in meters per second


def read_gamic(filename, field_names=None, additional_metadata=None,
               file_field_names=False, exclude_fields=None,
               include_fields=None, valid_range_from_file=True,
               units_from_file=True, pulse_width=None, **kwargs):
    """
    Read a GAMIC hdf5 file.

    Parameters
    ----------
    filename : str
        Name of GAMIC HDF5 file to read data from.
    field_names : dict, optional
        Dictionary mapping field names in the file names to radar field names.
        Unlike other read functions, fields not in this dictionary or having a
        value of None are still included in the radar.fields dictionary, to
        exclude them use the `exclude_fields` parameter. Fields which are
        mapped by this dictionary will be renamed from key to value.
    additional_metadata : dict of dicts, optional
        This parameter is not used, it is included for uniformity.
    file_field_names : bool, optional
        True to force the use of the field names from the file in which
        case the `field_names` parameter is ignored. False will use to
        `field_names` parameter to rename fields.
    exclude_fields : list or None, optional
        List of fields to exclude from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters. Set
        to None to include all fields specified by include_fields.
    include_fields : list or None, optional
        List of fields to include from the radar object. This is applied
        after the `file_field_names` and `field_names` parameters. Set
        to None to include all fields not specified by exclude_fields.
    valid_range_from_file : bool, optional
        True to extract valid range (valid_min and valid_max) for all
        field from the file when they are present.  False will not extract
        these parameters.
    units_from_file : bool, optional
        True to extract the units for all fields from the file when available.
        False will not extract units using the default units for the fields.
    pulse_width : list or None,
        Mandatory for gamic radar processors which have pulsewidth enums.
        pulse_width should contain the pulsewidth' in us.

    Returns
    -------
    radar : Radar
        Radar object.

    """
    # check that h5py is available
    if not _H5PY_AVAILABLE:
        raise MissingOptionalDependency(
            "h5py is required to use read_gamic but is not installed")

    # test for non empty kwargs
    _test_arguments(kwargs)

    # create metadata retrieval object
    filemetadata = FileMetadata('gamic', field_names, additional_metadata,
                                file_field_names, exclude_fields,
                                include_fields)

    # Open HDF5 file and get handle
    gfile = GAMICFile(filename)

    # verify that all scans are present in file
    assert gfile.is_file_complete()

    # latitude, longitude and altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')
    latitude['data'] = gfile.where_attr('lat', 'float64')
    longitude['data'] = gfile.where_attr('lon', 'float64')
    altitude['data'] = gfile.where_attr('height', 'float64')

    # metadata
    metadata = filemetadata('metadata')
    metadata['original_container'] = 'GAMIC-HDF5'
    what_mapping = {
        'version': 'gamic_version', 'sets_scheduled': 'sets_scheduled',
        'object': 'gamic_object', 'date': 'gamic_date'}
    for gamic_key, metadata_key in what_mapping.items():
        if gfile.is_attr_in_group('what', gamic_key):
            metadata[metadata_key] = gfile.raw_group_attr('what', gamic_key)

    how_keys = ['software', 'template_name', 'site_name', 'host_name',
                'azimuth_beam', 'elevation_beam', 'sdp_name', 'sw_version',
                'sdp_version', 'simulated']
    for key in how_keys:
        if gfile.is_attr_in_group('how', key):
            metadata[key] = gfile.raw_group_attr('how', key)

    # sweep_start_ray_index, sweep_end_ray_index
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')
    sweep_start_ray_index['data'] = gfile.start_ray.astype('int32')
    sweep_end_ray_index['data'] = gfile.end_ray.astype('int32')

    # sweep number
    sweep_number = filemetadata('sweep_number')
    try:
        sweep_number['data'] = gfile.what_attrs('set_idx', 'int32')
    except KeyError:
        sweep_number['data'] = np.arange(gfile.nsweeps, dtype='int32')

    # sweep_type
    scan_type = gfile.raw_scan0_group_attr('what', 'scan_type').lower()
    scan_type = scan_type.decode('utf-8')
    # check that all scans in the volume are the same type
    if not gfile.is_file_single_scan_type():
        raise NotImplementedError('Mixed scan_type volume.')
    if scan_type not in ['ppi', 'rhi']:
        message = "Unknown scan type: %s, reading as RHI scans." % (scan_type)
        warnings.warn(message)
        scan_type = 'rhi'

    # sweep_mode, fixed_angle
    sweep_mode = filemetadata('sweep_mode')
    fixed_angle = filemetadata('fixed_angle')
    if scan_type == 'rhi':
        sweep_mode['data'] = np.array(gfile.nsweeps * ['rhi'])
        fixed_angle['data'] = gfile.how_attrs('azimuth', 'float32')
    elif scan_type == 'ppi':
        sweep_mode['data'] = np.array(gfile.nsweeps * ['azimuth_surveillance'])
        fixed_angle['data'] = gfile.how_attrs('elevation', 'float32')

    # time
    time = filemetadata('time')
    t_data = gfile.ray_header('timestamp', 'int64')
    start_epoch = t_data[0] // 1.e6     # truncate to second resolution
    start_time = datetime.datetime.utcfromtimestamp(start_epoch)
    time['units'] = make_time_unit_str(start_time)
    time['data'] = ((t_data - start_epoch * 1.e6) / 1.e6).astype('float64')

    # range
    _range = filemetadata('range')
    ngates = int(gfile.raw_scan0_group_attr('how', 'bin_count'))
    range_start = float(gfile.raw_scan0_group_attr('how', 'range_start'))
    #n_samples insertion
    range_samples = int(gfile.raw_scan0_group_attr('how', 'range_samples')) 
    range_step = float(gfile.raw_scan0_group_attr('how', 'range_step'))*range_samples
    # range_step may need to be scaled by range_samples
    # XXX This gives distances to start of gates not center, this matches
    # Radx but may be incorrect, add range_step / 2. for center
    _range['data'] = (np.arange(ngates, dtype='float32') * range_step +
                      range_start)
    _range['meters_to_center_of_first_gate'] = range_start
    _range['meters_between_gates'] = range_step

    # elevation
    elevation = filemetadata('elevation')
    start_angle = gfile.ray_header('elevation_start', 'float32')
    stop_angle = gfile.ray_header('elevation_stop', 'float32')
    elevation['data'] = _avg_radial_angles(start_angle, stop_angle)

    # azimuth
    azimuth = filemetadata('azimuth')
    start_angle = gfile.ray_header('azimuth_start', 'float32')
    stop_angle = gfile.ray_header('azimuth_stop', 'float32')
    azimuth['data'] = _avg_radial_angles(start_angle, stop_angle) % 360.

    # fields
    fields = {}
    moment_groups = gfile.moment_groups()
    moment_names = gfile.moment_names(moment_groups)

    for moment_name, group in zip(moment_names, moment_groups):

        field_name = filemetadata.get_field_name(moment_name)
        if field_name is None:
            continue

        field_dic = filemetadata(field_name)
        field_dic['data'] = gfile.moment_data(group, 'float32')
        field_dic['_FillValue'] = get_fillvalue()

        if valid_range_from_file:
            try:
                valid_min = gfile.raw_scan0_group_attr(group, 'dyn_range_min')
                valid_max = gfile.raw_scan0_group_attr(group, 'dyn_range_max')
                field_dic['valid_min'] = valid_min.decode('utf-8')
                field_dic['valid_max'] = valid_max.decode('utf-8')
            except:
                pass

        if units_from_file:
            try:
                units = gfile.raw_scan0_group_attr(group, 'unit')
                field_dic['units'] = units.decode('utf-8')
            except:
                pass
        fields[field_name] = field_dic

    # ray_angle_res
    ray_angle_res = filemetadata('ray_angle_res')
    ray_angle_res['data'] = gfile.how_attrs('angle_step', 'float32')

    # rays_are_indexed
    rays_are_indexed = filemetadata('rays_are_indexed')
    rays_are_indexed['data'] = np.array(
        [['false', 'true'][i] for i in gfile.how_attrs('angle_sync', 'uint8')])

    # target_scan_rate
    target_scan_rate = filemetadata('target_scan_rate')
    target_scan_rate['data'] = gfile.how_attrs('scan_speed', 'float32')

    # scan_rate
    scan_rate = filemetadata('scan_rate')
    if scan_type == 'ppi':
        azs_names = ['az_speed', 'azimuth_speed']
        azs_name = azs_names[0]
        for azs_name in azs_names:
            if gfile.is_field_in_ray_header(azs_name):
                break
        scan_rate['data'] = gfile.ray_header(azs_name, 'float32')
    elif scan_type == 'rhi':
        els_names = ['el_speed', 'elevation_speed']
        els_name = els_names[0]
        for els_name in els_names:
            if gfile.is_field_in_ray_header(els_name):
                break
        scan_rate['data'] = gfile.ray_header(els_name, 'float32')
    else:
        scan_rate = None

    # instrument_parameters
    instrument_parameters = _get_instrument_params(gfile, filemetadata,
                                                   pulse_width)

    gfile.close()

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters,
        ray_angle_res=ray_angle_res,
        rays_are_indexed=rays_are_indexed,
        scan_rate=scan_rate,
        target_scan_rate=target_scan_rate)


def _get_instrument_params(gfile, filemetadata, pulse_width):
    """ Return a dictionary containing instrument parameters. """

    instrument_params = {}

    dic = filemetadata('frequency')
    dic['data'] = np.array(
        [LIGHT_SPEED / gfile.raw_scan0_group_attr('how', 'radar_wave_length')],
        dtype='float32')
    instrument_params['frequency'] = dic

    dic = filemetadata('radar_beam_width_h')
    dic['data'] = gfile.how_attr('azimuth_beam', 'float32')
    instrument_params['radar_beam_width_h'] = dic

    dic = filemetadata('radar_beam_width_v')
    dic['data'] = gfile.how_attr('elevation_beam', 'float32')
    instrument_params['radar_beam_width_v'] = dic

    dic = filemetadata('pulse_width')
    pw_names = ['pulse_width_us', 'pulse_width_mks', 'pulse_width']
    pw_name = 'pulse_width_us'
    for pw_name in pw_names:
        if gfile.is_attr_in_group('/scan0/how', pw_name):
            break
    if pw_name == 'pulse_width':
        if not pulse_width:
            message = ("read_gamic() is missing 'pulse_width' "
                       "keyword argument")
            raise TypeError(message)
        dic['data'] = gfile.sweep_expand(
            pulse_width[gfile.how_attrs(pw_name, 'int')[0]] * 1e-6)
    else:
        dic['data'] = gfile.sweep_expand(
            gfile.how_attrs(pw_name, 'float32') * 1e-6)
    instrument_params['pulse_width'] = dic

    dic = filemetadata('prt')
    dic['data'] = gfile.sweep_expand(1. / gfile.how_attrs('PRF', 'float32'))
    instrument_params['prt'] = dic

    unfolding = gfile.how_attrs('unfolding', 'int32')
    dic = filemetadata('prt_mode')
    dic['data'] = np.array([_prt_mode_from_unfolding(i) for i in unfolding])
    instrument_params['prt_mode'] = dic

    dic = filemetadata('prt_ratio')
    dic['data'] = gfile.sweep_expand(
        [[1, 2./3., 3./4., 4./5.][i] for i in unfolding])
    instrument_params['prt_ratio'] = dic

    dic = filemetadata('unambiguous_range')
    dic['data'] = gfile.sweep_expand(gfile.how_attrs('range', 'float32'))
    instrument_params['unambiguous_range'] = dic

    if gfile.is_attr_in_group('/scan0/how/extended', 'nyquist_velocity'):
        dic = filemetadata('nyquist_velocity')
        dic['data'] = gfile.sweep_expand(
            gfile.how_ext_attrs('nyquist_velocity'))
        instrument_params['nyquist_velocity'] = dic

    dic = filemetadata('n_samples')
    dic['data'] = gfile.sweep_expand(
        gfile.how_attrs('range_samples', 'int32') *
        gfile.how_attrs('time_samples', 'int32'), dtype='int32')
    instrument_params['n_samples'] = dic

    return instrument_params


def _avg_radial_angles(angle1, angle2):
    """ Return the average angle between two radial angles. """
    return np.angle(
        (np.exp(1.j*np.deg2rad(angle1)) +
         np.exp(1.j*np.deg2rad(angle2))) / 2., deg=True)


def _prt_mode_from_unfolding(unfolding):
    """ Return 'fixed' or 'staggered' depending on unfolding flag """
    if unfolding == 0:
        return 'fixed'
    else:
        return 'staggered'
