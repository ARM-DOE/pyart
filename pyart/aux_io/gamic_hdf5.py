"""
pyart.aux_io.read_gamic
=======================

Utilities for reading gamic hdf5 files.

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    _GAMICFileHelper

.. autosummary::
    :toctree: generated/

    read_gamic
    _get_instrument_params
    _avg_radial_angles
    _prt_mode_from_unfolding
    _get_gamic_sweep_data


"""

# TODO to move out of aux_io namespace:
# * unit tests - need small sample files.
# * auto-detect file type with pyart.io.read function
# * move to pyart.io namespace

import datetime
import warnings

import h5py
import numpy as np

from ..config import FileMetadata, get_fillvalue
from ..io.common import make_time_unit_str
from ..core.radar import Radar


LIGHT_SPEED = 2.99792458e8  # speed of light in meters per second


def read_gamic(filename, field_names=None, additional_metadata=None,
               file_field_names=False, exclude_fields=None,
               valid_range_from_file=True, units_from_file=True):
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
        after the `file_field_names` and `field_names` parameters.
    valid_range_from_file : bool, optional
        True to extract valid range (valid_min and valid_max) for all
        field from the file when they are present.  False will not extract
        these parameters.
    units_from_file : bool, optional
        True to extract the units for all fields from the file when available.
        False will not extract units using the default units for the fields.

    Returns
    -------
    radar : Radar
        Radar object.

    """
    # create metadata retrieval object
    filemetadata = FileMetadata('gamic', field_names, additional_metadata,
                                file_field_names, exclude_fields)

    # Open HDF5 file and get handle
    hfile = h5py.File(filename, 'r')
    ghelp = _GAMICFileHelper(hfile)

    # verify that scan0 ... scan[nsweeps-1] are present in file
    for scan in ghelp.scans:
        assert scan in hfile

    # latitude, longitude and altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')
    latitude['data'] = ghelp.where_attr('lat', 'float64')
    longitude['data'] = ghelp.where_attr('lon', 'float64')
    altitude['data'] = ghelp.where_attr('height', 'float64')

    # metadata
    metadata = filemetadata('metadata')
    metadata['original_container'] = 'GAMIC-HDF5'
    what_mapping = {
        'version': 'gamic_version', 'sets_scheduled': 'sets_scheduled',
        'object': 'gamic_object', 'date': 'gamic_date'}
    for hkey, mkey in what_mapping.items():
        if hkey in hfile['what'].attrs:
            metadata[mkey] = hfile['what'].attrs[hkey]

    how_keys = ['software', 'template_name', 'site_name', 'host_name',
                'azimuth_beam', 'elevation_beam', 'sdp_name', 'sw_version',
                'sdp_version', 'simulated']
    for key in how_keys:
        if key in hfile['how'].attrs:
            metadata[key] = hfile['how'].attrs[key]

    # sweep_start_ray_index, sweep_end_ray_index
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')
    sweep_start_ray_index['data'] = ghelp.start_ray.astype('int32')
    sweep_end_ray_index['data'] = ghelp.end_ray.astype('int32')

    # sweep number
    sweep_number = filemetadata('sweep_number')
    try:
        sweep_number['data'] = ghelp.what_attrs('set_idx', 'int32')
    except KeyError:
        sweep_number['data'] = np.arange(len(ghelp.scans), dtype='int32')

    # sweep_type
    scan_type = hfile['scan0']['what'].attrs['scan_type'].lower()
    # check that all scans in the volume are the same type
    for scan in ghelp.scans:
        if hfile[scan]['what'].attrs['scan_type'].lower() != scan_type:
            raise NotImplementedError('Mixed scan_type volume in scan:' + scan)
    if scan_type not in ['ppi', 'rhi']:
        message = "Unknown scan type: %s, reading as RHI scans." % (scan_type)
        warnings.warn(message)
        scan_type = 'rhi'

    # sweep_mode, fixed_angle
    sweep_mode = filemetadata('sweep_mode')
    fixed_angle = filemetadata('fixed_angle')
    if scan_type == 'rhi':
        sweep_mode['data'] = np.array(ghelp.nsweeps * ['rhi'])
        fixed_angle['data'] = ghelp.how_attrs('azimuth', 'float32')
    elif scan_type == 'ppi':
        sweep_mode['data'] = np.array(ghelp.nsweeps * ['azimuth_surveillance'])
        fixed_angle['data'] = ghelp.how_attrs('elevation', 'float32')

    # time
    time = filemetadata('time')
    t_data = ghelp.ray_header('timestamp', 'int64')
    start_epoch = t_data[0] // 1.e6     # truncate to second resolution
    start_time = datetime.datetime.utcfromtimestamp(start_epoch)
    time['units'] = make_time_unit_str(start_time)
    time['data'] = ((t_data - start_epoch * 1.e6) / 1.e6).astype('float64')

    # range
    _range = filemetadata('range')
    ngates = int(hfile['/scan0/how'].attrs['bin_count'])
    range_start = float(hfile['/scan0/how'].attrs['range_start'])
    range_step = float(hfile['/scan0/how'].attrs['range_step'])
    # range_step may need to be scaled by range_samples
    # XXX This gives distances to start of gates not center, this matches
    # Radx but may be incorrect, add range_step / 2. for center
    _range['data'] = (np.arange(ngates, dtype='float32') * range_step +
                      range_start)
    _range['meters_to_center_of_first_gate'] = range_start
    _range['meters_between_gates'] = range_step

    # elevation
    elevation = filemetadata('elevation')
    start_angle = ghelp.ray_header('elevation_start', 'float32')
    stop_angle = ghelp.ray_header('elevation_stop', 'float32')
    elevation['data'] = _avg_radial_angles(start_angle, stop_angle)

    # azimuth
    azimuth = filemetadata('azimuth')
    start_angle = ghelp.ray_header('azimuth_start', 'float32')
    stop_angle = ghelp.ray_header('azimuth_stop', 'float32')
    azimuth['data'] = _avg_radial_angles(start_angle, stop_angle) % 360.

    # fields
    fields = {}
    h_keys = [k for k in hfile['/scan0'] if k.startswith('moment_')]
    h_fields = [hfile['/scan0'][k].attrs['moment'] for k in h_keys]

    for h_field, h_key in zip(h_fields, h_keys):

        field_name = filemetadata.get_field_name(h_field)
        if field_name is None:
            continue

        field_dic = filemetadata(field_name)
        field_dic['data'] = ghelp.moment_data(h_key, 'float32')
        field_dic['_FillValue'] = get_fillvalue()

        if valid_range_from_file:
            try:
                valid_min = hfile['/scan0'][h_key].attrs['dyn_range_min']
                valid_max = hfile['/scan0'][h_key].attrs['dyn_range_max']
                field_dic['valid_min'] = valid_min
                field_dic['valid_max'] = valid_max
            except:
                pass

        if units_from_file:
            try:
                field_dic['units'] = hfile['/scan0'][h_key].attrs['unit']
            except:
                pass
        fields[field_name] = field_dic

    # ray_angle_res
    ray_angle_res = filemetadata('ray_angle_res')
    ray_angle_res['data'] = ghelp.how_attrs('angle_step', 'float32')

    # rays_are_indexed
    rays_are_indexed = filemetadata('rays_are_indexed')
    rays_are_indexed['data'] = np.array(
        [['false', 'true'][i] for i in ghelp.how_attrs('angle_sync', 'uint8')])

    # target_scan_rate
    target_scan_rate = filemetadata('target_scan_rate')
    target_scan_rate['data'] = ghelp.how_attrs('scan_speed', 'float32')

    # scan_rate
    scan_rate = filemetadata('scan_rate')
    if scan_type == 'ppi':
        scan_rate['data'] = ghelp.ray_header('az_speed', 'float32')
    elif scan_type == 'rhi':
        scan_rate['data'] = ghelp.ray_header('el_speed', 'float32')
    else:
        scan_rate = None

    # instrument_parameters
    instrument_parameters = _get_instrument_params(ghelp, filemetadata)

    hfile.close()

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


def _get_instrument_params(ghelp, filemetadata):
    """ Return a dictionary containing instrument parameters. """

    instrument_params = {}

    dic = filemetadata('frequency')
    dic['data'] = np.array(
        [LIGHT_SPEED / ghelp.hfile['scan0']['how'].attrs['radar_wave_length']],
        dtype='float32')
    instrument_params['frequency'] = dic

    dic = filemetadata('radar_beam_width_h')
    dic['data'] = ghelp.how_attr('azimuth_beam', 'float32')
    instrument_params['radar_beam_width_h'] = dic

    dic = filemetadata('radar_beam_width_v')
    dic['data'] = ghelp.how_attr('elevation_beam', 'float32')
    instrument_params['radar_beam_width_v'] = dic

    dic = filemetadata('pulse_width')
    dic['data'] = ghelp.sweep_expand(
        ghelp.how_attrs('pulse_width_us', 'float32') * 1e-6)
    instrument_params['pulse_width'] = dic

    dic = filemetadata('prt')
    dic['data'] = ghelp.sweep_expand(1. / ghelp.how_attrs('PRF', 'float32'))
    instrument_params['prt'] = dic

    unfolding = ghelp.how_attrs('unfolding', 'int32')
    dic = filemetadata('prt_mode')
    dic['data'] = np.array([_prt_mode_from_unfolding(i) for i in unfolding])
    instrument_params['prt_mode'] = dic

    dic = filemetadata('prt_ratio')
    dic['data'] = ghelp.sweep_expand(
        [[1, 2./3., 3./4., 4./5.][i] for i in unfolding])
    instrument_params['prt_ratio'] = dic

    dic = filemetadata('unambiguous_range')
    dic['data'] = ghelp.sweep_expand(ghelp.how_attrs('range', 'float32'))
    instrument_params['unambiguous_range'] = dic

    dic = filemetadata('nyquist_velocity')
    dic['data'] = ghelp.sweep_expand(ghelp.how_ext_attrs('nyquist_velocity'))
    instrument_params['nyquist_velocity'] = dic

    dic = filemetadata('n_samples')
    dic['data'] = ghelp.sweep_expand(
        ghelp.how_attrs('range_samples', 'int32') *
        ghelp.how_attrs('time_samples', 'int32'), dtype='int32')
    instrument_params['n_samples'] = dic

    return instrument_params


class _GAMICFileHelper(object):
    """
    A class to help read GAMIC files.

    This class has methods which performs looping routines which are
    useful when reading GAMIC files.
    """

    def __init__(self, hfile):
        """ initialize object. """

        self.hfile = hfile
        self.nsweeps = hfile['what'].attrs['sets']
        self.scans = ['scan%i' % (i) for i in range(self.nsweeps)]
        self.rays_per_sweep = self.how_attrs('ray_count', 'int32')
        self.total_rays = sum(self.rays_per_sweep)
        # starting and ending ray for each sweep
        self.start_ray = np.cumsum(np.append([0], self.rays_per_sweep[:-1]))
        self.end_ray = np.cumsum(self.rays_per_sweep) - 1

        return

    # attribute look up
    def where_attr(self, attr, dtype):
        """ Return an array containing a attribute from the where group. """
        return np.array([self.hfile['where'].attrs[attr]], dtype=dtype)

    def how_attr(self, attr, dtype):
        """ Return an array containing a attribute from the how group. """
        return np.array([self.hfile['how'].attrs[attr]], dtype=dtype)

    # scan/sweep based attribute lookup
    def how_attrs(self, attr, dtype):
        """ Return an array of an attribute for each scan's how group. """
        return np.array([self.hfile[s]['how'].attrs[attr]
                         for s in self.scans], dtype=dtype)

    def how_ext_attrs(self, attr):
        """
        Return a list of an attribute in each scan's how/extended group.
        """
        return [float(self.hfile[s]['how']['extended'].attrs[attr])
                for s in self.scans]

    def what_attrs(self, attr, dtype):
        """ Return a list of an attribute for each scan's what group. """
        return np.array([self.hfile[s]['what'].attrs[attr]
                         for s in self.scans], dtype=dtype)

    # misc looping
    def ray_header(self, field, dtype):
        """ Return an array containing a ray_header field for each sweep. """
        data = np.empty((self.total_rays, ), dtype=dtype)
        for scan, start, end in zip(self.scans, self.start_ray, self.end_ray):
            data[start:end+1] = self.hfile[scan]['ray_header'][field]
        return data

    def moment_data(self, moment, dtype):
        """ Read in moment data from all sweeps. """
        ngates = int(self.hfile['/scan0/how'].attrs['bin_count'])
        data = np.ma.zeros((self.total_rays, ngates), dtype=dtype)
        data[:] = np.ma.masked      # volume data initially all masked
        for scan, start, end in zip(self.scans, self.start_ray, self.end_ray):
            # read in sweep data if field exists in scan.
            if moment in self.hfile[scan]:
                sweep_data = _get_gamic_sweep_data(self.hfile[scan][moment])
                data[start:end+1, :sweep_data.shape[1]] = sweep_data[:]
        return data

    def sweep_expand(self, arr, dtype='float32'):
        """ Expand an sweep indexed array to be ray indexed """
        return np.repeat(arr, self.rays_per_sweep).astype(dtype)


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


def _get_gamic_sweep_data(group):
    """ Get GAMIC HDF5 sweep data from an HDF5 group. """
    dyn_range_min = group.attrs['dyn_range_min']
    dyn_range_max = group.attrs['dyn_range_max']
    raw_data = group[:]
    fmt = group.attrs['format']
    if fmt == 'UV16':
        # unsigned 16-bit integer data, 0 indicates a masked value
        assert raw_data.dtype == np.uint16
        scale = (dyn_range_max - dyn_range_min) / 65535.
        offset = dyn_range_min
        sweep_data = np.ma.masked_array(
            raw_data * scale + offset, mask=(raw_data == 0), dtype='float32')
    elif fmt == 'UV8':
        # unsigned 8-bit integer data, 0 indicates a masked value
        assert raw_data.dtype == np.uint8
        scale = (dyn_range_max - dyn_range_min) / 255.
        offset = dyn_range_min
        sweep_data = np.ma.masked_array(
            raw_data * scale + offset, mask=(raw_data == 0), dtype='float32')
    else:
        raise NotImplementedError('GAMIC data format: %s', fmt)
    return sweep_data
