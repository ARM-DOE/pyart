"""
pyart.io.uf
===========

Reading of Universal format (UF) files

.. autosummary::
    :toctree: generated/

    read_uf
    _get_scan_type
    _get_instrument_parameters

"""

import warnings

import numpy as np
from netCDF4 import date2num

from ..config import FileMetadata, get_fillvalue
from ..core.radar import Radar
from .common import make_time_unit_str, _test_arguments, prepare_for_read
from .uffile import UFFile

_LIGHT_SPEED = 2.99792458e8  # speed of light in meters per second
_UF_SWEEP_MODES = {
    0: 'calibration',
    1: 'ppi',
    2: 'coplane',
    3: 'rhi',
    4: 'vpt',
    5: 'target',
    6: 'manual',
    7: 'idle',
    8: 'ppi',   # RadX used this to indicate surveillance PPI scans
}

_SWEEP_MODE_STR = {
    'calibration': 'calibration',
    'ppi': 'azimuth_surveillance',
    'coplane': 'coplane',
    'rhi': 'rhi',
    'vpt': 'vertical_pointing',
    'target': 'pointing',
    'manual': 'manual',
    'idle': 'idle',
}


def read_uf(filename, field_names=None, additional_metadata=None,
            file_field_names=False, exclude_fields=None,
            include_fields=None, delay_field_loading=False, **kwargs):
    """
    Read a UF File.

    Parameters
    ----------
    filename : str or file-like
        Name of Universal format file to read data from.
    field_names : dict, optional
        Dictionary mapping UF data type names to radar field names. If a
        data type found in the file does not appear in this dictionary or has
        a value of None it will not be placed in the radar.fields dictionary.
        A value of None, the default, will use the mapping defined in the
        Py-ART configuration file.
    additional_metadata : dict of dicts, optional
        Dictionary of dictionaries to retrieve metadata from during this read.
        This metadata is not used during any successive file reads unless
        explicitly included.  A value of None, the default, will not
        introduce any addition metadata and the file specific or default
        metadata as specified by the Py-ART configuration file will be used.
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
    delay_field_loading : bool
        This option is not implemented in the function but included for
        compatibility.

    Returns
    -------
    radar : Radar
        Radar object.

    """
    # test for non empty kwargs
    _test_arguments(kwargs)

    # create metadata retrieval object
    filemetadata = FileMetadata('uf', field_names, additional_metadata,
                                file_field_names, exclude_fields,
                                include_fields)

    # Open UF file and get handle
    ufile = UFFile(prepare_for_read(filename))
    first_ray = ufile.rays[0]

    # time
    dts = ufile.get_datetimes()
    units = make_time_unit_str(min(dts))
    time = filemetadata('time')
    time['units'] = units
    time['data'] = date2num(dts, units).astype('float32')

    # range
    _range = filemetadata('range')
    # assume that the number of gates and spacing from the first ray is
    # representative of the entire volume
    field_header = first_ray.field_headers[0]
    ngates = field_header['nbins']
    step = field_header['range_spacing_m']
    # this gives distances to the center of each gate, remove step/2 for start
    start = (field_header['range_start_km'] * 1000. +
             field_header['range_start_m'] + step / 2.)
    _range['data'] = np.arange(ngates, dtype='float32') * step + start
    _range['meters_to_center_of_first_gate'] = start
    _range['meters_between_gates'] = step

    # latitude, longitude and altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')
    lat, lon, height = first_ray.get_location()
    latitude['data'] = np.array([lat], dtype='float64')
    longitude['data'] = np.array([lon], dtype='float64')
    altitude['data'] = np.array([height], dtype='float64')

    # metadata
    metadata = filemetadata('metadata')
    metadata['original_container'] = 'UF'
    metadata['site_name'] = first_ray.mandatory_header['site_name']
    metadata['radar_name'] = first_ray.mandatory_header['radar_name']

    # sweep_start_ray_index, sweep_end_ray_index
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')
    sweep_start_ray_index['data'] = ufile.first_ray_in_sweep
    sweep_end_ray_index['data'] = ufile.last_ray_in_sweep

    # sweep number
    sweep_number = filemetadata('sweep_number')
    sweep_number['data'] = np.arange(ufile.nsweeps, dtype='int32')

    # scan_type
    scan_type = _get_scan_type(first_ray)

    # sweep_mode
    sweep_mode = filemetadata('sweep_mode')
    sweep_mode['data'] = np.array(
        ufile.nsweeps * [_SWEEP_MODE_STR[scan_type]], dtype='S')

    # elevation
    elevation = filemetadata('elevation')
    elevation['data'] = ufile.get_elevations()

    # azimuth
    azimuth = filemetadata('azimuth')
    azimuth['data'] = ufile.get_azimuths()

    # fixed_angle
    fixed_angle = filemetadata('fixed_angle')
    fixed_angle['data'] = ufile.get_sweep_fixed_angles()

    # fields
    fields = {}
    for uf_field_number, uf_field_dic in enumerate(first_ray.field_positions):
        uf_field_name = uf_field_dic['data_type'].decode('ascii')
        field_name = filemetadata.get_field_name(uf_field_name)
        if field_name is None:
            continue
        field_dic = filemetadata(field_name)
        field_dic['data'] = ufile.get_field_data(uf_field_number)
        field_dic['_FillValue'] = get_fillvalue()
        fields[field_name] = field_dic

    # instrument_parameters
    instrument_parameters = _get_instrument_parameters(ufile, filemetadata)

    # scan rate
    scan_rate = filemetadata('scan_rate')
    scan_rate['data'] = ufile.get_sweep_rates()

    ufile.close()
    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        scan_rate=scan_rate,
        instrument_parameters=instrument_parameters)


def _get_scan_type(ufray):
    """ Ruturn the scan type of a UF ray. """
    uf_sweep_mode = ufray.mandatory_header['sweep_mode']
    if uf_sweep_mode in _UF_SWEEP_MODES:
        return _UF_SWEEP_MODES[uf_sweep_mode]
    else:
        warnings.warn('Unknown sweep mode, defaulting to ppi')
        return 'ppi'


def _get_instrument_parameters(ufile, filemetadata):
    """ Return a dictionary containing instrument parameters. """

    # pulse width
    pulse_width = filemetadata('pulse_width')
    pulse_width['data'] = ufile.get_pulse_widths() / _LIGHT_SPEED  # m->sec

    # assume that the parameters in the first ray represent the beam widths,
    # bandwidth and frequency in the entire volume
    first_ray = ufile.rays[0]
    field_header = first_ray.field_headers[0]
    beam_width_h = field_header['beam_width_h'] / 64.
    beam_width_v = field_header['beam_width_v'] / 64.
    bandwidth = field_header['bandwidth'] / 16. * 1.e6
    wavelength_cm = field_header['wavelength_cm'] / 64.
    if wavelength_cm == 0:
        warnings.warn('Invalid wavelength, frequency set to default value.')
        wavelength_hz = 9999.0
    else:
        wavelength_hz = _LIGHT_SPEED / (wavelength_cm / 100.)

    # radar_beam_width_h
    radar_beam_width_h = filemetadata('radar_beam_width_h')
    radar_beam_width_h['data'] = np.array([beam_width_h], dtype='float32')

    # radar_beam_width_v
    radar_beam_width_v = filemetadata('radar_beam_width_v')
    radar_beam_width_v['data'] = np.array([beam_width_v], dtype='float32')

    # radar_receiver_bandwidth
    radar_receiver_bandwidth = filemetadata('radar_receiver_bandwidth')
    radar_receiver_bandwidth['data'] = np.array([bandwidth], dtype='float32')

    # polarization_mode
    polarization_mode = filemetadata('polarization_mode')
    polarization_mode['data'] = ufile.get_sweep_polarizations()

    # frequency
    frequency = filemetadata('frequency')
    frequency['data'] = np.array([wavelength_hz], dtype='float32')

    # prt
    prt = filemetadata('prt')
    prt['data'] = ufile.get_prts() / 1e6  # us->sec

    instrument_parameters = {
        'pulse_width': pulse_width,
        'radar_beam_width_h': radar_beam_width_h,
        'radar_beam_width_v': radar_beam_width_v,
        'radar_receiver_bandwidth': radar_receiver_bandwidth,
        'polarization_mode': polarization_mode,
        'frequency': frequency,
        'prt': prt,
    }

    # nyquist velocity if defined
    nyquist_velocity = filemetadata('nyquist_velocity')
    nyquist_velocity['data'] = ufile.get_nyquists()
    if nyquist_velocity['data'] is not None:
        instrument_parameters['nyquist_velocity'] = nyquist_velocity

    return instrument_parameters
