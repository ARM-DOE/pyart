"""
pyart.io.sigmet
===============

Reading and writing of Sigmet (raw format) files

.. autosummary::
    :toctree: generated/

    read_sigmet
    SigmetFile
    convert_sigmet_data
    bin2_to_angle
    bin4_to_angle
    ymds_time_to_datetime
    _data_types_from_mask
    _is_bit_set
    _parse_ray_headers
    _unpack_structure
    _unpack_key
    _unpack_ingest_data_headers
    _unpack_ingest_data_header
    _unpack_raw_prod_bhdr
    _unpack_product_hdr
    _unpack_ingest_header

"""

from collections import OrderedDict
import struct
import datetime

import numpy as np

from .common import get_metadata, make_time_unit_str
from .radar import Radar

SPEED_OF_LIGHT = 299793000.0


def read_sigmet(filename, field_names=None, field_metadata=None,
                sigmet_field_names=False, debug=False):
    """
    Read a Sigmet (IRIS) product file.

    Parameters
    ----------
    filename : str
        Name of Sigmet (IRIS) product file to read.
    field_names : dict, optional
        Dictionary mapping Sigmet data type names to radar field names. If a
        data type found in the file does not appear in this dictionary it will
        not be placed in the radar.fields dictionary.  If None (default) a
        standard mapping will be used (SIGMET_TO_STANDARD)
    field_metadata : dict of dicts, optional
        Dictionary of dictionaries to retrieve field metadata from, if a mapped
        field does not appear in this list standard metadata will be used, or
        none if these is no metadata.  None (default) uses a blank dictionary.
    sigmet_field_names : bool, optional
        True to use the Sigmet data type names for the field names. If this
        case the field_metadata and field_names parameters are ignored and the
        returned radar object has a fields attribute filled with the sigmet
        data type names with no metadata.
    debug : bool, optional

    Returns
    -------
    radar : Radar
        Radar object

    """
    # parse parameters
    if field_names is None:
        field_names = SIGMET_TO_STANDARD
    if field_metadata is None:
        field_metadata = {}

    # open the file, read data
    sigmetfile = SigmetFile(filename, debug=debug)
    sigmet_data, sigmet_metadata = sigmetfile.read_data()
    sigmetfile.close()

    ingest_config = sigmetfile.ingest_header['ingest_configuration']
    task_config = sigmetfile.ingest_header['task_configuration']
    first_data_type = sigmetfile.data_type_names[0]
    nsweeps, nrays, nbins = sigmet_data[first_data_type].shape

    # time
    time = get_metadata('time')

    dts = [ymds_time_to_datetime(d['sweep_start_time'])
           for d in sigmetfile.ingest_data_headers[first_data_type]]

    tdata = sigmet_metadata[first_data_type]['time'].astype('float64')
    for i, dt in enumerate(dts):
        tdata[i] += (dt - dts[0]).total_seconds()
    time['data'] = tdata.flatten()
    time['units'] = make_time_unit_str(dts[0])

    # _range
    _range = get_metadata('range')
    range_info = task_config['task_range_info']
    gate_0 = range_info['first_bin_range'] / 100.       # meters
    gate_nbin = range_info['last_bin_range'] / 100.     # meters
    gate_size = round((gate_nbin - gate_0) / (nbins))
    _range['data'] = gate_0 + gate_size * np.arange(nbins, dtype='float32')
    _range['meters_to_center_of_first_gate'] = np.array([gate_0],
                                                        dtype='float32')
    _range['meters_between_gates'] = np.array([gate_size], dtype='float32')

    # fields
    fields = {}
    for data_type_name, fdata in sigmet_data.iteritems():
        if sigmet_field_names:
            fields[data_type_name] = {'data': fdata.reshape(-1, nbins)}
        elif data_type_name in field_names:

            field_name = field_names[data_type_name]

            if field_name in field_metadata:
                field_dic = field_metadata[field_name].copy()
            else:
                field_dic = get_metadata(field_name)

            field_dic['data'] = fdata.reshape(-1, nbins)
            field_dic['_FillValue'] = -9999.0
            fields[field_name] = field_dic

    # metadata
    metadata = {'title': '', 'institution': '', 'references': '',
                'source': '', 'history': '', 'comment': ''}
    metadata['original_container'] = 'sigmet'
    metadata['instrument_name'] = ingest_config['site_name'].strip()

    # scan_type
    if task_config['task_scan_info']['antenna_scan_mode'] == 2:
        scan_type = 'rhi'
    else:
        scan_type = 'ppi'

    # latitude
    latitude = get_metadata('latitude')
    lat = bin4_to_angle(ingest_config['latitude_radar'])
    if lat > 180.0:
        lat -= 360.0
    latitude['data'] = np.array([lat], dtype='float64')

    # longitude
    longitude = get_metadata('longitude')
    lon = bin4_to_angle(ingest_config['longitude_radar'])
    if lon > 180.0:
        lon -= 360.0
    longitude['data'] = np.array([lon], dtype='float64')

    # altitude
    altitude = get_metadata('altitude')
    alt = sigmetfile.product_hdr['product_end']['ground_height']
    altitude['data'] = np.array([alt], dtype='float64')

    # sweep_number
    sweep_number = get_metadata('sweep_number')
    sweep_number['data'] = np.arange(nsweeps, dtype='int32')

    # sweep_mode
    sweep_mode = get_metadata('sweep_mode')
    if scan_type == 'ppi':
        sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])
    else:
        sweep_mode['data'] = np.array(nsweeps * ['rhi'])

    # fixed_angle
    fixed_angle = get_metadata('fixed_angle')
    fa = [d['fixed_angle'] for d in
          sigmetfile.ingest_data_headers[first_data_type]]
    fixed_angle['data'] = bin2_to_angle(np.array(fa, dtype='float32'))

    # sweep_start_ray_index
    sweep_start_ray_index = get_metadata('sweep_start_ray_index')
    sweep_start_ray_index['data'] = nrays * np.arange(nsweeps, dtype='int32')

    # sweep_end_ray_index
    sweep_end_ray_index = get_metadata('sweep_end_ray_index')
    seri = nrays * np.arange(1, nsweeps + 1, dtype='int32') - 1
    sweep_end_ray_index['data'] = seri

    # azimuth
    az0 = sigmet_metadata[first_data_type]['azimuth_0']
    az1 = sigmet_metadata[first_data_type]['azimuth_1']
    el0 = sigmet_metadata[first_data_type]['elevation_0']
    el1 = sigmet_metadata[first_data_type]['elevation_1']

    azimuth = get_metadata('azimuth')
    az = (az0 + az1) / 2.
    az[np.where(np.abs(az0 - az1) > 180.0)] += 180.
    az[az > 360.0] -= 360.
    azimuth['data'] = az.astype('float32').flatten()

    # elevation
    elevation = get_metadata('elevation')
    elevation['data'] = ((el0 + el1) / 2.).astype('float32').flatten()

    # instrument_parameters
    prt = get_metadata('prt')
    prt_mode = get_metadata('prt_mode')
    nyquist_velocity = get_metadata('nyquist_velocity')
    unambiguous_range = get_metadata('unambiguous_range')

    trays = nsweeps * nrays
    prt_value = 1. / sigmetfile.product_hdr['product_end']['prf']
    prt['data'] = prt_value * np.ones(trays, dtype='float32')

    ur_value = SPEED_OF_LIGHT * prt_value / 2.
    unambiguous_range['data'] = ur_value * np.ones(trays, dtype='float32')

    # TODO Multi PRF mode when
    # task_config['task_dsp_info']['multi_prf_flag'] != 0
    prt_mode['data'] = np.array(nsweeps * ['fixed'])

    wavelength_cm = sigmetfile.product_hdr['product_end']['wavelength']
    nv_value = wavelength_cm / (10000.0 * 4.0 * prt_value)
    nyquist_velocity['data'] = nv_value * np.ones(trays, dtype='float32')

    instrument_parameters = {'unambiguous_range': unambiguous_range,
                             'prt_mode': prt_mode, 'prt': prt,
                             'nyquist_velocity': nyquist_velocity}

    return Radar(
        time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)

# This dictionary maps sigmet data types -> radar field names
# Users can pass their own version of this dictionary to the read_sigmet
# function as the field_names parameter.
SIGMET_TO_STANDARD = {
    #'DBT': 'DBT',                       # (1) Total Power
    #'DBZ': 'DBZ',                       # (2) Reflectivity
    #'VEL': 'VEL',                       # (3) Velocity
    #'WIDTH': 'WIDTH',                   # (4) Width
    #'ZDR': 'ZDR',                       # (5) Differential reflectivity
    #'DBZC': 'DBZC',                     # (7) Corrected reflectivity
    #'DBT2': 'DBT2',                     # (8) Total Power
    'DBZ2': 'reflectivity_horizontal',  # (9) Reflectivity
    'VEL2': 'mean_doppler_velocity',    # (10) Velocity
    #'WIDTH2': 'WIDTH2',                 # (11) Width
    'ZDR2': 'diff_reflectivity',        # (12) Differential reflectivity
    #'RAINRATE2': 'RAINRATE2',           # (13) Rainfall rate
    #'KDP': 'KDP',                       # (14) KDP (differential phase)
    'KDP2': 'diff_phase',               # (15) KDP (differential phase)
    #'PHIDP': 'PHIDP',                   # (16) PhiDP (differential phase)
    #'VELC': 'VELC',                     # (17) Corrected velocity
    #'SQI': 'SQI',                       # (18) SQI
    #'RHOHV': 'RHOHV',                   # (19) RhoHV
    'RHOHV2': 'copol_coeff',            # (20) RhoHV
    'DBZC2': 'reflectivity_horizontal_filtered',    # (21) Corrected Reflec.
    #'VELC2': 'VELC2',                   # (21) Corrected Velocity
    'SQI2': 'norm_coherent_power',      # (23) SQI
    'PHIDP2': 'dp_phase_shift',         # (24) PhiDP (differential phase)
    #'LDRH': 'LDRH',                     # (25) LDR xmt H, rcv V
    #'LDRH2': 'LDRH2',                   # (26) LDR xmt H, rcv V
    #'LDRV': 'LDRV',                     # (27) LDR xmt V, rcv H
    #'LDRV2': 'LDRV2',                   # (28) LDR xmt V, rcv H
    #'HEIGHT': 'HEIGHT',                 # (32) Height (1/10 km)
    #'VIL2': 'VIL2',                     # (33) Linear Liquid
    #'RAW': 'RAW',                       # (34) Raw Data
    #'SHEAR': 'SHEAR',                   # (35) Wind Shear
    #'DIVERGE2': 'DIVERGE2',             # (36) Divergence
    #'FLIQUID2': 'FLIQUID2',             # (37) Floated liquid
    #'USER': 'USER',                     # (38) User type
    #'OTHER': 'OTHER',                   # (39) Unspecified
    #'DEFORM2': 'DEFORM2',               # (40) Deformation
    #'VVEL2': 'VVEL2',                   # (41) Vertical velocity
    #'HVEL2': 'HVEL2',                   # (42) Horizontal velocity
    #'HDIR2': 'HDIR2',                   # (43) Horizontal wind direction
    #'AXDIL2': 'AXDIL2',                 # (44) Axis of dilation
    #'TIME2': 'TIME2',                   # (45) Time in seconds
    #'RHOH': 'RHOH',                     # (46) Rho, xmt H, rcv V
    #'RHOH2': 'RHOH2',                   # (47) Rho, xmt H, rcv V
    #'RHOV': 'RHOV',                     # (48) Rho, xmt V, rcv H
    #'RHOV2': 'RHOV2',                   # (49) Rho, xmt V, rcv H
    #'PHIH': 'PHIH',                     # (50) Phi, xmt H, rcv V
    #'PHIH2': 'PHIH2',                   # (51) Phi, xmt H, rcv V
    #'PHIV': 'PHIV',                     # (52) Phi, xmt V, rcv H
    #'PHIV2': 'PHIV2',                   # (53) Phi, xmt V, rcv H
    #'USER2': 'USER2',                   # (54) User type
    #'HCLASS': 'HCLASS',                 # (55) Hydrometeor class
    #'HCLASS2': 'HCLASS2',               # (56) Hydrometeor class
    #'ZDRC': 'ZDRC',                     # (57) Corrected diff. refl.
    #'ZDRC2': 'ZDRC2'                    # (58) Corrected diff. refl.
}


RECORD_SIZE = 6144      # Raw product file blocked into 6144 byte records


class SigmetFile():
    """
    A class for accessing data from Sigmet (IRIS) product files.

    Parameters
    ----------
    filename : str
        Filename.

    Attributes
    ----------
    debug : bool
        Set True to print out debugging information, False otherwise.
    product_hdr : dict
        Product_hdr structure.
    ingest_header : dict
        Ingest_header structure.
    ingest_data_headers : list of dict
        Ingest_data_header structures for each data type.  Indexed by the
        data type name (str).  None when data has not yet been read.
    data_types : list
        List of data types (int) in the file.
    data_type_names : list
        List of data type names (stR) in the file.
    ndata_types : int
        Number of data types in the file.
    _fh : file
        Open file being read.
    _record_number : int
        Last record number read.
    _rbuf : ndarray of dtype 'int16'
        Record buffer, holds current working record.  None when no record has
        been loaded.
    _rbuf_pos : int
        Working location in rbuf.  None when no record has been loaded.
    _raw_product_bhdrs : list
        List of raw_product_bhdr structure dictionaries seperated by sweep.
        None when data has not yet been read.

    """

    def __init__(self, filename, debug=True):
        """ initalize the object. """

        self.debug = debug

        # open the file
        fh = open(filename, 'rb')

        # read the headers from the first 2 records.
        self.product_hdr = _unpack_product_hdr(fh.read(RECORD_SIZE))
        self.ingest_header = _unpack_ingest_header(fh.read(RECORD_SIZE))

        # determine data types contained in the file
        self.data_types = self._determine_data_types()
        self.ndata_types = len(self.data_types)
        self.data_type_names = [SIGMET_DATA_TYPES[i] for i in self.data_types]

        # set attributes
        self.ingest_data_headers = None
        self._fh = fh
        self._record_number = 2
        self._rbuf = None
        self._rbuf_pos = None
        self._raw_product_bhdrs = None

    def _determine_data_types(self):
        """ Determine the available data types in the file. """
        # determine the available fields
        task_config = self.ingest_header['task_configuration']
        task_dsp_info = task_config['task_dsp_info']
        word1 = task_dsp_info['current_data_type_mask']['mask_word_0']
        word2 = task_dsp_info['current_data_type_mask']['mask_word_1']
        return _data_types_from_mask(word1, word2)

    def close(self):
        """ Close the file. """
        self._fh.close()

    def read_data(self):
        """
        Read all data from the file.

        Returns
        -------
        data : dict of ndarrays
            Data arrays of shape=(nsweeps, nrays, nbins) for each data type.
            Indexed by data type name (str).
        metadata : dict of dicts
            Arrays of 'azimuth_0', 'azimuth_1', 'elevation_0', 'elevation_1',
            'nbins', and 'time' for each data type.  Indexed by data type name
            (str).

        """

        # determine size of data
        nsweeps = self.ingest_header['task_configuration'][
            'task_scan_info']['number_sweeps']
        nbins = self.product_hdr['product_end']['number_bins']
        nrays = self.ingest_header['ingest_configuration'][
            'number_rays_sweep']

        # create empty outputs
        shape = (nsweeps, nrays, nbins)
        data = dict([(name, np.ma.empty(shape, dtype='float32'))
                    for name in self.data_type_names])

        metadata = {}
        for name in self.data_type_names:
            header_dic = {
                'azimuth_0': np.empty((nsweeps, nrays), dtype='float32'),
                'elevation_0': np.empty((nsweeps, nrays), dtype='float32'),
                'azimuth_1': np.empty((nsweeps, nrays), dtype='float32'),
                'elevation_1': np.empty((nsweeps, nrays), dtype='float32'),
                'nbins': np.empty((nsweeps, nrays), dtype='int16'),
                'time': np.empty((nsweeps, nrays), dtype='uint16')}
            metadata[name] = header_dic

        self.ingest_data_headers = dict([(name, []) for name in
                                         self.data_type_names])

        self._raw_product_bhdrs = []

        # read in data sweep by sweep
        for i in xrange(nsweeps):
            ingest_data_hdrs, sweep_data, sweep_metadata = self._get_sweep()

            for j, name in enumerate(self.data_type_names):
                az0, el0, az1, el1, ray_nbins, ray_time = sweep_metadata[j]

                data[name][i] = sweep_data[j]
                metadata[name]['azimuth_0'][i] = az0
                metadata[name]['azimuth_1'][i] = az1
                metadata[name]['elevation_0'][i] = el0
                metadata[name]['elevation_1'][i] = el1
                metadata[name]['nbins'][i] = ray_nbins
                metadata[name]['time'][i] = ray_time
                self.ingest_data_headers[name].append(ingest_data_hdrs[j])

        return data, metadata

    def _get_sweep(self):
        """
        Get the data and metadata from the next sweep.

        Returns
        -------
        ingest_data_headers : list of dict
            List of ingest_data_header structures for each data type.
        sweep_data : list of arrays
            Sweep data for each data types in the order they appear in the
            file.
        sweep_metadata : list of tuples
            Sweep metadata for each data type in the same order as sweep_data.

        """

        # get the next record
        lead_record = self._fh.read(RECORD_SIZE)
        self._record_number += 1

        # unpack structures
        raw_prod_bhdr = _unpack_raw_prod_bhdr(lead_record)
        self._raw_product_bhdrs.append([raw_prod_bhdr])
        ingest_data_headers = _unpack_ingest_data_headers(
            lead_record, self.ndata_types)

        # determine size of data
        nray_data_types = [d['number_rays_file_actual']
                           for d in ingest_data_headers]
        nrays = sum(nray_data_types)    # total rays
        nbins = self.product_hdr['product_end']['number_bins']

        # prepare to read rays
        self._rbuf = np.fromstring(lead_record, dtype='int16')
        self._rbuf_pos = int((12 + 76 * self.ndata_types) / 2)
        raw_sweep_data = np.empty((nrays, nbins + 6), dtype='int16')

        # get the raw data ray-by-ray
        for ray_i in xrange(nrays):
            if self.debug:
                print "Reading ray: %i of %i" % (ray_i, nrays),
                print "self._rbuf_pos is", self._rbuf_pos
            raw_sweep_data[ray_i] = self._get_ray(nbins)

        # convert the data and parse the metadata
        sweep_data = []
        sweep_metadata = []
        for i, data_type in enumerate(self.data_types):
            sweep_data.append(convert_sigmet_data(
                data_type, raw_sweep_data[i::self.ndata_types, 6:]))
            sweep_metadata.append(_parse_ray_headers(
                raw_sweep_data[i::self.ndata_types, :6]))
        return ingest_data_headers, sweep_data, sweep_metadata

    def _get_ray(self, nbins):
        """
        Get the next ray, loading new records as needed.

        Parameters
        ----------
        nbins : int
            Number of bins in the ray.

        Returns
        -------
        out : array
            Raw data for the ray.

        """
        compression_code = self._rbuf[self._rbuf_pos]
        self._incr_rbuf_pos()
        out = np.empty((nbins + 6,), dtype='int16')
        out_pos = 0

        while compression_code != 1:

            if compression_code < 0:
                words = compression_code + 32768    # last 7 bits give size
                if self._rbuf_pos + words <= 3072:
                    # all compressed data is in the current record
                    out[out_pos:out_pos + words] = \
                        self._rbuf[self._rbuf_pos:self._rbuf_pos + words]
                    self._incr_rbuf_pos(words)
                    out_pos += words
                else:
                    # data is split between current and next record
                    out1 = self._rbuf[self._rbuf_pos:]
                    remain = words - (3072 - self._rbuf_pos)
                    self._load_record()
                    out2 = self._rbuf[self._rbuf_pos:self._rbuf_pos + remain]
                    out[out_pos:out_pos + words] = np.append(out1, out2)
                    self._incr_rbuf_pos(remain)
                    out_pos += words
            else:
                # add zeros to out
                out[out_pos: out_pos + compression_code] = 0
                out_pos += compression_code
            compression_code = self._rbuf[self._rbuf_pos]
            self._incr_rbuf_pos()

        return out

    def _incr_rbuf_pos(self, incr=1):
        """
        Increment the record buffer position, load a new record if needed.
        """
        self._rbuf_pos += incr
        if self._rbuf_pos >= 3072:
            self._load_record()

    def _load_record(self):
        """ Load the next record. """
        record = self._fh.read(RECORD_SIZE)
        self._record_number += 1
        if self.debug:
            print "Finished loading record:", self._record_number
        self._raw_product_bhdrs[-1].append(_unpack_raw_prod_bhdr(record))
        self._rbuf = np.fromstring(record, dtype='int16')
        self._rbuf_pos = 6

# functions used by the SigmetFile class


def _data_types_from_mask(word0, word1):
    """
    Return a list of the data types from the words in the data_type mask.
    """
    data_types = [i for i in range(32) if _is_bit_set(word0, i)]
    data_types += [i+32 for i in range(32) if _is_bit_set(word1, i)]
    return data_types


def _is_bit_set(number, bit):
    """ Return True if bit is set in number. """
    return number >> bit & 1 == 1


def _parse_ray_headers(ray_headers):
    """
    Parse the metadata from Sigmet ray headers.

    Parameters
    ----------
    ray_headers : array, shape=(..., 6)
        Ray headers to parse.

    Returns
    -------
    az0 : array
        Azimuth angles (in degrees) at beginning of the rays.
    el0 : array
        Elevation angles at the beginning of the rays.
    az1 : array
        Azimuth angles at the end of the rays.
    el1 : array
        Elevation angles at the end of the rays.
    nbins : array
        Number of bins in the rays.
    time : array
        Seconds since the start of the sweep for the rays.

    """
    az0 = bin2_to_angle(ray_headers.view('uint16')[..., 0])
    el0 = bin2_to_angle(ray_headers.view('uint16')[..., 1])
    az1 = bin2_to_angle(ray_headers.view('uint16')[..., 2])
    el1 = bin2_to_angle(ray_headers.view('uint16')[..., 3])
    nbins = ray_headers.view('int16')[..., 4]
    time = ray_headers.view('uint16')[..., 5]
    return (az0, el0, az1, el1, nbins, time)


###################
# format converts #
###################

SIGMET_DATA_TYPES = {
    1: 'DBT',
    2: 'DBZ',
    3: 'VEL',
    4: 'WIDTH',
    5: 'ZDR',
    7: 'DBZC',
    8: 'DBT2',
    9: 'DBZ2',
    10: 'VEL2',
    11: 'WIDTH2',
    12: 'ZDR2',
    13: 'RAINRATE2',
    14: 'KDP',
    15: 'KDP2',
    16: 'PHIDP',
    17: 'VELC',
    18: 'SQI',
    19: 'RHOHV',
    20: 'RHOHV2',
    21: 'DBZC2',
    22: 'VELC2',
    23: 'SQI2',
    24: 'PHIDP2',
    25: 'LDRH',
    26: 'LDRH2',
    27: 'LDRV',
    28: 'LDRV2',
    32: 'HEIGHT',
    33: 'VIL2',
    34: 'RAW',
    35: 'SHEAR',
    36: 'DIVERGE2',
    37: 'FLIQUID2',
    38: 'USER',
    39: 'OTHER',
    40: 'DEFORM2',
    41: 'VVEL2',
    42: 'HVEL2',
    43: 'HDIR2',
    44: 'AXDIL2',
    45: 'TIME2',
    46: 'RHOH',
    47: 'RHOH2',
    48: 'RHOV',
    49: 'RHOV2',
    50: 'PHIH',
    51: 'PHIH2',
    52: 'PHIV',
    53: 'PHIV2',
    54: 'USER2',
    55: 'HCLASS',
    56: 'HCLASS2',
    57: 'ZDRC',
    58: 'ZDRC2'
}


def convert_sigmet_data(data_type, data):
    """ Convert sigmet data. """
    out = np.ma.empty_like(data, dtype='float32')

    data_type_name = SIGMET_DATA_TYPES[data_type]

    like_dbt2 = [
        'DBT2',     # 2-byte Reflectivity Format, section 4.3.4
        'DBZ2',     # " "
        'KDP2',     # 2-byte KDP Format, section 4.3.13
        'LDRH2',    # 2-byte LDR Format, section 4.3.15
        'LDRV2',    # " "
        'VEL2',     # 2-byte Velocity Format, section 4.3.30
        'VELC2',    # 2-byte Unfolded Velocity Format, section 4.3.32
        'ZDR2',     # 2-byte ZDR Format, section 4.3.38
        'DBZC2',    # Corrected reflectivity, XXX not certain of format
        'ZDRC2',    # Corrected differential reflectivity, XXX not certain
    ]

    like_sqi2 = [
        'RHOV2',   # 2-byte Rho Format, section 4.3.22
        'RHOH2',   # " "
        'RHOHV2',  # 2-byte RhoHV Format, section 4.3.24
        'SQI2',    # 2-byte Signal Quality Index Format, section 4.3.27
    ]

    if data_type_name in like_dbt2:
        # value = (N - 32768) / 100.
        # 0 : no data available (mask)
        # 65535 Reserved for area not scanned in product file (nothing)
        out[:] = (data.view('uint16') - 32768.) / 100.
        out[data.view('uint16') == 0] = np.ma.masked

    elif data_type_name in like_sqi2:
        # value = (N - 1) / 65533
        # 0 : no data available (mask)
        # 65535 Area not scanned
        out[:] = (data.view('uint16') - 1.) / 65533.
        out[data.view('uint16') == 0] = np.ma.masked

    elif data_type_name == 'WIDTH2':
        # DB_WIDTH2, 11, Width (2 byte)
        # 2-byte Width Format, section 4.3.36
        out[:] = data.view('uint16') / 100.
        out[data.view('uint16') == 0] = np.ma.masked

    elif data_type_name == 'PHIDP2':
        # DB_PHIDP2, 24, PhiDP (Differential Phase) (2 byte)
        # 2-byte PhiDP format, section 4.3.19
        out[:] = 360. * (data.view('uint16') - 1.) / 65534.
        out[data.view('uint16') == 0] = np.ma.masked

    elif data_type_name == 'HCLASS2':
        # DB_HCLASS2, 56, Hydrometeor class (2 byte)
        # 2-byte HydroClass Format, section 4.3.9
        out[:] = data.view('uint16')

    else:
        # TODO implement conversions for additional formats.
        raise NotImplementedError

    out.set_fill_value(-9999.0)
    return out


def bin2_to_angle(bin2):
    """ Return an angle from Sigmet bin2 encoded value (or array). """
    return 360. * bin2 / 65536


def bin4_to_angle(bin4):
    """ Return an angle from Sigmet bin4 encoded value (or array). """
    return 360. * bin4 / 4294967296


def ymds_time_to_datetime(ymds):
    """ Return a datetime object from a Sigmet ymds_time dictionary. """
    dt = datetime.datetime(ymds['year'], ymds['month'], ymds['day'])
    delta = datetime.timedelta(seconds=ymds['seconds'])
    return dt + delta


#####################
# get/put functions #
#####################


def _unpack_structure(string, structure_dict):
    """ Unpack a structure """
    fmt = ''.join(structure_dict.values())
    l = struct.unpack(fmt, string)
    return dict(zip(structure_dict, l))


def _unpack_key(dic, key, structure_dic):
    """ Unpack a key. """
    dic[key] = _unpack_structure(dic[key], structure_dic)


def _unpack_ingest_data_headers(record, ndata_types):
    """
    Unpack one or more ingest_data_header from a record.

    Returns a list of dictionaries.

    """
    return [_unpack_ingest_data_header(record, i) for i in range(ndata_types)]


def _unpack_ingest_data_header(record, number):
    """
    Unpack a single ingest_data_header from record.

    """
    offset = 12 + 76 * number
    string = record[offset:offset + 76]
    idh = _unpack_structure(string, INGEST_DATA_HEADER)
    _unpack_key(idh, 'structure_header', STRUCTURE_HEADER)
    _unpack_key(idh, 'sweep_start_time', YMDS_TIME)
    if idh['structure_header']['structure_identifier'] != 24:
        raise ValueError
    return idh


def _unpack_raw_prod_bhdr(record):
    """ Return a dict with the unpacked raw_prod_bhdr from a record. """
    return _unpack_structure(record[:12], RAW_PROD_BHDR)


def _unpack_product_hdr(record):
    """
    Return a dict with the unpacked product_hdr from the first record.
    """

    # unpack the product_hdr structure from the first record
    product_hdr = _unpack_structure(record[:640], PRODUCT_HDR)

    # product_hdr substructure
    _unpack_key(product_hdr, 'structure_header', STRUCTURE_HEADER)
    _unpack_key(product_hdr, 'product_configuration',
                PRODUCT_CONFIGURATION)
    _unpack_key(product_hdr, 'product_end', PRODUCT_END)

    # product_config substructure
    product_config = product_hdr['product_configuration']
    _unpack_key(product_config, 'structure_header', STRUCTURE_HEADER)
    _unpack_key(product_config, 'generation_time', YMDS_TIME)
    _unpack_key(product_config, 'sweep_ingest_time', YMDS_TIME)
    _unpack_key(product_config, 'file_ingest_time', YMDS_TIME)
    _unpack_key(product_config, 'color_scale_def', COLOR_SCALE_DEF)

    # product_end substructure
    product_end = product_hdr['product_end']
    _unpack_key(product_end, 'ingest_time', YMDS_TIME)

    return product_hdr


def _unpack_ingest_header(record):
    """
    Return a dict with the unpacked ingest_header from the second record.
    """

    # unpack the ingest_header structure from the second_record
    ingest_header = _unpack_structure(record[:4884], INGEST_HEADER)

    # ingest_header substructure
    _unpack_key(ingest_header, 'structure_header', STRUCTURE_HEADER)
    _unpack_key(ingest_header, 'ingest_configuration',
                INGEST_CONFIGURATION)
    _unpack_key(ingest_header, 'task_configuration', TASK_CONFIGURATION)

    # ingest_configuration substructure
    ingest_configuration = ingest_header['ingest_configuration']
    _unpack_key(ingest_configuration, 'volume_scan_start_time', YMDS_TIME)

    # task_configuration substructure
    task_configuration = ingest_header['task_configuration']
    _unpack_key(task_configuration, 'structure_header', STRUCTURE_HEADER)
    _unpack_key(task_configuration, 'task_sched_info', TASK_SCHED_INFO)
    _unpack_key(task_configuration, 'task_dsp_info', TASK_DSP_INFO)
    _unpack_key(task_configuration, 'task_calib_info', TASK_CALIB_INFO)
    _unpack_key(task_configuration, 'task_range_info', TASK_RANGE_INFO)
    _unpack_key(task_configuration, 'task_scan_info', TASK_SCAN_INFO)
    _unpack_key(task_configuration, 'task_misc_info', TASK_MISC_INFO)
    _unpack_key(task_configuration, 'task_end_info', TASK_END_INFO)

    # task_dsp_info substructure
    task_dsp_info = task_configuration['task_dsp_info']
    _unpack_key(task_dsp_info, 'current_data_type_mask', DSP_DATA_MASK)
    _unpack_key(task_dsp_info, 'original_data_type_mask', DSP_DATA_MASK)
    _unpack_key(task_dsp_info, 'task_dsp_mode', TASK_DSP_MODE_BATCH)

    # task_scan_info substructure
    # TODO unpack task_scan_type_scan_info based on scan type
    #task_scan_info = task_configuration['task_scan_info']
    #    scan_type_struct =
    #_unpack_key(task_scan_info, 'task_scan_type_scan_info',
    #            scan_type_struct)

    # task_end_info substructure
    task_end_info = task_configuration['task_end_info']
    _unpack_key(task_end_info, 'task_data_time', YMDS_TIME)

    return ingest_header


##############
# structures #
##############

# scalar defitions, section 4.1, table 7, and corresponding
SINT1 = 'b'
UINT1 = 'B'
SINT2 = 'h'
UINT2 = 'H'
SINT4 = 'i'
UINT4 = 'I'
FLT4 = 'f'
FLT8 = 'd'
BIN1 = 'B'
BIN2 = 'H'      # these values need to be decoded with _bin2_to_angle
BIN4 = 'I'      # these values need to be decoded with _bin4_to_angle
MESSAGE = 'I'
UINT16_T = 'H'

# structures are taken from Vaisala PROGRAMMER'S MANUAL IRIS
# M21131EN-B

# Chapter 4 deals with data formats
# section 4.5.4 deals with the RAW Product format

# 640 bytes: product_hdr (section 4.2.25, page 47)
PRODUCT_HDR = OrderedDict([
    ('structure_header', '12s'),        # 12 bytes
    ('product_configuration', '320s'),  # 320 bytes
    ('product_end', '308s'),            # 308 bytes
])

# 12 bytes : structure_header (section 4.2.47)
STRUCTURE_HEADER = OrderedDict([
    ('structure_identifier', SINT2),
    ('format_version', SINT2),
    ('bytes_in_structure', SINT4),
    ('reserved', SINT2),
    ('flag', SINT2),
])

# 320 bytes: product_configuration (section 4.2.23, page 43) 320 bytes
PRODUCT_CONFIGURATION = OrderedDict([
    ('structure_header', '12s'),    # 12 bytes: structure_header
    ('product_type_code', UINT2),
    ('scheduling_code', UINT2),
    ('seconds_between_runs', SINT4),
    ('generation_time', '12s'),     # 12 bytes: ymds_time
    ('sweep_ingest_time', '12s'),   # 12 bytes: ymds_time
    ('file_ingest_time', '12s'),    # 12 bytes: ymds_time
    ('spare_0', '6s'),              # 6 bytes
    ('product_name', '12s'),
    ('task_name', '12s'),
    ('flag', UINT2),
    ('x_scale', SINT4),
    ('y_scale', SINT4),
    ('z_scale', SINT4),
    ('x_size', SINT4),
    ('y_size', SINT4),
    ('z_size', SINT4),
    ('x_location', SINT4),
    ('y_location', SINT4),
    ('z_location', SINT4),
    ('maximum_range', SINT4),
    ('data_type', UINT2),
    ('projection_name', '12s'),
    ('input_data_type', UINT2),
    ('projection_type', UINT1),
    ('spare_1', '1s'),              # 1 bytes
    ('radial_smoother', SINT2),
    ('times_run', SINT2),
    ('zr_constant', SINT4),
    ('zr_exponent', SINT4),
    ('x_smoother', SINT2),
    ('y_smoother', SINT2),
    ('product_specific_bytes', '80s'),  # 80 bytes:
    ('minor_task_suffix', '16s'),
    ('spare_2', '12s'),             # 12 bytes
    ('color_scale_def', '48s')      # 48 bytes: color_scale_def
])

# 12 bytes: ymds_time Structure (section 4.2.76, page 72)
YMDS_TIME = OrderedDict([
    ('seconds', SINT4),
    ('milliseconds', UINT2),    # milliseconds in lowest 10 bits,
    ('year', SINT2),
    ('month', SINT2),
    ('day', SINT2),
])

# 48 bytes: color_scale_def (section 4.2.5, page 34)
COLOR_SCALE_DEF = OrderedDict([
    ('iflags', UINT4),
    ('istart', SINT4),
    ('istep', SINT4),
    ('icolcnt', SINT2),
    ('iset_and_scale', UINT2),
    ('ilevel_seams', '32s')     # 32 bytes: UINT2[16]
])

# 308 bytes : product_end (section 4.2.24)
PRODUCT_END = OrderedDict([
    ('site_name', '16s'),
    ('iris_version_created', '8s'),
    ('ingest_iris_version', '8s'),
    ('ingest_time', '12s'),         # 12 bytes: ymds_time
    ('spare_0', '28s'),             # 28 bytes
    ('GMT_minute_offset_local', SINT2),
    ('ingest_hardware_name_', '16s'),
    ('ingest_site_name_', '16s'),
    ('GMT_minute_offset_standard', SINT2),
    ('latitude', BIN4),
    ('longitude', BIN4),
    ('ground_height', SINT2),
    ('radar_height', SINT2),
    ('prf', SINT4),
    ('pulse_width', SINT4),
    ('signal_processor_type', UINT2),
    ('trigger_rate', UINT2),
    ('samples_used', SINT2),
    ('clutter_filter', '12s'),
    ('number_linear_filter', UINT2),
    ('wavelength', SINT4),
    ('truncation_height', SINT4),
    ('first_bin_range', SINT4),
    ('last_bin_range', SINT4),
    ('number_bins', SINT4),
    ('flag', UINT2),
    ('number_ingest', SINT2),
    ('polarization', UINT2),
    ('horizontal_calibration_i0', SINT2),
    ('horizontal_calibration_noise', SINT2),
    ('horizontal_radar_constant', SINT2),
    ('reciever_bandwidth', UINT2),
    ('horizontal_current_noise', SINT2),
    ('vertical_current_noise', SINT2),
    ('ldr_offset', SINT2),
    ('zdr_offset', SINT2),
    ('tcf_cal_flags_1', UINT16_T),
    ('tcf_cal_flags_2', UINT16_T),
    ('spare_1', '18s'),             # 18 bytes
    ('standard_parallel_1', BIN4),
    ('standard_parallel_2', BIN4),
    ('earth_radius', UINT4),
    ('inverse_flatting', UINT4),
    ('fault_status', UINT4),
    ('input_mask', UINT4),
    ('number_log_filter', UINT2),
    ('cluttermap', UINT2),
    ('latitude_projection', BIN4),
    ('longitude_projection', BIN4),
    ('product_sequence_number', SINT2),
    ('spare_2', '32s'),             # 32 bytes
    ('melting_level', SINT2),
    ('radar_height_above_reference', SINT2),
    ('number_elements', SINT2),
    ('mean_wind_speed', UINT1),
    ('mean_wind_direction', BIN1),
    ('spare_3', '2s'),              # 2 bytes
    ('tz_name', '8s'),
    ('extended_product_header_offset', UINT4),
    ('spare_4', '4s'),              # 4 bytes
])

# 4884 bytes ingest_header Structure (section 4.2.16, page 40)
INGEST_HEADER = OrderedDict([
    ('structure_header', '12s'),        # 12 bytes: structure_header
    ('ingest_configuration', '480s'),   # 480 bytes: ingest_configuration
    ('task_configuration', '2612s'),    # 2612 bytes: task_configuration
    ('spare_0', '732s'),                # 732 bytes
    ('gparm', '128s'),                  # 128 bytes
    ('reserved', '920s'),               # 920 bytes
])

# 480 bytes ingest_configuration Structure (section 4.2.14, page 38)
INGEST_CONFIGURATION = OrderedDict([
    ('filename', '80s'),
    ('number_files', SINT2),
    ('number_sweeps_completed', SINT2),
    ('total_size', SINT4),
    ('volume_scan_start_time', '12s'),  # 12 bytes: ymds_time
    ('spare_0', '12s'),                 # 12 bytes
    ('ray_header_bytes', SINT2),
    ('extended_ray_header_bytes', SINT2),
    ('number_task_config_table', SINT2),
    ('playback_version', SINT2),
    ('spare_1', '4s'),                  # 4 bytes
    ('iris_version', '8s'),
    ('hardware_site', '16s'),
    ('gmt_offset_minutes_local', SINT2),
    ('site_name', '16s'),
    ('gmt_offset_minutes_standard', SINT2),
    ('latitude_radar', BIN4),
    ('longitude_radar', BIN4),
    ('height_site', SINT2),
    ('height_radar', SINT2),
    ('resolution_rays', UINT2),
    ('first_ray_index', UINT2),
    ('number_rays_sweep', UINT2),
    ('gparam_bytes', SINT2),
    ('altitude_radar', SINT4),
    ('velocity_east', SINT4),
    ('velocity_north', SINT4),
    ('velocity_up', SINT4),
    ('antenna_offset_starboard', SINT4),
    ('antenna_offset_bow', SINT4),
    ('antenna_offset_up', SINT4),
    ('fault_status', UINT4),
    ('melting_layer', SINT2),
    ('spare_2', '2s'),              # 2 bytes
    ('local_timezone', '8s'),
    ('flags', UINT4),
    ('configuration_name', '16s'),
    ('spare_3', '228s')
])

# 2612 bytes: task_configuration Structure (section 4.2.50, page 61)
TASK_CONFIGURATION = OrderedDict([
    ('structure_header', '12s'),    # 12 bytes: structure_header
    ('task_sched_info', '120s'),    # 120 bytes: task_sched_info
    ('task_dsp_info', '320s'),      # 320 bytes: task_dsp_info
    ('task_calib_info', '320s'),    # 320 bytes: task_calib_info
    ('task_range_info', '160s'),    # 160 bytes: task_range_info
    ('task_scan_info', '320s'),     # 320 bytes: task_scan_info
    ('task_misc_info', '320s'),     # 320 bytes: task_misc_info
    ('task_end_info', '320s'),      # 320 bytes: task_end_info
    ('comments', '720s'),
])

# 120 bytes: task_sched_info Structure (section 4.2.61, page 65)
TASK_SCHED_INFO = OrderedDict([
    ('start_time', SINT4),
    ('stop_time', SINT4),
    ('skip_time', SINT4),
    ('last_run_time', SINT4),
    ('time_used_last_run', SINT4),
    ('last_run_day', SINT4),
    ('flag', UINT2),
    ('spare_0', '94s'),
])

# 320 bytes: task_dsp_info Structure (section 4.2.51, page 61)
TASK_DSP_INFO = OrderedDict([
    ('major_mode', UINT2),
    ('dsp_type', UINT2),
    ('current_data_type_mask', '24s'),      # 24 bytes: dsp_data_mask
    ('original_data_type_mask', '24s'),     # 24 bytes: dsp_data_mask
    ('task_dsp_mode', '32s'),               # 32 bytes: task_dsp_mode_batch
    ('spare_0', '52s'),
    ('prf', SINT4),
    ('pulse_width', SINT4),
    ('multi_prf_flag', UINT2),
    ('dual_prf_delay', SINT2),
    ('agc_feedback_code', UINT2),
    ('sample_size', SINT2),
    ('gain_control_flag', UINT2),
    ('clutter_filter_name', '12s'),
    ('linear_filter_first_bin', UINT1),
    ('log_filter_first_bin', UINT1),
    ('attenuation', SINT2),
    ('gas_attenuation', UINT2),
    ('cluttermap_flag', UINT2),
    ('tranmitter_phase_sequence', UINT2),
    ('ray_header_mask', UINT4),
    ('playback_flag', UINT2),
    ('spare_1', '2s'),
    ('custom_ray_header_name', '16s'),
    ('spare_2', '120s')
])

# 24 bytes: dsp_data_mask Structure (section 4.2.7, page 36)
DSP_DATA_MASK = OrderedDict([
    ('mask_word_0', UINT4),
    ('extended_header_type', UINT4),
    ('mask_word_1', UINT4),
    ('mask_word_2', UINT4),
    ('mask_word_3', UINT4),
    ('mask_word_4', UINT4),
])

# 32 bytes: task_dsp_mode_batch (section 4.2.52, page 62)
TASK_DSP_MODE_BATCH = OrderedDict([
    ('low_prf_hz', UINT2),
    ('low_prf_factional', UINT2),
    ('low_prf_sample_size', SINT2),
    ('low_prf_range_averaging', SINT2),
    ('reflectivity_unfolding_threshold', SINT2),
    ('velocity_unfolding_threshold', SINT2),
    ('width_unfolding_threshold', SINT2),
    ('spare_0', '18s'),
])

# 320 bytes: task_calib_info Structure (section 4.2.49, page 59)
TASK_CALIB_INFO = OrderedDict([
    ('reflectivity_slope', SINT2),
    ('reflectivity_noise_threshold', SINT2),
    ('clutter_correction_threshold', SINT2),
    ('sqi_threshold', SINT2),
    ('power_threshold', SINT2),
    ('spare_0', '8s'),
    ('reflectivity_calibration', SINT2),
    ('uncorrected_reflectivity_threshold_flags', UINT2),
    ('corrected_reflectivity_threshold_flags', UINT2),
    ('velocity_threshold_flags', UINT2),
    ('width_threshold_flags', UINT2),
    ('zdr_threshold_flags', UINT2),
    ('spare_1', '6s'),
    ('flags', UINT2),
    ('spare_2', '2s'),
    ('ldr_bias', SINT2),
    ('zdr_bias', SINT2),
    ('nexrad_clutter_threshold', SINT2),
    ('nexrad_clutter_skip', UINT2),
    ('horizontal_i0_calibration', SINT2),
    ('vertical_i0_calibration', SINT2),
    ('horizontal_noise_calibration', SINT2),
    ('vertical_noise_calibration', SINT2),
    ('horizontal_radar_constant', SINT2),
    ('vertical_radar_constant', SINT2),
    ('reciever_bandwidth', UINT2),
    ('flags2', UINT16_T),
    ('spare_3', '256s'),
])

# 160 bytes: task_range_info Structure (section 4.2.58, page 64)
TASK_RANGE_INFO = OrderedDict([
    ('first_bin_range', SINT4),
    ('last_bin_range', SINT4),
    ('number_input_bins', SINT2),
    ('number_output_bins', SINT2),
    ('step_input_bins', SINT4),
    ('step_output_bins', SINT4),
    ('variable_range_bin_flag', UINT2),
    ('range_bin_averaging_flag', SINT2),
    ('spare_0', '136s'),
])

# 320 bytes: task_scan_info Structure (section 4.2.60, page 65)
TASK_SCAN_INFO = OrderedDict([
    ('antenna_scan_mode', UINT2),
    ('angular_resolution_desired', SINT2),
    ('spare_0', '2s'),
    ('number_sweeps', SINT2),
    ('task_scan_type_scan_info', '200s'),   # 200 bytes: task_foo_scan_info
    ('spare_1', '112s'),
])

# 200 bytes: task_rhi_scan_info Structure (section 4.2.59, page 64)
TASK_RHI_SCAN_INFO = OrderedDict([
    ('lower_elevation_limit', UINT2),
    ('upper_elevation_limit', UINT2),
    ('azimuth_list', '80s'),            # UINT2[40]
    ('spare_0', '115s'),
    ('start_first_sector_flag', 'c'),   # unknown type
])

# 200 bytes: task_ppi_scan_info (section 4.2.57, page 64)
TASK_PPI_SCAN_INFO = OrderedDict([
    ('left_azimuth_limit', BIN2),
    ('right_azimuth_limit', BIN2),
    ('elevation_list', '80s'),          # UINT2[40]
    ('spare_0', '115s'),
    ('start_first_section_flag', 'c'),  # unknown type
])

# 200 bytes: task_file_scan_info (section 4.2.54, page 63)
TASK_FILE_SCAN_INFO = OrderedDict([
    ('first_azimuth', UINT2),
    ('first_elevation', UINT2),
    ('filename', '12s'),
    ('spare_0', '184s'),
])

# 200 bytes: task_manual_scan_info (section 4.2.55, page 63)
TASK_MANUAL_SCAN_INFO = OrderedDict([
    ('flags', UINT2),
    ('spare_0', '198s'),
])

# 320 bytes: task_misc_info Structure (section 4.2.55, page 63)
TASK_MISC_INFO = OrderedDict([
    ('wavelength', SINT4),
    ('tr_serial_number', '16s'),
    ('transmit_power', SINT4),
    ('flags', UINT2),
    ('polarization_type', UINT2),
    ('trucation_height', SINT4),
    ('spare_0', '18s'),             # 18 bytes
    ('spare_1', '12s'),             # 12 bytes
    ('comment_bytes', SINT2),
    ('horizontal_beamwidth', BIN4),
    ('vertical_beamwidth', BIN4),
    ('customer_storage', '40s'),    # 40 bytes, uint4[10]
    ('spare_2', '208s'),
])

# 320 bytes: task_end_info Structure (section 4.2.53, page 62)
TASK_END_INFO = OrderedDict([
    ('task_major_number', SINT2),
    ('task_minor_number', SINT2),
    ('task_configuration_file_name', '12s'),
    ('task_description', '80s'),
    ('number_tasks', SINT4),
    ('task_state', UINT2),
    ('spare_0', '2s'),
    ('task_data_time', '12s'),      # 12 bytes: ymds_time
    ('spare_1', '204s'),
])


# 12 bytes raw_prod_bhdr structure (section 4.2.30, page 50)
RAW_PROD_BHDR = OrderedDict([
    ('record_number', SINT2),
    ('sweep_number', SINT2),
    ('first_ray_offset', SINT2),
    ('ray_number', SINT2),
    ('flags', UINT2),
    ('spare_0', '2s'),
])

# 76 bytes ingest_data_header (section 4.2.15, pages 40)
INGEST_DATA_HEADER = OrderedDict([
    ('structure_header', '12s'),    # 12 bytes: structure_header
    ('sweep_start_time', '12s'),    # 12 bytes: ymds_time
    ('sweep_number', SINT2),
    ('number_rays_sweep', SINT2),
    ('first_ray_index', SINT2),
    ('number_rays_file_expected', SINT2),
    ('number_rays_file_actual', SINT2),
    ('fixed_angle', BIN2),
    ('bit_per_bin', SINT2),
    ('data_type', UINT2),
    ('spare_0', '36s')      # 36 bytes
])
