"""
pyart.io._sigmetfile
====================

A class and supporting functions for reading Sigmet (raw format) files.

.. autosummary::
    :toctree: generated/

    SigmetFile
    convert_sigmet_data
    bin2_to_angle
    bin4_to_angle
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
import struct
import datetime
import warnings

import numpy as np
cimport numpy as np
cimport cython

RECORD_SIZE = 6144      # Raw product file blocked into 6144 byte records


cdef class SigmetFile:
    """
    A class for accessing data from Sigmet (IRIS) product files.

    Parameters
    ----------
    filename : str
        Filename or file-like object.

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
    _raw_product_bhdrs : list
        List of raw_product_bhdr structure dictionaries seperated by sweep.
        None when data has not yet been read.

    """
    cdef public debug, product_hdr, ingest_header, ingest_data_headers, \
        data_types, data_type_names, ndata_types,
    cdef public _fh, _raw_product_bhdrs

    cdef np.ndarray _rbuf
    cdef np.int16_t * _rbuf_p   # hack for fast indexing of _rbuf
    cdef public int _rbuf_pos, _record_number

    def __init__(self, filename, debug=False):
        """ initalize the object. """

        self.debug = debug

        # open the file
        if hasattr(filename, 'read'):
            fh = filename
        else:
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
        self._raw_product_bhdrs = []

    def _determine_data_types(self):
        """ Determine the available data types in the file. """
        # determine the available fields
        task_config = self.ingest_header['task_configuration']
        task_dsp_info = task_config['task_dsp_info']
        word0 = task_dsp_info['current_data_type_mask']['mask_word_0']
        word1 = task_dsp_info['current_data_type_mask']['mask_word_1']
        word2 = task_dsp_info['current_data_type_mask']['mask_word_2']
        word3 = task_dsp_info['current_data_type_mask']['mask_word_3']
        return _data_types_from_mask(word0, word1, word2, word3)

    def close(self):
        """ Close the file. """
        self._fh.close()

    def read_data(self, full_xhdr=False):
        """
        Read all data from the file.

        Parameters
        ----------
        full_xhdr : bool
            True to return the full extended headers if they exist padded with
            ones.  False will return a length 1 extended header converted to
            int32.  This is useful when the file contains a customer specified
            extended header (for example aircraft radar).

        Returns
        -------
        data : dict of ndarrays
            Data arrays of shape=(nsweeps, nrays, nbins) for each data type.
            Indexed by data type name (str).
        metadata : dict of dicts
            Arrays of 'azimuth_0', 'azimuth_1', 'elevation_0', 'elevation_1',
            'nbins', and 'time' for each data type.  Indexed by data type name
            (str).  Rays which were not collected are marked with a value of
            -1 in the 'nbins' array.

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
        if 'XHDR' in self.data_type_names:
            if full_xhdr:
                data['XHDR'] = np.ones(shape, dtype='int16')
            else:
                data['XHDR'] = np.ones((nsweeps, nrays, 1), dtype='int32')

        metadata = {}
        for name in self.data_type_names:
            header_dic = {
                'azimuth_0': np.empty((nsweeps, nrays), dtype='float32'),
                'elevation_0': np.empty((nsweeps, nrays), dtype='float32'),
                'azimuth_1': np.empty((nsweeps, nrays), dtype='float32'),
                'elevation_1': np.empty((nsweeps, nrays), dtype='float32'),
                'nbins': np.empty((nsweeps, nrays), dtype='int16'),
                'time': np.empty((nsweeps, nrays), dtype='uint16'),
                'prf_flag': np.empty((nsweeps, nrays), dtype='int16')}
            metadata[name] = header_dic

        self.ingest_data_headers = dict([(name, []) for name in
                                         self.data_type_names])

        self._raw_product_bhdrs = []

        # read in data sweep by sweep
        for i in xrange(nsweeps):
            ingest_data_hdrs, sweep_data, sweep_metadata = self._get_sweep(
                full_xhdr=full_xhdr)

            # check for a truncated file, return sweep(s) read up until error
            if ingest_data_hdrs is None:

                mess = ('File truncated or corrupt, %i of %i sweeps read' %
                        (i, nsweeps))
                warnings.warn(mess)

                for name in self.data_type_names:
                    data[name] = data[name][:i]
                    for k in metadata[name]:
                        metadata[name][k] = metadata[name][k][:i]
                return data, metadata

            for j, name in enumerate(self.data_type_names):
                temp = sweep_metadata[j]
                (az0, el0, az1, el1, ray_nbins, ray_time, prf_flag) = temp

                data[name][i] = sweep_data[j]
                metadata[name]['azimuth_0'][i] = az0
                metadata[name]['azimuth_1'][i] = az1
                metadata[name]['elevation_0'][i] = el0
                metadata[name]['elevation_1'][i] = el1
                metadata[name]['nbins'][i] = ray_nbins
                metadata[name]['time'][i] = ray_time
                metadata[name]['prf_flag'][i] = prf_flag
                self.ingest_data_headers[name].append(ingest_data_hdrs[j])

        # scale 1-byte velocity by the Nyquist (section 4.3.29)
        # this conversion is kept in this method so that the
        # product_hdr does not need to be accessed at lower abstraction
        # layers.
        if 'VEL' in self.data_type_names:
            wavelength_cm = self.product_hdr['product_end']['wavelength']
            prt_value = 1. / self.product_hdr['product_end']['prf']
            task_config = self.ingest_header['task_configuration']
            multi_prf_flag = task_config['task_dsp_info']['multi_prf_flag']
            if multi_prf_flag > 3 or multi_prf_flag < 0:
                multiplier = 1  # multiplier not defined in IRIS manual
            else:
                multiplier = [1, 2, 3, 4][multi_prf_flag]
            nyquist = wavelength_cm / (10000.0 * 4.0 * prt_value) * multiplier
            data['VEL'] *= nyquist
        # scale 1-byte width by the Nyquist
        if 'WIDTH' in self.data_type_names:
            # The IRIS Programmer's Manual indicates 1-byte width format data
            # should be scaled by the unambiguous velocity, twice the nyquist,
            # (section 4.3.35) but both RSL and RadX scale this data by the
            # nyquist.  Therefore to agree with these two packages the width
            # is scaled by the nyquist.
            wavelength_cm = self.product_hdr['product_end']['wavelength']
            prt_value = 1. / self.product_hdr['product_end']['prf']
            nyquist = wavelength_cm / (10000.0 * 4.0 * prt_value)
            data['WIDTH'] *= nyquist

        return data, metadata

    def _get_sweep(self, full_xhdr=False, raw_data=False):
        """
        Get the data and metadata from the next sweep.

        If the file ends early None is returned for all values.

        Parameters
        ----------
        full_xhdr : bool
            True to return the full extended headers if they exist padded with
            ones.  False will return a length 1 extended header converted to
            int32.  This is useful when the file contains a customer specified
            extended header (for example aircraft radar).
        raw_data : bool, optional
            True to return the raw_data for the given sweep, False to
            convert the data to floating point representation.

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

        # check if the file ended early, if so return Nones
        if len(lead_record) != RECORD_SIZE:
            return None, None, None

        # unpack structures
        raw_prod_bhdr = _unpack_raw_prod_bhdr(lead_record)
        self._raw_product_bhdrs.append([raw_prod_bhdr])
        ingest_data_headers = _unpack_ingest_data_headers(
            lead_record, self.ndata_types)
        if ingest_data_headers is None:
            return None, None, None

        # determine size of data
        nray_data_types = [d['number_rays_file_expected']
                           for d in ingest_data_headers]
        nrays = sum(nray_data_types)    # total rays
        nbins = self.product_hdr['product_end']['number_bins']

        # prepare to read rays
        self._rbuf = np.frombuffer(lead_record, dtype='int16')
        self._rbuf_p = <np.int16_t*>self._rbuf.data
        self._rbuf_pos = int((12 + 76 * self.ndata_types) / 2) - 1
        # set data initially to ones so that missing data can be better
        # seen when debugging
        raw_sweep_data = np.ones((nrays, nbins + 6), dtype='int16')

        # get the raw data ray-by-ray
        for ray_i in xrange(nrays):
            if self.debug:
                print("Reading ray: %i of %i" % (ray_i, nrays), end='')
                print("self._rbuf_pos is", self._rbuf_pos)
            if self._get_ray(nbins, raw_sweep_data[ray_i]):
                return None, None, None

        # return raw data if requested
        if raw_data:
            return ingest_data_headers, raw_sweep_data

        # convert the data and parse the metadata
        sweep_data = []
        sweep_metadata = []
        for i, data_type in enumerate(self.data_types):
            if data_type == 0 and full_xhdr:
                sweep_data.append(raw_sweep_data[i::self.ndata_types, 6:])
            else:
                sweep_data.append(convert_sigmet_data(
                    data_type, raw_sweep_data[i::self.ndata_types, 6:],
                    raw_sweep_data[i::self.ndata_types, 4]))
            sweep_metadata.append(_parse_ray_headers(
                raw_sweep_data[i::self.ndata_types, :6]))
        return ingest_data_headers, sweep_data, sweep_metadata

    @cython.wraparound(False)
    cdef int _get_ray(self, int nbins, np.ndarray[np.int16_t, ndim=1] out):
        """
        Get the next ray, loading new records as needed.

        Parameters
        ----------
        nbins : int
            Number of bins in the ray.
        out : ndarray
            Array to load ray data into.

        Returns
        -------
        status : int
            0 on success, -1 if failed.

        """
        cdef short compression_code
        cdef int words, remain, out_pos, first_end, i

        if self._incr_rbuf_pos():
            return -1   # failed read
        compression_code = self._rbuf_p[self._rbuf_pos]
        out_pos = 0

        if compression_code == 1:
            # mark ray as missing by setting numbers of bins to -1
            out[4] = -1
            return 0

        while compression_code != 1:

            if self._incr_rbuf_pos():
                return -1   # failed read
            if compression_code < 0:
                words = compression_code + 32768    # last 7 bits give size
                if self._rbuf_pos + words <= 3072:
                    # all compressed data is in the current record
                    for i in range(words):
                        out[out_pos + i] = self._rbuf_p[self._rbuf_pos + i]
                    if self._incr_rbuf_pos(words):
                        return -1   # failed read
                    out_pos += words
                else:
                    # data is split between current and next record
                    # store data from current record
                    remain = words - (3072 - self._rbuf_pos)
                    first_end = out_pos + words - remain
                    for i in range(first_end - out_pos):
                        out[out_pos + i] = self._rbuf_p[self._rbuf_pos + i]

                    # read data from next record and store
                    self._load_record()
                    for i in range(out_pos + words - first_end):
                        out[first_end + i] = self._rbuf_p[self._rbuf_pos + i]

                    if self._incr_rbuf_pos(remain):
                        return -1   # failed read
                    out_pos += words
            else:
                # add zeros to out
                if compression_code + out_pos > nbins + 6:
                    return -1   # file is corrupt
                for i in range(compression_code):
                    out[out_pos + i] = 0
                out_pos += compression_code
            compression_code = self._rbuf_p[self._rbuf_pos]

        return 0

    cdef int _incr_rbuf_pos(self, int incr=1):
        """
        Increment the record buffer position, load a new record if needed.
        """
        self._rbuf_pos += incr
        if self._rbuf_pos >= 3072:
            if self._load_record():
                return -1   # failed read
        return 0

    cdef int _load_record(self):
        """ Load the next record. returns -1 on fail, 0 if success. """
        record = self._fh.read(RECORD_SIZE)
        if len(record) != RECORD_SIZE:
            return -1   # failed read
        self._record_number += 1
        if self.debug:
            print("Finished loading record:", self._record_number)
        self._raw_product_bhdrs[-1].append(_unpack_raw_prod_bhdr(record))
        self._rbuf = np.frombuffer(record, dtype='int16')
        self._rbuf_pos = 6
        self._rbuf_p = <np.int16_t*>self._rbuf.data
        return 0

# functions used by the SigmetFile class


def _data_types_from_mask(word0, word1, word2, word3):
    """
    Return a list of the data types from the words in the data_type mask.
    """
    data_types = [i for i in range(32) if _is_bit_set(word0, i)]
    data_types += [i+32 for i in range(32) if _is_bit_set(word1, i)]
    data_types += [i+64 for i in range(32) if _is_bit_set(word2, i)]
    data_types += [i+96 for i in range(32) if _is_bit_set(word3, i)]
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
    prf_flag : array
        Numerical indication of what PRF was used, 0 for high, 1 for low.
        Not applicable if dual-PRF is not used during collection.

    """
    az0 = bin2_to_angle(ray_headers.view('uint16')[..., 0])
    el0 = bin2_to_angle(ray_headers.view('uint16')[..., 1])
    az1 = bin2_to_angle(ray_headers.view('uint16')[..., 2])
    el1 = bin2_to_angle(ray_headers.view('uint16')[..., 3])
    nbins = ray_headers.view('int16')[..., 4]
    time = ray_headers.view('uint16')[..., 5]
    prf_flag = np.mod(ray_headers.view('int16')[..., 0], 2)
    return (az0, el0, az1, el1, nbins, time, prf_flag)


###################
# format converts #
###################

# Data type constants, table 13, section 4.8
SIGMET_DATA_TYPES = {
    0: 'XHDR',
    1: 'DBT',
    2: 'DBZ',
    3: 'VEL',
    4: 'WIDTH',
    5: 'ZDR',
    6: 'UNKNOWN_6',     # Not known
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
    29: 'UNKNOWN_29',   # Not known
    30: 'UNKNOWN_30',   # Not known
    31: 'UNKNOWN_31',   # Not known
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
    58: 'ZDRC2',
    59: 'TEMPERATURE16',
    60: 'VIR16',
    61: 'DBTV8',
    62: 'DBTV16',
    63: 'DBZV8',
    64: 'DBZV16',
    65: 'SNR8',
    66: 'SNR16',
    67: 'ALBEDO8',
    68: 'ALBEDO16',
    69: 'VILD16',
    70: 'TURB16',
    71: 'DBTE8',
    72: 'DBTE16',       # Total Power Enhanced
    73: 'DBZE8',
    74: 'DBZE16',       # Clutter Corrected Reflectivity Enhanced
    # Uknown fields, do not know internal names, some may be user defined.
    75: 'UNKNOWN_75',
    76: 'UNKNOWN_76',
    77: 'UNKNOWN_77',
    78: 'UNKNOWN_78',
    79: 'UNKNOWN_79',
    80: 'UNKNOWN_80',
    81: 'UNKNOWN_81',
    82: 'UNKNOWN_82',
    83: 'UNKNOWN_83',
    84: 'UNKNOWN_84',
    85: 'UNKNOWN_85',
    86: 'UNKNOWN_86',
    87: 'UNKNOWN_87',
    88: 'UNKNOWN_88',
    89: 'UNKNOWN_89',
    90: 'UNKNOWN_90',
    91: 'UNKNOWN_91',
    92: 'UNKNOWN_92',
    93: 'UNKNOWN_93',
    94: 'UNKNOWN_94',
    95: 'UNKNOWN_95',
    96: 'UNKNOWN_96',
    97: 'UNKNOWN_97',
    98: 'UNKNOWN_98',
    99: 'UNKNOWN_99',
    100: 'UNKNOWN_100',
    101: 'UNKNOWN_101',
    102: 'UNKNOWN_102',
    103: 'UNKNOWN_103',
    104: 'UNKNOWN_104',
    105: 'UNKNOWN_105',
    106: 'UNKNOWN_106',
    107: 'UNKNOWN_107',
    108: 'UNKNOWN_108',
    109: 'UNKNOWN_109',
    110: 'UNKNOWN_110',
    111: 'UNKNOWN_110',
    112: 'UNKNOWN_112',
    113: 'UNKNOWN_113',
    114: 'UNKNOWN_114',
    115: 'UNKNOWN_115',
    116: 'UNKNOWN_116',
    117: 'UNKNOWN_117',
    118: 'UNKNOWN_118',
    119: 'UNKNOWN_119',
    120: 'UNKNOWN_120',
    121: 'UNKNOWN_121',
    122: 'UNKNOWN_122',
    123: 'UNKNOWN_123',
    124: 'UNKNOWN_124',
    125: 'UNKNOWN_125',
    126: 'UNKNOWN_126',
    127: 'UNKNOWN_127',
}   # there may be more field, add as needed


# This function takes a majority of the time when reading data from a sigmet
# file. Rewriting the convertions/masking in Cython does not seem to improved
# performance likely since most of the routines are already vectorized.
def convert_sigmet_data(data_type, data, nbins):
    """ Convert sigmet data. """
    out = np.empty_like(data, dtype='float32')
    mask = np.zeros_like(data, dtype=np.bool8)

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
        'DBTV16',   # Total V power, 2-byte
        'DBZV16',   # Clutter corrected V reflectivity, 2-byte
        'SNR16',    # Signal to noise ratio, 2-byte
        'DBTE16',   # Total Power Enhanced, 2-byte
        'DBZE16',   # Clutter corrected reflectivity enhanced, 2-byte
    ]

    like_sqi = [
        'RHOH',     # 1-byte Rho Format, section 4.3.21
        'RHOV',     # " "
        'RHOHV',    # 1-byte RhoHV Format, section 4.3.23
        'SQI',      # 1-byte Signal Quality Index Format, section 4.3.26
    ]

    like_sqi2 = [
        'RHOV2',    # 2-byte Rho Format, section 4.3.22
        'RHOH2',    # " "
        'RHOHV2',   # 2-byte RhoHV Format, section 4.3.24
        'SQI2',     # 2-byte Signal Quality Index Format, section 4.3.27
    ]

    like_dbt = [
        'DBT',      # 1-bytes Reflectivity Format, section 4.3.3
        'DBZ',      # " "
        'DBTV8',    # Total V power, 1-byte
        'DBZV8',    # Clutter corrected V reflectivity, 1-byte
        'SNR8',     # Signal to noise ratio, 1-byte
        'DBTE8',    # Total power enhanced, 1-byte
        'DBZE8',    # Clutter corrected reflectivity enhanced, 1-byte
    ]

    if data_type_name in like_dbt2:
        # value = (N - 32768) / 100.
        # 0 : no data available (mask)
        # 65535 Reserved for area not scanned in product file (nothing)
        out[:] = (data.view('uint16') - 32768.) / 100.
        mask[data.view('uint16') == 0] = True

    elif data_type_name in like_sqi2:
        # value = (N - 1) / 65533
        # 0 : no data available (mask)
        # 65535 Area not scanned
        out[:] = (data.view('uint16') - 1.) / 65533.
        mask[data.view('uint16') == 0] = True

    elif data_type_name == 'WIDTH2':
        # DB_WIDTH2, 11, Width (2 byte)
        # 2-byte Width Format, section 4.3.36
        out[:] = data.view('uint16') / 100.
        mask[data.view('uint16') == 0] = True

    elif data_type_name == 'PHIDP2':
        # DB_PHIDP2, 24, PhiDP (Differential Phase) (2 byte)
        # 2-byte PhiDP format, section 4.3.19
        out[:] = 360. * (data.view('uint16') - 1.) / 65534.
        mask[data.view('uint16') == 0] = True

    elif data_type_name == 'HCLASS2':
        # DB_HCLASS2, 56, Hydrometeor class (2 byte)
        # 2-byte HydroClass Format, section 4.3.9
        out[:] = data.view('uint16')

    elif data_type_name == 'XHDR':
        # Extended Headers, 0
        # extended_header_v0, _v1, _v2, section 4.2.8-4.2.10
        # Here we return an array with the times in milliseconds.
        return data[..., :2].copy().view('i4')

    # one byte data types
    elif data_type_name[-1] != '2':
        # make a view of left half of the data as uint8,
        # this is the actual ray data collected, the right half is blank.
        nrays, nbin = data.shape
        ndata = data.view('(2,) uint8').reshape(nrays, -1)[:, :nbin]

        if data_type_name in like_dbt:
            # DB_DBT, 1, Total Power (1 byte)
            # 1-byte Reflectivity Format, section 4.3.3
            out[:] = (ndata - 64.) / 2.
            mask[ndata == 0] = True

        elif data_type_name in like_sqi:
            # value = sqrt((N - 1) / 253)
            # 0 : no data available (mask)
            # 255 Area not scanned
            out[:] = np.sqrt((ndata - 1.) / 253.)
            mask[ndata == 0] = True
            mask[ndata == 255] = True

        elif data_type_name == 'VEL':
            # VEL, 3, Velocity (1 byte)
            # 1-byte Velocity Format, section 4.3.29
            # Note that this data should be multiplied by Nyquist,
            # this is done in the get_data method of the SigmetFile class.
            out[:] = (ndata - 128.) / 127.
            mask[ndata == 0] = True

        elif data_type_name == 'WIDTH':
            # WIDTH, 4, Width (1 byte)
            # 1-byte Width format, section 4.3.25
            # Note that this data should be multiplied by the unambiguous
            # velocity
            out[:] = ndata / 256.
            mask[ndata == 0] = True

        elif data_type_name == 'ZDR':
            # ZDR, 5, Differential reflectivity (1 byte)
            # 1-byte ZDR format, section 4.3.37
            out[:] = (ndata - 128.) / 16.
            mask[ndata == 0] = True

        elif data_type_name == 'KDP':
            # KDP, 14, KDP (Differential phase) (1 byte)
            # 1-byte KDP format, section 4.3.12
            # Note that this data should be divided by the wavelength in cm
            # as is the units are deg * cm / km

            # above 128 use positive value equation
            exp = np.power(600., (ndata[ndata > 128] - 129.) / 126.)
            out[ndata > 128] = 0.25 * exp
            # below 128, use negative value equation
            exp = np.power(600., (127. - ndata[ndata < 128]) / 126.)
            out[ndata < 128] = -0.25 * exp
            # equal to 128, zero
            out[ndata == 128] = 0

            mask[ndata == 0] = True
            mask[ndata == 255] = True

        elif data_type_name == 'PHIDP':
            # PHIDP, 16, PhiDP(Differential phase) (1 byte)
            # 1-byte PhiDP format, section 4.3.18
            out[:] = 180. * ((ndata - 1.) / 254.)
            mask[ndata == 0] = True
            mask[ndata == 255] = True

        elif data_type_name == "HCLASS":
            # HCLASS, 55, Hydrometeor class (1 byte)
            # 1-byte HydroClass format, section 4.3.8
            out[:] = ndata[:]
            mask[ndata == 0] = True     # No data available
            mask[ndata == 255] = True   # Area not scanned

        else:
            # TODO implement conversions for addition 1-byte formats
            warnings.warn('Unknown type: %s, returning raw data' % data_type)
            out[:] = np.ma.masked_array(data)
            return out
    else:
        # TODO implement conversions for additional formats.
        warnings.warn('Unknown type: %s, returning raw data' % data_type)
        out[:] = data
        return np.ma.masked_array(out)

    # mask any gates which are beyond the number of gates in that ray.
    _mask_gates_not_collected(mask.view(np.uint8), nbins)

    return np.ma.masked_array(out, mask=mask, fill_value=-9999.0,
                              shrink=False)


@cython.boundscheck(False)
cdef _mask_gates_not_collected(
        np.ndarray[np.uint8_t, ndim=2] mask,
        np.ndarray[np.int16_t, ndim=1] nbins):
    """ Add gates not collected (beyond nbin) to the mask. """
    cdef int i, j, nrays, nbin, full_nbins
    nrays = mask.shape[0]
    full_nbins = mask.shape[1]
    for i in range(nrays):
        nbin = nbins[i]
        for j in range(nbin, full_nbins):
            mask[i, j] = 1
    return


def bin2_to_angle(bin2):
    """ Return an angle from Sigmet bin2 encoded value (or array). """
    return 360. * bin2 / 65536


def bin4_to_angle(bin4):
    """ Return an angle from Sigmet bin4 encoded value (or array). """
    return 360. * bin4 / 4294967296


#####################
# get/put functions #
#####################


def _unpack_structure(string, structure):
    """ Unpack a structure """
    fmt = ''.join([i[1] for i in structure])
    l = struct.unpack(fmt, string)
    return dict(zip([i[0] for i in structure], l))


def _unpack_key(dic, key, structure):
    """ Unpack a key. """
    dic[key] = _unpack_structure(dic[key], structure)


def _unpack_ingest_data_headers(record, ndata_types):
    """
    Unpack one or more ingest_data_header from a record.

    Returns a list of dictionaries or None when an error occurs.

    """
    idh = [_unpack_ingest_data_header(record, i) for i in range(ndata_types)]
    if None in idh:
        return None
    else:
        return idh


def _unpack_ingest_data_header(record, number):
    """
    Unpack a single ingest_data_header from record.  Return None on error.
    """
    offset = 12 + 76 * number
    string = record[offset:offset + 76]
    idh = _unpack_structure(string, INGEST_DATA_HEADER)
    _unpack_key(idh, 'structure_header', STRUCTURE_HEADER)
    _unpack_key(idh, 'sweep_start_time', YMDS_TIME)
    if idh['structure_header']['structure_identifier'] != 24:
        return None
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
    # task_scan_info = task_configuration['task_scan_info']
    #    scan_type_struct =
    # _unpack_key(task_scan_info, 'task_scan_type_scan_info',
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
PRODUCT_HDR = (
    ('structure_header', '12s'),        # 12 bytes
    ('product_configuration', '320s'),  # 320 bytes
    ('product_end', '308s'),            # 308 bytes
)

# 12 bytes : structure_header (section 4.2.47)
STRUCTURE_HEADER = (
    ('structure_identifier', SINT2),
    ('format_version', SINT2),
    ('bytes_in_structure', SINT4),
    ('reserved', SINT2),
    ('flag', SINT2),
)

# 320 bytes: product_configuration (section 4.2.23, page 43) 320 bytes
PRODUCT_CONFIGURATION = (
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
)

# 12 bytes: ymds_time Structure (section 4.2.76, page 72)
YMDS_TIME = (
    ('seconds', SINT4),
    ('milliseconds', UINT2),    # milliseconds in lowest 10 bits,
    ('year', SINT2),
    ('month', SINT2),
    ('day', SINT2),
)

# 48 bytes: color_scale_def (section 4.2.5, page 34)
COLOR_SCALE_DEF = (
    ('iflags', UINT4),
    ('istart', SINT4),
    ('istep', SINT4),
    ('icolcnt', SINT2),
    ('iset_and_scale', UINT2),
    ('ilevel_seams', '32s')     # 32 bytes: UINT2[16]
)

# 308 bytes : product_end (section 4.2.24)
PRODUCT_END = (
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
)

# 4884 bytes ingest_header Structure (section 4.2.16, page 40)
INGEST_HEADER = (
    ('structure_header', '12s'),        # 12 bytes: structure_header
    ('ingest_configuration', '480s'),   # 480 bytes: ingest_configuration
    ('task_configuration', '2612s'),    # 2612 bytes: task_configuration
    ('spare_0', '732s'),                # 732 bytes
    ('gparm', '128s'),                  # 128 bytes
    ('reserved', '920s'),               # 920 bytes
)

# 480 bytes ingest_configuration Structure (section 4.2.14, page 38)
INGEST_CONFIGURATION = (
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
)

# 2612 bytes: task_configuration Structure (section 4.2.50, page 61)
TASK_CONFIGURATION = (
    ('structure_header', '12s'),    # 12 bytes: structure_header
    ('task_sched_info', '120s'),    # 120 bytes: task_sched_info
    ('task_dsp_info', '320s'),      # 320 bytes: task_dsp_info
    ('task_calib_info', '320s'),    # 320 bytes: task_calib_info
    ('task_range_info', '160s'),    # 160 bytes: task_range_info
    ('task_scan_info', '320s'),     # 320 bytes: task_scan_info
    ('task_misc_info', '320s'),     # 320 bytes: task_misc_info
    ('task_end_info', '320s'),      # 320 bytes: task_end_info
    ('comments', '720s'),
)

# 120 bytes: task_sched_info Structure (section 4.2.61, page 65)
TASK_SCHED_INFO = (
    ('start_time', SINT4),
    ('stop_time', SINT4),
    ('skip_time', SINT4),
    ('last_run_time', SINT4),
    ('time_used_last_run', SINT4),
    ('last_run_day', SINT4),
    ('flag', UINT2),
    ('spare_0', '94s'),
)

# 320 bytes: task_dsp_info Structure (section 4.2.51, page 61)
TASK_DSP_INFO = (
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
)

# 24 bytes: dsp_data_mask Structure (section 4.2.7, page 36)
DSP_DATA_MASK = (
    ('mask_word_0', UINT4),
    ('extended_header_type', UINT4),
    ('mask_word_1', UINT4),
    ('mask_word_2', UINT4),
    ('mask_word_3', UINT4),
    ('mask_word_4', UINT4),
)

# 32 bytes: task_dsp_mode_batch (section 4.2.52, page 62)
TASK_DSP_MODE_BATCH = (
    ('low_prf_hz', UINT2),
    ('low_prf_factional', UINT2),
    ('low_prf_sample_size', SINT2),
    ('low_prf_range_averaging', SINT2),
    ('reflectivity_unfolding_threshold', SINT2),
    ('velocity_unfolding_threshold', SINT2),
    ('width_unfolding_threshold', SINT2),
    ('spare_0', '18s'),
)

# 320 bytes: task_calib_info Structure (section 4.2.49, page 59)
TASK_CALIB_INFO = (
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
)

# 160 bytes: task_range_info Structure (section 4.2.58, page 64)
TASK_RANGE_INFO = (
    ('first_bin_range', SINT4),
    ('last_bin_range', SINT4),
    ('number_input_bins', SINT2),
    ('number_output_bins', SINT2),
    ('step_input_bins', SINT4),
    ('step_output_bins', SINT4),
    ('variable_range_bin_flag', UINT2),
    ('range_bin_averaging_flag', SINT2),
    ('spare_0', '136s'),
)

# 320 bytes: task_scan_info Structure (section 4.2.60, page 65)
TASK_SCAN_INFO = (
    ('antenna_scan_mode', UINT2),
    ('angular_resolution_desired', SINT2),
    ('spare_0', '2s'),
    ('number_sweeps', SINT2),
    ('task_scan_type_scan_info', '200s'),   # 200 bytes: task_foo_scan_info
    ('spare_1', '112s'),
)

# 200 bytes: task_rhi_scan_info Structure (section 4.2.59, page 64)
TASK_RHI_SCAN_INFO = (
    ('lower_elevation_limit', UINT2),
    ('upper_elevation_limit', UINT2),
    ('azimuth_list', '80s'),            # UINT2[40]
    ('spare_0', '115s'),
    ('start_first_sector_flag', 'c'),   # unknown type
)

# 200 bytes: task_ppi_scan_info (section 4.2.57, page 64)
TASK_PPI_SCAN_INFO = (
    ('left_azimuth_limit', BIN2),
    ('right_azimuth_limit', BIN2),
    ('elevation_list', '80s'),          # UINT2[40]
    ('spare_0', '115s'),
    ('start_first_section_flag', 'c'),  # unknown type
)

# 200 bytes: task_file_scan_info (section 4.2.54, page 63)
TASK_FILE_SCAN_INFO = (
    ('first_azimuth', UINT2),
    ('first_elevation', UINT2),
    ('filename', '12s'),
    ('spare_0', '184s'),
)

# 200 bytes: task_manual_scan_info (section 4.2.55, page 63)
TASK_MANUAL_SCAN_INFO = (
    ('flags', UINT2),
    ('spare_0', '198s'),
)

# 320 bytes: task_misc_info Structure (section 4.2.55, page 63)
TASK_MISC_INFO = (
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
)

# 320 bytes: task_end_info Structure (section 4.2.53, page 62)
TASK_END_INFO = (
    ('task_major_number', SINT2),
    ('task_minor_number', SINT2),
    ('task_configuration_file_name', '12s'),
    ('task_description', '80s'),
    ('number_tasks', SINT4),
    ('task_state', UINT2),
    ('spare_0', '2s'),
    ('task_data_time', '12s'),      # 12 bytes: ymds_time
    ('spare_1', '204s'),
)


# 12 bytes raw_prod_bhdr structure (section 4.2.30, page 50)
RAW_PROD_BHDR = (
    ('record_number', SINT2),
    ('sweep_number', SINT2),
    ('first_ray_offset', SINT2),
    ('ray_number', SINT2),
    ('flags', UINT2),
    ('spare_0', '2s'),
)

# 76 bytes ingest_data_header (section 4.2.15, pages 40)
INGEST_DATA_HEADER = (
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
)
