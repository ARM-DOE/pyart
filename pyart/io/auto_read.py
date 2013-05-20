"""
pyart.io.auto_read
==================

Automatic reading of radar files by detecting format.

.. autosummary::
    :toctree: generated/

    read
    determine_filetype

"""

from .mdv import read_mdv
from .rsl import read_rsl
from .netcdf import read_netcdf


def read(filename, callid='KABR'):
    """
    Read a radar file and return a radar object.

    Parameters
    ----------
    filename : str
        Name of radar file to read
    callid : str
        Four letter NEXRAD radar call id, only used if format is determined to
        be 'WSR88D'.  The default value will set the location of the radar to
        Aberdeen, SD.  The fields will still be correct.

    Returns
    -------
    radar : Radar
        Radar object.  A TypeError is raised if the format cannot be
        determined.

    """
    filetype = determine_filetype(filename)

    if filetype == "MDV":
        return read_mdv(filename)
    if filetype == "NETCDF3" or filetype == "NETCDF4":
        return read_netcdf(filename)
    if filetype == 'WSR88D':
        read_rsl(filename, callid=callid)
    rsl_formats = ['UF', 'HDF4', 'RSL', 'DORAD', 'SIGMET']
    if filetype in rsl_formats:
        return read_rsl(filename)

    raise TypeError('file format is cannot be determined.')


def determine_filetype(filename):
    """
    Return the filetype of a given file by examining the first few bytes.

    The following filetypes are detected:

    * 'MDV'
    * 'NETCDF3'
    * 'NETCDF4'
    * 'WSR88D'
    * 'UF'
    * 'HDF4'
    * 'RSL'
    * 'DORAD'
    * 'SIGMET'
    * 'UNKNOWN'

    Parameters
    ----------
    filename : str
        Name of file to examine.

    Returns
    -------
    filetype : str
        Type of file.


    """
    # TODO
    # detect the following formats, those supported by RSL
    # 'RADTEC', the SPANDAR radar at Wallops Island, VA
    # 'LASSEN', Darwin Australia
    # 'MCGILL', McGill S-band
    # 'TOGA', DYNAMO project's radar
    # 'RAPIC', Berrimah Australia
    # 'RAINBOW'

    # read the first 12 bytes from the file
    f = open(filename, 'rb')
    begin = f.read(12)
    f.close()

    # MDV, read with read_mdv
    # MDV format signature from MDV FORMAT Interface Control Document (ICD)
    # recond_len1, struct_id, revision_number
    # 1016, 14142, 1
    # import struct
    # mdv_signature = struct.pack('>3i', 1016, 14142, 1)
    mdv_signature = '\x00\x00\x03\xf8\x00\x007>\x00\x00\x00\x01'
    if begin[:12] == mdv_signature:
        return "MDV"

    # NetCDF3, read with read_netcdf
    if begin[:3] == "CDF":
        return "NETCDF3"

    # NetCDF4, read with read_netcdf, contained in a HDF5 container
    # HDF5 format signature from HDF5 specification documentation
    hdf5_signature = '\x89\x48\x44\x46\x0d\x0a\x1a\x0a'
    if begin[:8] == hdf5_signature:
        return "NETCDF4"

    # Other files should be read with read_rsl
    # WSR-88D begin with ARCHIVE2. or AR2V000
    if begin[:9] == 'ARCHIVE2.' or begin[:7] == 'AR2V000':
        return "WSR88D"

    # Universal format has UF in bytes 0,1 or 2,3 or 4,5
    if begin[:2] == "UF" or begin[2:4] == "UF" or begin[4:6] == "UF":
        return "UF"

    # DORADE files
    if begin[:4] == "SSWB" or begin[:4] == "VOLD" or begin[:4] == "COMM":
        return "DORADE"

    # RSL file
    if begin[:3] == "RSL":
        return "RSL"

    # HDF4 file
    # HDF4 format signature from HDF4 specification documentation
    hdf4_signature = '\x0e\x03\x13\x01'
    if begin[:4] == hdf4_signature:
        return "HDF4"

    # SIGMET files
    # SIGMET format is a structure_header with a Product configuration
    # indicator (see section 4.2.47)
    # sigmet_signature = chr(27)
    sigmet_signature = '\x1b'
    if begin[0] == sigmet_signature:
        return "SIGMET"

    # Cannot determine filetype
    return "UNKNOWN"
