"""
Input/output routines common to many file formats.

"""

import bz2
import gzip

import fsspec
import netCDF4
import numpy as np


def prepare_for_read(filename, storage_options={"anon": True}):
    """
    Return a file like object read for reading.

    Open a file for reading in binary mode with transparent decompression of
    Gzip and BZip2 files. The resulting file-like object should be closed.

    Parameters
    ----------
    filename : str or file-like object
        Filename or file-like object which will be opened. File-like objects
        will not be examined for compressed data.

    storage_options : dict, optional
        Parameters passed to the backend file-system such as Google Cloud Storage,
        Amazon Web Service S3.

    Returns
    -------
    file_like : file-like object
        File like object from which data can be read.

    """
    # if a file-like object was provided, return
    if hasattr(filename, "read"):  # file-like object
        return filename

    # look for compressed data by examining the first few bytes
    fh = fsspec.open(filename, mode="rb", compression="infer", **storage_options).open()
    magic = fh.read(3)
    fh.close()

    # If the data is still compressed, use gunzip/bz2 to uncompress the data
    if magic.startswith(b"\x1f\x8b"):
        return gzip.GzipFile(filename, "rb")

    if magic.startswith(b"BZh"):
        return bz2.BZ2File(filename, "rb")

    return fsspec.open(
        filename, mode="rb", compression="infer", **storage_options
    ).open()


def stringarray_to_chararray(arr, numchars=None):
    """
    Convert an string array to a character array with one extra dimension.

    Parameters
    ----------
    arr : array
        Array with numpy dtype 'SN', where N is the number of characters
        in the string.

    numchars : int
        Number of characters used to represent the string. If numchar > N
        the results will be padded on the right with blanks. The default,
        None will use N.

    Returns
    -------
    chararr : array
        Array with dtype 'S1' and shape = arr.shape + (numchars, ).

    """
    carr = netCDF4.stringtochar(arr)
    if numchars is None:
        return carr

    arr_numchars = carr.shape[-1]
    if numchars <= arr_numchars:
        raise ValueError("numchars must be >= %i" % (arr_numchars))
    chararr = np.zeros(arr.shape + (numchars,), dtype="S1")
    chararr[..., :arr_numchars] = carr[:]
    return chararr


def _test_arguments(dic):
    """Issue a warning if receive non-empty argument dict."""
    if dic:
        import warnings

        warnings.warn("Unexpected arguments: %s" % dic.keys())


def make_time_unit_str(dtobj):
    """Return a time unit string from a datetime object."""
    return "seconds since " + dtobj.strftime("%Y-%m-%dT%H:%M:%SZ")
