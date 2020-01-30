"""
Reading files using Radx to first convert the file to Cf.Radial format

"""

import os
import tempfile
import subprocess

from ..io.cfradial import read_cfradial
from ..io.common import _test_arguments


def read_radx(filename, radx_dir=None, **kwargs):
    """
    Read a file by first converting it to Cf/Radial using RadxConvert.

    Parameters
    ----------
    filename : str
        Name of file to read using RadxConvert.

    radx_dir : str, optional
        path to the radx install

    Returns
    -------
    radar : Radar
        Radar object.

    """
    # test for non empty kwargs
    _test_arguments(kwargs)
    if radx_dir is not None:
        executable = os.path.join(radx_dir, 'RadxConvert')
    else:
        executable = 'RadxConvert'

    tmpfile = tempfile.mkstemp(suffix='.nc', dir='.')[1]
    head, tail = os.path.split(tmpfile)
    try:
        subprocess.check_call(
            [executable, '-const_ngates',
             '-outdir', head, '-outname', tail, '-f', filename])
        if not os.path.isfile(tmpfile):
            raise IOError(
                'RadxConvert failed to create a file, upgrading to the '
                ' latest version of Radx may be necessary.')
        radar = read_cfradial(tmpfile)
    finally:
        os.remove(tmpfile)
    return radar
