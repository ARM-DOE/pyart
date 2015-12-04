"""
pyart.retrieve.echo_class
=========================

steiner_conv_strat

.. autosummary::
    :toctree: generated/

    steiner_conv_strat

"""

import numpy as np

from ..config import get_fillvalue, get_field_name
from ..exceptions import MissingOptionalDependency
try:
    from . import _echo_steiner
    _F90_EXTENSIONS_AVAILABLE = True
except ImportError:
    _F90_EXTENSIONS_AVAILABLE = False


def steiner_conv_strat(grid, dx=None, dy=None, intense=42.0,
                       work_level=3000.0, peak_relation='default',
                       area_relation='medium', bkg_rad=11000.0,
                       use_intense=True, fill_value=None,
                       refl_field=None):
    """
    Partition reflectivity into convective-stratiform using the Steiner et
    al. (1995) algorithm.

    Parameters
    ----------
    grid : Grid
        Grid containing reflectivity field to partition.

    Other Parameters
    ----------------
    dx, dy : float
        The x- and y-dimension resolutions in meters, respectively.  If None
        the resolution is determined from the first two axes values.
    intense : float
        The intensity value in dBZ. Grid points with a reflectivity
        value greater or equal to the intensity are automatically
        flagged as convective. See reference for more information.
    work_level : float
        The working level (separation altitude) in meters. This is the height
        at which the partitioning will be done, and should minimize bright band
        contamination. See reference for more information.
    peak_relation : 'default' or 'sgp'
        The peakedness relation. See reference for more information.
    area_relation : 'small', 'medium', 'large', or 'sgp'
        The convective area relation. See reference for more information.
    bkg_rad : float
        The background radius in meters. See reference for more information.
    use_intense : bool
        True to use the intensity criteria.
    fill_value : float
         Missing value used to signify bad data points. A value of None
         will use the default fill value as defined in the Py-ART
         configuration file.
    refl_field : str
         Field in grid to use as the reflectivity during partitioning. None
         will use the default reflectivity field name from the Py-ART
         configuration file.

    Returns
    -------
    eclass : dict
        Steiner convective-stratiform classification dictionary.

    References
    ----------
    Steiner, M. R., R. A. Houze Jr., and S. E. Yuter, 1995: Climatological
    Characterization of Three-Dimensional Storm Structure from Operational
    Radar and Rain Gauge Data. J. Appl. Meteor., 34, 1978-2007.
    """
    # check that Fortran extensions is available
    if not _F90_EXTENSIONS_AVAILABLE:
        raise MissingOptionalDependency(
            "Py-ART must be build on a system with a Fortran compiler to use" +
            "the steiner_conv_strat function.")

    # Get fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field parameters
    if refl_field is None:
        refl_field = get_field_name('reflectivity')

    # parse dx and dy
    if dx is None:
        dx = grid.regular_x['data'][1] - grid.regular_x['data'][0]
    if dy is None:
        dy = grid.regular_y['data'][1] - grid.regular_y['data'][0]

    # Get coordinates
    x = grid.regular_x['data']
    y = grid.regular_y['data']
    z = grid.regular_z['data']

    # Get reflectivity data
    ze = np.ma.copy(grid.fields[refl_field]['data'])
    ze = np.ma.filled(ze, fill_value).astype(np.float64)

    # Call Fortran routine
    eclass = _echo_steiner.classify(
        ze, x, y, z, dx=dx, dy=dy, bkg_rad=bkg_rad, work_level=work_level,
        intense=intense, peak_relation=peak_relation,
        area_relation=area_relation, use_intense=use_intense,
        fill_value=fill_value)

    return {'data': eclass.astype(np.int32),
            'standard_name': 'echo_classification',
            'long_name': 'Steiner echo classification',
            'valid_min': 0,
            'valid_max': 2,
            'comment_1': ('Convective-stratiform echo '
                          'classification based on '
                          'Steiner et al. (1995)'),
            'comment_2': ('0 = Undefined, 1 = Stratiform, '
                          '2 = Convective')}
