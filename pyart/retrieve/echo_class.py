"""
pyart.retrieve.echo_class
=========================

"""

import numpy as np

from ..config import get_fillvalue, get_field_name
from ..retrieve import echo_steiner


def steiner_conv_strat(grid, dx=500.0, dy=500.0, intense=42.0,
                       work_level=3000.0, peak_relation='default',
                       area_relation='medium', bkg_rad=11000.0,
                       use_intense=True, fill_value=None,
                       refl_field=None):
    """
    
    Parameters
    ----------
    grid : Grid
    
    Optional parameters
    -------------------
    dx, dy : float
        The x- and y-dimension resolutions, respectively.
    intense : float
    
    work_level : float
    
    peak_relation : 'default' or 'sgp'
    
    area_relation : 'small', 'medium', 'large', or 'sgp'
    
    bkg_rad : float
        Background radius.
    use_intense : bool
        True to use the intensity criteria.
    fill_value : float
    
    refl_field : str
    
    Returns
    -------
    eclass : dict
        Steiner convective-stratiform classification dictionary.
    """
    
    # Get fill value
    if fill_value is None:
        fill_value = get_fillvalue()
        
    # Parse field parameters
    if refl_field is None:
        refl_field = get_field_name('corrected_reflectivity')
        
    # Get axes
    x = grid.axes['x_disp']['data']
    y = grid.axes['y_disp']['data']
    z = grid.axes['z_disp']['data']
    
    # Get reflectivity data
    ze = np.copy(grid.fields[refl_field]['data'])
    ze = np.ma.filled(ze, fill_value).astype(np.float64)
    
    # Call Fortran routine
    eclass = echo_steiner.classify(ze, x, y, z, dx=dx, dy=dy, bkg_rad=bkg_rad,
                                   work_level=work_level, intense=intense,
                                   peak_relation=peak_relation,
                                   area_relation=area_relation,
                                   use_intense=use_intense,
                                   fill_value=fill_value)
    
    
    return {'data': eclass.astype(np.int32),
            'standard_name': 'echo_classification',
            'long_name': 'Steiner echo classification',
            'valid_min': 0,
            'valid_max': 2,
            'comment': ('Convective-stratiform echo '
                        'classification based on '
                        'Steiner et al. (1995)')}
    