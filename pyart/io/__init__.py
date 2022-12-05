"""
Functions to read and write radar and grid data to and from a number of file
formats.

In most cases the :py:func:`pyart.io.read` function should be used to read
in radar data from a file. In certain cases the function the read function
for the format in question should be used.

"""

from .arm_sonde import read_arm_sonde, read_arm_sonde_vap
from .auto_read import read
from .cfradial import read_cfradial, write_cfradial
from .chl import read_chl
from .common import prepare_for_read
from .grid_io import read_grid, write_grid
from .mdv_grid import read_grid_mdv, write_grid_mdv
from .mdv_radar import read_mdv
from .nexrad_archive import read_nexrad_archive
from .nexrad_cdm import read_nexrad_cdm
from .nexradl3_read import read_nexrad_level3
from .output_to_geotiff import write_grid_geotiff
from .rsl import read_rsl
from .sigmet import read_sigmet
from .uf import read_uf
from .uf_write import write_uf

__all__ = [s for s in dir() if not s.startswith("_")]
