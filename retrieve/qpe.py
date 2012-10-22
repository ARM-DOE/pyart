#!/usr/bin/python
__author__ = "Scott Collis <scollis@anl.gov>"
__version__ = "Evaluation 0.1"
import matplotlib as mpl
mpl.use('Agg')
import sys
import os
import numpy as np
pyart_dir=os.environ.get('PYART_DIR',os.environ['HOME']+'/python')
sys.path.append(pyart_dir)

import pyart.util.met as met
#import ballsy_masked as ballsy
from time import time
import netCDF4
from pyart import sounding
from pyart.io.nc_utils import save_netcdf_cube, dms_to_d, corner_to_point, rsl_to_arm_netcdf
import pylab as pl

#development start May 29 2012
#E0.1

def get_lowest_elevation(mmcg_ncf, req_moments, corrections):
	#for each of the moments in req_moments fetch the lowest valid return
	



