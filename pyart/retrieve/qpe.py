#! /usr/bin/env python

import os
import sys
from time import time

import numpy as np
import netCDF4
import matplotlib as mpl
import pylab as pl

import pyart.util.met as met
#import ballsy_masked as ballsy
from pyart import sounding
from pyart.io.nc_utils import save_netcdf_cube, dms_to_d, corner_to_point, \
    rsl_to_arm_netcdf


def get_lowest_elevation(mmcg_ncf, req_moments, corrections):
    pass
    #for each of the moments in req_moments fetch the lowest valid return
