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


def get_lowest_elevation(mmcg_ncf, req_moments, corrections):
    pass
    #for each of the moments in req_moments fetch the lowest valid return
