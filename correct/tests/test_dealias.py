""" Tests for the dealias module in pyart.correct """

import datetime
import os.path

import netCDF4
import numpy as np
from numpy.testing import assert_array_equal

import pyart
from pyart.correct import dealias


DIR = os.path.dirname(__file__)
RSLNAME = os.path.join(DIR, "sample.sigmet")
SONDNAME = os.path.join(DIR, 'sgpinterpolatedsondeC1.c1.20110510.000000.cdf')
NCFNAME = os.path.join(DIR, 'sample.nc')
REFNAME = os.path.join(DIR, 'dealias_reference.npy')

#################
# Dealias Tests #
#################


def test_find_time():
    """ find_time_in_interp_sonde test """
    target = datetime.datetime(2011, 5, 10, 11, 30, 8)
    interp_sounde = netCDF4.Dataset(SONDNAME)
    t = dealias.find_time_in_interp_sonde(interp_sounde, target)
    height, speed, direction = t

    assert height.shape == (316,)
    assert speed.shape == (316,)
    assert direction.shape == (316,)

    assert height.dtype == 'float32'
    assert speed.dtype == 'float32'
    assert direction.dtype == 'float32'

    assert round(height[100], 2) == 2.32
    assert round(speed[100], 2) == 15.54
    assert round(direction[100], 2) == 231.8


def test_dealias_rsl():
    """ Dealias data from a RSL object """

    # read in the data
    radarobj = pyart.io.py4dd.RSL_anyformat_to_radar(RSLNAME)
    radar = pyart.io.radar.Radar(radarobj)

    # find and extract sonde data
    target = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    interp_sounde = netCDF4.Dataset(SONDNAME)
    t = dealias.find_time_in_interp_sonde(interp_sounde, target)
    height, speed, direction = t

    # perform dealiasing
    dealias_data = dealias.dealias(radar, height * 1000.0, speed, direction,
                                   target)

    # compare against known good data
    reference_data = np.load(REFNAME)
    assert_array_equal(reference_data, dealias_data['data'].data)


def test_dealias_ncf():
    """ Dealias data from a NetCDF file """

    # read in the data
    ncf = netCDF4.Dataset(NCFNAME)
    radar = pyart.io.radar.Radar(ncf)

    # find and extract the sonde data
    target = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    interp_sounde = netCDF4.Dataset(SONDNAME)
    t = dealias.find_time_in_interp_sonde(interp_sounde, target)
    height, speed, direction = t

    # perform dealiasing
    dealias_data = dealias.dealias(radar, height * 1000.0, speed, direction,
                                   target)

    # compare against know good data
    reference_data = np.load(REFNAME)
    assert_array_equal(reference_data, dealias_data['data'].data)


def test_dealias_compare_rsl():
    """ Compare dealiasing with and without a rsl_radar """

    # read in the data
    radarobj = pyart.io.py4dd.RSL_anyformat_to_radar(RSLNAME)
    radar = pyart.io.radar.Radar(radarobj)

    # find and extract the sonde data
    target = netCDF4.num2date(radar.time['data'][0], radar.time['units'])
    interp_sounde = netCDF4.Dataset(SONDNAME)
    t = dealias.find_time_in_interp_sonde(interp_sounde, target)
    height, speed, direction = t

    # perform dealiasing, with recreation of a new RSL Radar
    dealias_data = dealias.dealias(radar, height * 1000.0, speed,
                                   direction, target)

    # perform dealiasing, without creating a new RSL Radar

    # the original RSL file does not set the altitude correct...
    radarobj.contents.volumes[0].sweeps[0].rays[0].h.alt = 340
    radarobj.contents.volumes[1].sweeps[0].rays[0].h.alt = 340
    dealias_data2 = dealias.dealias(radar, height * 1000.0, speed,
                                    direction, target, rsl_radar=radarobj)

    assert_array_equal(dealias_data2['data'].data, dealias_data['data'].data)
