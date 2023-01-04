""" Unit Tests for Py-ART's retrieve/comp_z.py module. """

import copy

import numpy as np
from numpy.testing import assert_array_equal, assert_equal

import pyart


def test_composite_z():
    # initalize test radar object
    radar = pyart.io.read(pyart.testing.NEXRAD_ARCHIVE_MSG31_FILE)
    ref_field = "reflectivity"
    gatefilter = pyart.filters.GateFilter(radar)
    gatefilter.exclude_transition()
    gatefilter.exclude_below("cross_correlation_ratio", 0.95)

    #################################
    # Trivial first test, all 0s
    #################################

    # make fake radar data with 0 everywhere
    z = np.zeros(radar.fields["reflectivity"]["data"].shape)
    z_old = radar.fields["reflectivity"]

    z_new = copy.deepcopy(z_old)
    z_new["data"] = z.astype("float32")
    radar.add_field("reflectivity", z_new, replace_existing=True)

    compz = pyart.retrieve.composite_reflectivity(
        radar, field=ref_field, gatefilter=gatefilter
    )
    assert_equal(compz.fields["composite_reflectivity"]["data"].max(), 0)

    #################################
    # Insert 1 layer of all 40 dBZ
    #################################

    z = np.zeros(radar.fields["reflectivity"]["data"].shape)

    # choose random sweep, it shouldnt matter which one.
    random_sweep = np.random.randint(0, radar.nsweeps)
    # loop over all measured sweeps
    for sweep in radar.sweep_number["data"]:
        if sweep == random_sweep:
            # get start and stop index numbers
            s_idx = radar.sweep_start_ray_index["data"][sweep]
            e_idx = radar.sweep_end_ray_index["data"][sweep] + 1
            z[s_idx:e_idx, :] = 40

    z_new = copy.deepcopy(z_old)
    z_new["data"] = z.astype("float32")
    radar.add_field("reflectivity", z_new, replace_existing=True)
    compz = pyart.retrieve.composite_reflectivity(
        radar, field=ref_field, gatefilter=gatefilter
    )
    assert_equal(compz.fields["composite_reflectivity"]["data"].max(), 40)

    ##########################################
    # Have dBZ increase according to range bin
    ##########################################

    z = np.zeros(radar.fields["reflectivity"]["data"].shape)

    # choose random sweep, it shouldnt matter which one.
    random_sweep = np.random.randint(0, radar.nsweeps)

    # loop over all measured sweeps
    for sweep in radar.sweep_number["data"]:
        if sweep == random_sweep:
            # get start and stop index numbers
            s_idx = radar.sweep_start_ray_index["data"][sweep]
            e_idx = radar.sweep_end_ray_index["data"][sweep] + 1
            z[s_idx:e_idx, :] = np.tile(
                np.arange(0, z.shape[1])[np.newaxis, :], (e_idx - s_idx, 1)
            )

    z_new = copy.deepcopy(z_old)
    z_new["data"] = z.astype("float32")
    radar.add_field("reflectivity", z_new, replace_existing=True)
    compz = pyart.retrieve.composite_reflectivity(
        radar, field=ref_field, gatefilter=gatefilter
    )

    # choose a random az
    random_az = np.random.randint(0, 720)
    assert_array_equal(
        compz.fields["composite_reflectivity"]["data"][random_az, :],
        np.arange(0, z.shape[1]),
    )
