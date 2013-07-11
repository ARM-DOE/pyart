""" Unit Tests for Py-ART's correct/attenuation.py module. """

# python test_phase_attenuation.py
# to recreate the reference_rays.npz file.

import os

import pyart
import numpy as np
from numpy.testing import assert_allclose

PATH = os.path.dirname(__file__)
REFERENCE_RAYS_FILE = os.path.join(PATH, 'attenuation_rays.npz')


def test_attenuation():
    spec_at, cor_z = perform_attenuation()
    ref = np.load(REFERENCE_RAYS_FILE)
    assert_allclose(ref['spec_at'], spec_at['data'])
    assert_allclose(ref['cor_z'], cor_z['data'].data)


def perform_attenuation():
    """ Perform attenuation correction on a single ray radar. """
    radar = pyart.testing.make_single_ray_radar()
    # hack to get test to run, TODO fix the attenuation code to not need this
    radar.fields['proc_dp_phase_shift'] = radar.fields['dp_phase_shift']
    a = radar.fields['reflectivity_horizontal']['data']
    radar.fields['reflectivity_horizontal']['data'] = np.ma.array(a)
    spec_at, cor_z = pyart.correct.calculate_attenuation(radar, 0.0)
    return spec_at, cor_z


def save_reference_rays(spec_at, cor_z):
    """ Save the phase processed rays to REFERENCE_RAY_FILE. """
    np.savez(REFERENCE_RAYS_FILE,
             spec_at=spec_at['data'], cor_z=cor_z['data'].data)


if __name__ == "__main__":
    spec_at, cor_z = perform_attenuation()
    save_reference_rays(spec_at, cor_z)
    print("Reference rays saved")
