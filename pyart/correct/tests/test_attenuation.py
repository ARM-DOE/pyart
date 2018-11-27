""" Unit Tests for Py-ART's correct/attenuation.py module. """

# python test_phase_attenuation.py
# to recreate the reference_rays.npz file.

import os

import pyart
import numpy as np
from numpy.testing import assert_allclose

PATH = os.path.dirname(__file__)
REFERENCE_RAYS_FILE = os.path.join(PATH, 'attenuation_rays.npz')
REFERENCE_RAYS_FILE_ZPHI = os.path.join(PATH, 'attenuation_rays_zphi.npz')
REFERENCE_RAYS_FILE_PHILINEAR = os.path.join(
    PATH, 'attenuation_rays_philinear.npz')


def test_attenuation():
    spec_at, cor_z = perform_attenuation()
    ref = np.load(REFERENCE_RAYS_FILE)
    assert_allclose(ref['spec_at'], spec_at['data'], rtol=1e-2, atol=1e-3)
    assert_allclose(ref['cor_z'], cor_z['data'].data, rtol=1e-2, atol=1e-3)


def perform_attenuation_zphi():
    """ Perform attenuation correction on a single ray radar. """
    radar = pyart.testing.make_single_ray_radar()
    a = radar.fields['reflectivity']['data']
    radar.fields['reflectivity']['data'] = np.ma.array(a)
    spec_at, pia_dict, cor_z, spec_diff_at, pida_dict, cor_zdr = (
        pyart.correct.calculate_attenuation_zphi(radar, a_coef=0.06, beta=0.8,
                                                 fzl=4000.0, c=0.15917,
                                                 d=1.0804, doc=0.0,
                                                 temp_ref='fixed_fzl'))
    return spec_at, pia_dict, cor_z, spec_diff_at, pida_dict, cor_zdr


def perform_attenuation_philinear():
    """ Perform attenuation correction on a single ray radar. """
    radar = pyart.testing.make_single_ray_radar()
    a = radar.fields['reflectivity']['data']
    radar.fields['reflectivity']['data'] = np.ma.array(a)
    spec_at, pia_dict, cor_z, spec_diff_at, pida_dict, cor_zdr = (
        pyart.correct.calculate_attenuation_philinear(radar, pia_coef=0.06,
                                                      pida_coef=0.8,
                                                      fzl=4000.0, doc=0.0,
                                                      temp_ref='fixed_fzl'))
    return spec_at, pia_dict, cor_z, spec_diff_at, pida_dict, cor_zdr


def perform_attenuation():
    """ Perform attenuation correction on a single ray radar. """
    radar = pyart.testing.make_single_ray_radar()
    a = radar.fields['reflectivity']['data']
    radar.fields['reflectivity']['data'] = np.ma.array(a)
    spec_at, cor_z = pyart.correct.calculate_attenuation(radar, 0.0)
    return spec_at, cor_z


def save_reference_rays(spec_at, cor_z):
    """ Save the phase processed rays to REFERENCE_RAY_FILE. """
    np.savez(REFERENCE_RAYS_FILE,
             spec_at=spec_at['data'], cor_z=cor_z['data'].data)


def save_reference_rays_zphi(spec_at, pia_dict, cor_z, spec_diff_at,
                             pida_dict, cor_zdr):
    """ Save the phase processed rays to REFERENCE_RAY_FILE_ZPHI. """
    np.savez(REFERENCE_RAYS_FILE_ZPHI,
             spec_at=spec_at['data'], cor_z=cor_z['data'].data,
             pia_dict=pia_dict['data'], spec_diff_at=spec_diff_at['data'],
             pida_dict=pida_dict['data'], cor_zdr=cor_zdr['data'])


def save_reference_rays_philinear(spec_at, pia_dict, cor_z, spec_diff_at,
                                  pida_dict, cor_zdr):
    """ Save the phase processed rays to REFERENCE_RAY_FILE_PHILINEAR. """
    np.savez(REFERENCE_RAYS_FILE_PHILINEAR,
             spec_at=spec_at['data'], cor_z=cor_z['data'].data,
             pia_dict=pia_dict['data'], spec_diff_at=spec_diff_at['data'],
             pida_dict=pida_dict['data'], cor_zdr=cor_zdr['data'])


if __name__ == "__main__":
    spec_at, cor_z = perform_attenuation()
    save_reference_rays(spec_at, cor_z)
    spec_at, pia_dict, cor_z, spec_diff_at, pida_dict, cor_zdr = (
        perform_attenuation_zphi())
    save_reference_rays_zphi(spec_at, pia_dict, cor_z, spec_diff_at,
                             pida_dict, cor_zdr)
    spec_at, pia_dict, cor_z, spec_diff_at, pida_dict, cor_zdr = (
        perform_attenuation_philinear())
    save_reference_rays_philinear(spec_at, pia_dict, cor_z, spec_diff_at,
                                  pida_dict, cor_zdr)
    print("Reference rays saved")


def test_specific_diff_attenuation():
    """ Perform attenuation correction on a single ray radar. """
    spec_at, pia_dict, cor_z, spec_diff_at, pida_dict, cor_zdr = (
        perform_attenuation_zphi())
    ref = np.load(REFERENCE_RAYS_FILE_ZPHI)
    assert_allclose(ref['spec_at'], spec_at['data'], rtol=1e-2, atol=1e-3)
    assert_allclose(ref['cor_z'], cor_z['data'], rtol=1e-2, atol=1e-3)
    assert_allclose(ref['pia_dict'], pia_dict['data'], rtol=1e-2, atol=1e-3)
    assert_allclose(ref['spec_diff_at'], spec_diff_at['data'], 
                    rtol=1e-2, atol=1e-3)
    assert_allclose(ref['pida_dict'], pida_dict['data'], rtol=1e-2, atol=1e-3)
    assert_allclose(ref['cor_zdr'], cor_zdr['data'], rtol=1e-2, atol=1e-3)

    spec_at, pia_dict, cor_z, spec_diff_at, pida_dict, cor_zdr = (
        perform_attenuation_philinear())
    ref = np.load(REFERENCE_RAYS_FILE_PHILINEAR)
    assert_allclose(ref['spec_at'], spec_at['data'], rtol=1e-2, atol=1e-3)
    assert_allclose(ref['cor_z'], cor_z['data'], rtol=1e-2, atol=1e-3)
    assert_allclose(ref['pia_dict'], pia_dict['data'], rtol=1e-2, atol=1e-3)
    assert_allclose(ref['spec_diff_at'], spec_diff_at['data'], 
                    rtol=1e-2, atol=1e-3)
    assert_allclose(ref['pida_dict'], pida_dict['data'], rtol=1e-2, atol=1e-3)
    assert_allclose(ref['cor_zdr'], cor_zdr['data'], rtol=1e-2, atol=1e-3)
