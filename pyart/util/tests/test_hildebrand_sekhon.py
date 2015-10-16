""" Unit tests for the hildebrand_sekhon module. """

from pyart.util import hildebrand_sekhon

import numpy as np
from numpy.testing import assert_almost_equal

SAMPLE_SPECTRUM = np.array(
    [11.089, 10.386, 11.081, 9.696, 9.21, 11.103, 9.425,
     10.222, 8.393, 11.67, 9.632, 10.096, 9.807, 10.374,
     8.065, 8.682, 10.164, 320.206, 480.898, 640.818, 490.483,
     318.578, 9.682, 10.035, 9.134, 8.554, 9.274, 9.463,
     9.501, 10.036, 9.211, 8.265])


def test_hildebrand_sekhon():
    noise_params = hildebrand_sekhon.estimate_noise_hs74(SAMPLE_SPECTRUM)
    mean, thresh, var, nnoise = noise_params
    assert_almost_equal(mean, 10, 0)
    assert_almost_equal(thresh, 12, 0)
    assert_almost_equal(var, 0.79, 2)
    assert nnoise == 27
