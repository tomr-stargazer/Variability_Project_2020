import pytest
import numpy as np

from ..variability_selection import delta


def test_delta_1():
    # base case: no variability

    m = np.ones(5)
    sigma_m = np.ones(5) * 0.2

    desired_delta = np.zeros(5)

    np.testing.assert_allclose(delta(m, sigma_m), desired_delta)


def test_delta_2():
    # simple case, a little contrived; desired values computed on command line

    m = np.ones(5)
    m[2] = 2
    sigma_m = np.ones(5) * 0.2

    desired_delta = np.array(
        [-1.11803399, -1.11803399, 4.47213595, -1.11803399, -1.11803399]
    )

    np.testing.assert_allclose(delta(m, sigma_m), desired_delta)


def test_delta_3():
    # make sure that nans work as desired

    m = np.ones(5)
    m[2] = np.nan
    m[3] = 2
    sigma_m = np.ones(5) * 0.2

    desired_delta = np.sqrt(4 / 3) * (m - np.nanmean(m)) / sigma_m

    np.testing.assert_allclose(delta(m, sigma_m), desired_delta)
