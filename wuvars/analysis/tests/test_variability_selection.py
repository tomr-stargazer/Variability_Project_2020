import pytest
import numpy as np

from ..variability_selection import delta
from ..variability_selection import S_threeband
from ..variability_selection import S_twoband
from ..variability_selection import reduced_chisq


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


def test_S_threeband_1():
    # no variability
    j = np.ones(5)
    sigma_j = np.ones(5) * 0.2
    h = np.ones(5)
    sigma_h = np.ones(5) * 0.2
    k = np.ones(5)
    sigma_k = np.ones(5) * 0.2

    desired_S_threeband = 0

    assert S_threeband(j, sigma_j, h, sigma_h, k, sigma_k) == desired_S_threeband


def test_S_threeband_2():
    # positively correlated "variability"
    j = np.ones(5)
    sigma_j = np.ones(5) * 0.2
    h = np.ones(5)
    sigma_h = np.ones(5) * 0.2
    k = np.ones(5)
    sigma_k = np.ones(5) * 0.2

    j[2] = 2
    h[2] = 2
    k[2] = 2

    desired_S_threeband = 5.366563145999494

    np.testing.assert_equal(
        S_threeband(j, sigma_j, h, sigma_h, k, sigma_k), desired_S_threeband
    )


def test_S_threeband_3():
    # negatively correlated "variability"
    j = np.ones(5)
    sigma_j = np.ones(5) * 0.2
    h = np.ones(5)
    sigma_h = np.ones(5) * 0.2
    k = np.ones(5)
    sigma_k = np.ones(5) * 0.2

    j[2] = 2
    h[2] = -2

    desired_S_threeband = -3.098386676965933

    np.testing.assert_equal(
        S_threeband(j, sigma_j, h, sigma_h, k, sigma_k), desired_S_threeband
    )


def test_S_threeband_4():
    # are we dealing properly with nans?
    j = np.ones(5)
    sigma_j = np.ones(5) * 0.2
    h = np.ones(5)
    sigma_h = np.ones(5) * 0.2
    k = np.ones(5)
    sigma_k = np.ones(5) * 0.2

    j[2] = 2
    h[2] = 2
    k[3] = np.nan

    desired_S_threeband = 2.3094010767585025

    np.testing.assert_equal(
        S_threeband(j, sigma_j, h, sigma_h, k, sigma_k), desired_S_threeband
    )


def test_S_threeband_5():
    # are we dealing properly with nans when there's only one valid measurement?
    j = np.ones(5)
    sigma_j = np.ones(5) * 0.2
    h = np.ones(5)
    sigma_h = np.ones(5) * 0.2
    k = np.ones(5) * np.nan
    sigma_k = np.ones(5) * 0.2

    k[2] = 2

    desired_S_threeband = np.nan

    np.testing.assert_equal(
        S_threeband(j, sigma_j, h, sigma_h, k, sigma_k), desired_S_threeband
    )


def test_S_twoband_1():
    # zero variability
    j = np.ones(5)
    sigma_j = np.ones(5) * 0.2
    h = np.ones(5)
    sigma_h = np.ones(5) * 0.2

    desired_S_twoband = 0

    np.testing.assert_equal(S_twoband(j, sigma_j, h, sigma_h), desired_S_twoband)


def test_S_twoband_2():
    # one 'errant' datapoint
    j = np.ones(5)
    sigma_j = np.ones(5) * 0.2
    h = np.ones(5)
    sigma_h = np.ones(5) * 0.2
    j[2] = 2

    desired_S_twoband = 0

    np.testing.assert_equal(S_twoband(j, sigma_j, h, sigma_h), desired_S_twoband)


def test_S_twoband_3():
    # positively correlated variability
    j = np.ones(5)
    sigma_j = np.ones(5) * 0.2
    h = np.ones(5)
    sigma_h = np.ones(5) * 0.2
    j[2] = 2
    h[2] = 2

    desired_S_twoband = 1.788854381999832

    np.testing.assert_equal(S_twoband(j, sigma_j, h, sigma_h), desired_S_twoband)


def test_reduced_chisq_1():
    mag = np.ones(5)
    err = np.ones(5) * 0.2

    desired_reduced_chisq = 0

    np.testing.assert_equal(reduced_chisq(mag, err), desired_reduced_chisq)


def test_reduced_chisq_2():
    mag = np.ones(5)
    err = np.ones(5) * 0.2

    mag[2] = 2

    desired_reduced_chisq = 5

    np.testing.assert_equal(reduced_chisq(mag, err), desired_reduced_chisq)


def test_reduced_chisq_3():
    mag = np.arange(5) * 1.0
    err = np.ones(5) * 0.2

    mag[0] = np.nan
    mag[2] = np.nan

    desired_reduced_chisq = 58.33333333333332

    np.testing.assert_equal(reduced_chisq(mag, err), desired_reduced_chisq)


def test_reduced_chisq_3():
    mag = np.ones(5)
    err = np.ones(5) * 0.2

    mag[0] = np.nan
    mag[2] = np.nan
    mag[3] = np.nan
    mag[4] = np.nan

    desired_reduced_chisq = np.nan

    np.testing.assert_equal(reduced_chisq(mag, err), desired_reduced_chisq)
