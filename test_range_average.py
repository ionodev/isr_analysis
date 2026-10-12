"""
Tests of fit_lpi.range_average, the r^2 weighted range average of the ACFs.

Only gates with data at a lag count in that lag's weights, the variance is
sum(w^2 v)/sum(w)^2, an empty or masked centre gate gives NaN, and the
window is cut at the ends of the range axis.

Run with: python -m pytest test_range_average.py
"""
import numpy as np

import fit_lpi

RGS = np.array([100.0, 110.0, 120.0, 130.0, 140.0])


def _data():
    acf = np.array([[1.0 + 1j, 2.0], [2.0, 3.0], [3.0, 4.0], [4.0, 5.0], [5.0, 6.0]], dtype=complex)
    var = np.array([[1.0, 2.0], [2.0, 1.0], [3.0, 2.0], [4.0, 1.0], [5.0, 3.0]])
    return acf, var


def _expect(acf, var, ri, ra, lag):
    # direct sum over the finite gates of the window
    g = [i for i in range(max(0, ri - ra), min(len(RGS), ri + ra + 1)) if np.isfinite(acf[i, lag])]
    w = RGS[g]**2
    return np.sum(w * acf[g, lag]) / np.sum(w), np.sum(w**2 * var[g, lag]) / np.sum(w)**2


def test_ra0_is_identity():
    acf, var = _data()
    a, v = fit_lpi.range_average(acf, var, RGS, 0)
    np.testing.assert_array_equal(a, acf)
    np.testing.assert_array_equal(v, var)


def test_weights_and_variance():
    acf, var = _data()
    a, v = fit_lpi.range_average(acf, var, RGS, 1)
    for ri in range(len(RGS)):
        for lag in range(2):
            ea, ev = _expect(acf, var, ri, 1, lag)
            assert np.isclose(a[ri, lag], ea)
            assert np.isclose(v[ri, lag], ev)
    # equal weights and variances: the variance of a mean of 3 is v/3
    a, v = fit_lpi.range_average(np.ones((5, 1)), np.full((5, 1), 3.0), np.ones(5), 1)
    np.testing.assert_allclose(v[1:4, 0], 1.0)
    np.testing.assert_allclose(a[:, 0], 1.0)


def test_window_edges():
    acf, var = _data()
    a, v = fit_lpi.range_average(acf, var, RGS, 2)
    # first gate: gates 0-2 only; last gate: gates 2-4 only
    w = RGS[0:3]**2
    assert np.isclose(a[0, 1], np.sum(w * acf[0:3, 1]) / np.sum(w))
    w = RGS[2:5]**2
    assert np.isclose(a[4, 1], np.sum(w * acf[2:5, 1]) / np.sum(w))
    assert np.isclose(v[4, 1], np.sum(w**2 * var[2:5, 1]) / np.sum(w)**2)
    # a window wider than the axis takes every gate
    a, v = fit_lpi.range_average(acf, var, RGS, 10)
    w = RGS**2
    np.testing.assert_allclose(a[:, 1], np.sum(w * acf[:, 1]) / np.sum(w))


def test_nan_neighbour_does_not_bias():
    # a constant ACF must stay constant next to a gate without data
    acf = np.ones((5, 2), dtype=complex)
    var = np.ones((5, 2))
    acf[1, :] = np.nan
    var[1, :] = np.nan
    a, v = fit_lpi.range_average(acf, var, RGS, 1)
    np.testing.assert_allclose(a[[0, 2, 3, 4], :], 1.0)
    # the gate below gate 2 does not count in its variance either
    w = RGS[2:4]**2
    assert np.isclose(v[2, 0], np.sum(w**2) / np.sum(w)**2)


def test_empty_or_masked_centre_is_nan():
    acf, var = _data()
    # an empty lowest gate, and a masked gate in the middle
    acf[0, :] = np.nan
    var[0, :] = np.inf
    acf[2, :] = np.nan
    var[2, :] = np.nan
    # a centre gate with data at one lag only
    acf[3, 1] = np.nan
    a, v = fit_lpi.range_average(acf, var, RGS, 1)
    assert np.all(np.isnan(a[[0, 2], :])) and np.all(np.isnan(v[[0, 2], :]))
    assert np.isnan(a[3, 1]) and np.isnan(v[3, 1])
    # lag 0 of gate 3 averages gates 3 and 4 only
    ea, ev = _expect(acf, var, 3, 1, 0)
    assert np.isclose(a[3, 0], ea) and np.isclose(v[3, 0], ev)
    # a gate whose neighbours are all empty keeps its own value
    assert np.isclose(a[1, 0], acf[1, 0]) and np.isclose(v[1, 0], var[1, 0])
