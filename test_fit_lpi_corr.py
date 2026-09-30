"""
The correlated-gate variance factor of fit_lpi's range average (memo 31).
"""
import numpy as np

import fit_lpi


def test_independent_gates_factor_one():
    w = np.linspace(1, 2, 7)[:, None]
    v = np.linspace(1, 3, 7)[:, None]
    assert fit_lpi._corr_factor(w, v, np.zeros((7, 1, 3)))[0] == 1.0


def test_equal_gates_matches_formula():
    # a mean of M equally noisy gates: c = 1 + (2/M) sum_k (M-k) rho_k
    M, r = 10, [-0.30, -0.05, -0.02]
    rho = np.zeros((M, 1, 3))
    rho[:, :, :] = r
    c = fit_lpi._corr_factor(np.ones((M, 1)), np.ones((M, 1)), rho)[0]
    assert np.isclose(c, 1 + 2 / M * sum((M - k) * rk for k, rk in zip([1, 2, 3], r)))


def test_matches_monte_carlo_with_weights():
    rng = np.random.default_rng(1)
    M = 8
    w = (1 + 0.1 * np.arange(M))[:, None]
    s = (1 + 0.2 * np.arange(M))
    C = np.eye(M)
    for i in range(M - 1):
        C[i, i + 1] = C[i + 1, i] = -0.3
    C = C * np.outer(s, s)
    x = rng.multivariate_normal(np.zeros(M), C, size=200000)
    est = np.var(x @ w[:, 0]) / np.sum(w[:, 0] ** 2 * s ** 2)
    rho = np.zeros((M, 1, 3))
    rho[:, :, 0] = -0.3
    assert np.isclose(fit_lpi._corr_factor(w, (s ** 2)[:, None], rho)[0], est, rtol=0.02)


def test_nan_gate_drops_out():
    rho = np.zeros((5, 1, 3))
    rho[:, :, 0] = -0.3
    v = np.ones((5, 1))
    v[2] = np.nan
    f = fit_lpi._corr_factor(np.ones((5, 1)), v, rho)[0]
    # gates 0-1 and 3-4 are the only neighbour pairs left: (4 - 2*0.3*... ) / 4
    assert np.isclose(f, (4 + 2 * 2 * (-0.3)) / 4)


def test_complex_typed_weights_as_in_fit_lpifiles():
    # range_weight in fit_lpifiles is a complex array holding r^2
    M = 6
    w = (np.arange(M) + 100.0)[:, None] ** 2
    rho = np.zeros((M, 1, 3))
    rho[:, :, 0] = -0.3
    v = np.ones((M, 1))
    f_real = fit_lpi._corr_factor(w, v, rho)
    f_cplx = fit_lpi._corr_factor(w.astype(np.complex64), v.astype(np.float32), rho)
    assert f_cplx.dtype.kind == "f"
    assert np.allclose(f_real, f_cplx, rtol=1e-6)
    avg_var = np.ones((1,), dtype=np.float64)
    avg_var *= f_cplx   # as fit_lpifiles applies it: must not raise
