"""
Tests of the switches added for product v3 (weeks 10-12): the correlated-gate
variance of fit_lpi.range_average (range_avg_rho) and the Monte Carlo check
of it; the defaults leave the range average as it was.

Run with: python -m pytest test_v3_switches.py
"""
import numpy as np

import fit_lpi

RGS = np.array([100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0])


def _data(seed=1):
    rng = np.random.default_rng(seed)
    acf = rng.normal(size=(7, 3)) + 1j * rng.normal(size=(7, 3))
    var = rng.uniform(0.5, 2.0, size=(7, 3))
    return acf, var


def test_rho_none_and_zero_unchanged():
    acf, var = _data()
    a0, v0 = fit_lpi.range_average(acf, var, RGS, 2)
    a1, v1 = fit_lpi.range_average(acf, var, RGS, 2, rho=None)
    a2, v2 = fit_lpi.range_average(acf, var, RGS, 2, rho=(0.0, 0.0))
    np.testing.assert_array_equal(a0, a1)
    np.testing.assert_array_equal(v0, v1)
    np.testing.assert_array_equal(a0, a2)
    np.testing.assert_allclose(v0, v2)


def test_equal_gates_follow_the_formula():
    # M equal gates with correlation rho_k: variance of the mean is
    # v/M * (1 + 2/M sum_k (M-k) rho_k)  (memo 31)
    M, rho = 7, (-0.30, -0.08, -0.06)
    a, v = fit_lpi.range_average(np.ones((M, 1)), np.full((M, 1), 2.0), np.ones(M), 3, rho=rho)
    c = 1 + 2.0 / M * sum((M - k) * r for k, r in enumerate(rho, start=1))
    assert np.isclose(v[3, 0], 2.0 / M * c)
    np.testing.assert_allclose(a[:, 0], 1.0)


def test_monte_carlo():
    # unequal weights and variances: the stated variance matches the scatter of
    # the average of correlated gates
    rng = np.random.default_rng(5)
    rho = (-0.30, -0.08, -0.06)
    n_g = 7
    sd = np.sqrt(np.linspace(0.5, 2.0, n_g))
    C = np.eye(n_g)
    for k, r in enumerate(rho, start=1):
        C += r * (np.eye(n_g, k=k) + np.eye(n_g, k=-k))
    L = np.linalg.cholesky(C)
    x = (L @ rng.normal(size=(n_g, 40000))) * sd[:, None]
    w = RGS**2
    avg = (w[:, None] * x).sum(axis=0) / w.sum()
    _, v = fit_lpi.range_average(np.zeros((n_g, 1)), (sd**2)[:, None], RGS, 3, rho=rho)
    assert abs(v[3, 0] / avg.var() - 1) < 0.03


def test_missing_gate_drops_out():
    acf, var = _data()
    acf[2, 0] = np.nan
    a, v = fit_lpi.range_average(acf, var, RGS, 1, rho=(-0.3,))
    assert np.isnan(v[2, 0]) and np.isnan(a[2, 0])
    # gate 3's window holds gates 3 and 4 at lag 0: one pair 1 apart
    w = RGS[3:5]**2
    s = w * np.sqrt(var[3:5, 0])
    expect = (np.sum(w**2 * var[3:5, 0]) + 2 * -0.3 * s[0] * s[1]) / np.sum(w)**2
    assert np.isclose(v[3, 0], expect)


def test_weights_width_leave_one_out():
    # impulse_rows.masked_running_sigma with a longer window (outlier_lpi's
    # weights_width): each row's level is the mean of |m|^2 over the window
    # around it, wrapping around the period, without the row itself
    import impulse_rows as improws
    rng = np.random.default_rng(3)
    p = rng.exponential(size=(200, 4))
    for width in (10, 30, 100):
        s = improws.masked_running_sigma(p, np.ones(p.shape, bool), width, exclude_self=True)
        for r in (0, 57, 199):
            idx = [(r + k) % 200 for k in range(width // 2 - width + 1, width // 2 + 1) if k != 0]
            assert len(idx) == width - 1
            np.testing.assert_allclose(s[r], np.sqrt(p[idx].mean(axis=0)), rtol=1e-9)


def test_density_constant_by_day(tmp_path):
    import h5py

    import fit_quality as fq

    d = tmp_path / "lpi_60" / "zenith-l"
    d.mkdir(parents=True)
    # 7 April 12 UTC, 8 April 12 UTC, 9 April 12 UTC (no constant of its own)
    for t in (1712491200, 1712577600, 1712664000):
        with h5py.File(d / ("pp-%d.h5" % t), "w") as h:
            h["t0"], h["t1"], h["rgs"] = float(t), float(t + 60), np.array([300.0])
            for k in ("ne", "dne", "Te", "dTe", "Ti", "dTi", "vi", "dvi"):
                h[k] = np.array([1.0])
    F = fq.load_fits(str(tmp_path))
    np.testing.assert_array_equal(F["ne"][:, 0], [5.90e8, 6.32e8, fq.DENSITY_CONSTANT])
    np.testing.assert_array_equal(fq.load_fits(str(tmp_path), constant=1.0)["ne"][:, 0], [1.0, 1.0, 1.0])


def test_row_cap():
    import impulse_rows as ir

    rng = np.random.default_rng(3)
    z = lambda s: (rng.standard_normal(s) + 1j * rng.standard_normal(s)) / np.sqrt(2)
    m = z((2000, 50)) * np.conj(z((2000, 50)))
    p = np.abs(m) ** 2
    m[100, 7] *= 30.0                        # one row with interference
    p[100, 7] = np.abs(m[100, 7]) ** 2
    s = ir.masked_running_sigma(p, np.ones(p.shape, bool), 100, exclude_self=True)
    c = ir.row_cap(s, p, 20.0)
    # the interfered row: its own power over the cap
    np.testing.assert_allclose(c[100, 7], np.sqrt(p[100, 7] / 20.0))
    # clean rows: unchanged except the few above twenty times their level
    changed = (c != s)
    changed[100, 7] = False
    assert changed.mean() < 0.003
    assert np.all(c >= s)
