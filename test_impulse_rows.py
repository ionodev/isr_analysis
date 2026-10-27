import numpy as n

import impulse_rows as ir
from outlier_lpi import fft_lpf, simple_decimator

N = 10000
L = 20               # filter_len of the eclipse configuration
RG = 60              # samples per row
M0, M1 = 13, 137     # the rows lpi_files uses (sample0=800, sample1=8200)
GC, LAST = 1000, 8000


def noise(rng, size, scale=1.0):
    return (scale * (rng.normal(size=size) + 1j * rng.normal(size=size)) / n.sqrt(2)).astype(n.complex64)


def rows(z, lpf, decim, lag):
    """The plain lagged products of one pulse in rows, as lpi_files forms them."""
    z = n.array(z, dtype=n.complex64)
    z[0:GC] = 0.0
    z[LAST:] = 0.0
    zf = lpf.lpf(z)
    zf[0:GC] = n.nan
    zf[LAST:] = n.nan
    return decim.decimate(zf[0:N - lag] * n.conj(zf[lag:N]))[M0:M1]


def test_detect_finds_an_impulse_but_not_an_extended_echo():
    rng = n.random.default_rng(1)
    z = noise(rng, N)
    z[3000] += 8.0                                   # 64 times the noise power
    z[5000:5480] += 10.0 * n.exp(1j * rng.uniform(0, 2 * n.pi, 480))   # a satellite-like echo, 100x
    hot = ir.detect(z, GC, LAST)
    assert hot[2998:3003].all()
    assert hot.sum() == 5
    assert not hot[4900:5600].any()


def test_widen_is_the_span_of_the_filter():
    lpf = fft_lpf(N, f0=1.2e5, L=L)
    for s in (2000, 4321):
        d = n.zeros(N, dtype=n.complex64)
        d[s] = 1.0
        reach = n.abs(lpf.lpf(d)) > 1e-6
        hot = n.zeros(N, dtype=bool)
        hot[s] = True
        w = ir.widen(hot, L)
        assert not (reach & ~w).any()               # every reached sample is flagged
        assert n.flatnonzero(w).tolist() == list(range(s - L, s + L))


def test_kept_rows_carry_no_impulse_and_rows_away_are_untouched():
    rng = n.random.default_rng(2)
    lpf = fft_lpf(N, f0=1.2e5, L=L)
    decim = simple_decimator(L=N, dec=RG)
    z = noise(rng, N)
    zi = z.copy()
    s = 4000
    zi[s] += 30.0
    touched = ir.widen(ir.detect(zi, GC, LAST), L)
    for lag in (10, 30, 100, 250, 450):
        hit = ir.hit_rows(touched, lag, decim, M0, M1)
        a, b = rows(z, lpf, decim, lag), rows(zi, lpf, decim, lag)
        scale = n.nanmedian(n.abs(a))
        ok = ~hit & n.isfinite(a)
        # fully removed: the kept rows are those of the data without the impulse
        assert n.max(n.abs(a[ok] - b[ok])) < 1e-3 * scale
        # and the removed rows are the ones the impulse reaches
        assert (n.abs(a - b)[hit] > 1e-3 * scale).any()
        # rows away from the impulse are untouched: only rows within reach of
        # s (products at t) and of s - lag (products at t + lag)
        rr = (n.flatnonzero(hit) + M0) * RG
        near = (n.abs(rr + RG / 2 - s) <= RG / 2 + L + 3) | (n.abs(rr + RG / 2 - (s - lag)) <= RG / 2 + L + 3)
        assert near.all()
        assert hit.sum() <= 2 * (2 + (2 * L + 5) // RG + 1)


def test_the_acf_is_unbiased_after_removal_and_biased_without():
    # a correlated process (complex AR(1)) in every pulse, and impulses at
    # random samples in most pulses; the lag-10 ACF as the mean of the rows
    rng = n.random.default_rng(3)
    lpf = fft_lpf(N, f0=1.2e5, L=L)
    decim = simple_decimator(L=N, dec=RG)
    lag = 10
    rho = 0.9
    est_clean, est_naive, est_removed = [], [], []
    for p in range(150):
        w = noise(rng, N)
        x = n.empty(N, dtype=n.complex64)
        x[0] = w[0]
        c = n.sqrt(1 - rho**2)
        for i in range(1, N):
            x[i] = rho * x[i - 1] + c * w[i]
        xi = x.copy()
        for s in rng.integers(GC + 100, LAST - 100, size=8):
            xi[s] += 25.0 * n.exp(1j * rng.uniform(0, 2 * n.pi))
        touched = ir.widen(ir.detect(xi, GC, LAST), L)
        hit = ir.hit_rows(touched, lag, decim, M0, M1)
        a, b = rows(x, lpf, decim, lag), rows(xi, lpf, decim, lag)
        ok = n.isfinite(a)
        est_clean.append(n.mean(a[ok]))
        est_naive.append(n.mean(b[ok]))
        est_removed.append(n.mean(b[ok & ~hit]))
    est_clean, est_naive, est_removed = map(n.array, (est_clean, est_naive, est_removed))
    # paired with the same pulses without impulses
    d_rem = (est_removed - est_clean).real
    d_naive = (est_naive - est_clean).real
    # unbiased within the error bar of the mean (3 sigma)
    assert abs(n.mean(d_rem)) < 3 * n.std(d_rem) / n.sqrt(len(d_rem))
    # without removal the impulses' own lagged products bias it
    assert n.mean(d_naive) > 5 * n.std(d_naive) / n.sqrt(len(d_naive))


def test_masked_running_sigma():
    rng = n.random.default_rng(4)
    p = rng.uniform(1, 2, size=(40, 5))
    keep = n.ones(p.shape, dtype=bool)
    # all rows kept: the plain running mean of lpi_files
    plain = n.empty_like(p)
    WF = n.fft.fft(n.repeat(0.1, 10), 40)
    for ri in range(5):
        plain[:, ri] = n.roll(n.sqrt(n.fft.ifft(WF * n.fft.fft(p[:, ri])).real), -5)
    assert n.allclose(ir.masked_running_sigma(p, keep, 10), plain)
    # a removed row with a huge value does not raise its neighbours' level
    q = p.copy()
    q[20, 2] = 1e6
    keep[20, 2] = False
    s = ir.masked_running_sigma(q, keep, 10)
    assert s[:, 2].max() < 2.0
    assert n.allclose(s[:, [0, 1, 3, 4]], plain[:, [0, 1, 3, 4]])
