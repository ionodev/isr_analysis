#!/usr/bin/env python3
"""
Monte Carlo of the bias of variance weights estimated from the data (Memo 51).

Each row is m = s + n, a constant complex signal s plus circular complex
Gaussian noise of unit variance (per-row SNR = |s|^2). The weight of a row is
1 / (running mean of |m|^2 over 10 rows), as outlier_lpi's variance weights,
either including the row itself (as in the products so far) or leaving it out
(weights_exclude_self); a share of the rows is removed at random (as impulse
rows are). The weighted mean of the rows estimates s.

First-order theory: with the row in its own weight, Cov(|m|^2, m) = s sigma^2,
so the estimate is pulled towards zero by 1 / ((1 + SNR) N), with N the kept
rows in the running mean; left out, the weight is independent of the row's
noise and there is no bias.

A lagged product is not "signal plus Gaussian noise" but the product of two
Gaussian samples, m = z1 z2* (plasma with correlation rho between them, plus
receiver noise). Its |m|^2 = |z1|^2 |z2|^2 is heavy tailed and more strongly
correlated with m, and the bias is larger: about -22 to -29 % for N = 10,
nearly independent of the SNR (the second table). In outlier_lpi every
weight below the median is raised to the median, so where the plasma is
weak most rows share one weight and the bias is small; that is why the
measured bias falls with altitude (Memo 51).

    python3 weight_bias_mc.py
"""
import numpy as n

import impulse_rows as ir


def main(pulses=2000, columns=200, width=10, seed=1):
    rng = n.random.default_rng(seed)
    print("per-row SNR  removed  bias (row included)  bias (row left out)  theory")
    for snr in (0.01, 0.1, 1.0):
        for frac in (0.0, 0.29):
            s = n.sqrt(snr)
            m = s + (rng.standard_normal((pulses, columns)) + 1j * rng.standard_normal((pulses, columns))) / n.sqrt(2)
            keep = rng.uniform(size=m.shape) > frac
            p = n.abs(m) ** 2
            out = []
            for ex in (False, True):
                w = keep / ir.masked_running_sigma(p, keep, width, exclude_self=ex) ** 2
                out.append(n.mean(((w * m).sum(0) / w.sum(0)).real) / s - 1)
            print("%10.2f  %7.2f  %+19.3f  %+19.3f  %+6.3f" % (
                snr, frac, out[0], out[1], -1 / ((1 + snr) * width * (1 - frac))))

    def cg(shape):
        return (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / n.sqrt(2)

    print("\nlagged products m = z1 z2*: plasma SNR, correlation rho  bias (row included)  bias (row left out)")
    for snr in (0.3, 3.0, 30.0):
        for rho in (0.8, 0.4):
            a = cg((pulses, columns))
            b = rho * a + n.sqrt(1 - rho ** 2) * cg((pulses, columns))
            m = (n.sqrt(snr) * a + cg(a.shape)) * n.conj(n.sqrt(snr) * b + cg(a.shape))
            keep = n.ones(m.shape, dtype=bool)
            out = []
            for ex in (False, True):
                w = 1 / ir.masked_running_sigma(n.abs(m) ** 2, keep, width, exclude_self=ex) ** 2
                out.append(n.mean(((w * m).sum(0) / w.sum(0)).real) / (snr * rho) - 1)
            print("%10.1f  %4.1f  %+19.3f  %+19.3f" % (snr, rho, out[0], out[1]))


if __name__ == "__main__":
    main()
