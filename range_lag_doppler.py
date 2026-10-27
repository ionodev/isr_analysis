#!/usr/bin/env python3
"""
Range-lag and range-Doppler images of one fit window, in two runs (for
example no satellite handling, and the tracked satellite columns), from the
measured ACF the fit sees (ion_line_fits.window_acf: the 1/variance weighted
mean over the window, receiver gain, range averaging).

    range-lag      the real part of the ACF, each gate divided by its
                   zero-lag power (as in fit_lpi's pp_fit_*.png)
    range-Doppler  the power spectrum of each gate, 10 log10, scaled by r^2
                   so that the topside is visible; the same colour scale in
                   both runs

The LPI has no zero lag (its first lag is 10 us). The zero lag is
extrapolated from the first three lags (quadratic, R0 = 3 R1 - 3 R2 + R3),
which also keeps the noise power, present only at zero lag, out of the
spectrum. The spectrum is the Fourier transform over lags -450..450 us with a
Hann taper (resolution about 2 kHz).

Usage:
    python3 range_lag_doppler.py <before_lpi_dir> <after_lpi_dir> <pp t0> <out.png> [--labels before after]
"""

import datetime
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as n

import ion_line_fits as ilf

R_MAX = 800.0
F_MAX = 15e3


def zero_lag(acf):
    return 3 * acf[:, 0] - 3 * acf[:, 1] + acf[:, 2]


def spectra(acf, lag, f):
    """Power spectrum of every gate on f (Hz), from the ACF at lags 0..lag[-1]."""
    tau = n.concatenate([[0.0], lag])
    r = n.concatenate([zero_lag(acf)[:, None], acf], axis=1)
    r = n.nan_to_num(r)
    w = 0.5 * (1 + n.cos(n.pi * tau / (tau[-1] + tau[1])))
    k = n.exp(-2j * n.pi * n.outer(tau, f))
    s = 2 * n.real((r * w) @ k) - n.real(r[:, :1])
    return s * (tau[1] - tau[0])


def main():
    args = sys.argv[1:]
    labels = ["before", "after"]
    if "--labels" in args:
        i = args.index("--labels")
        labels = args[i + 1:i + 3]
        args = args[:i]
    dirs, t0, out = args[0:2], int(args[2]), args[3]
    f = n.linspace(-F_MAX, F_MAX, 301)
    runs = []
    for d in dirs:
        with h5py.File("%s/pp-%d.h5" % (d, t0), "r") as pp:
            ta, tb, nf = pp["t0"][()], pp["t1"][()], pp["n_lpi_files"][()]
        rgs, lag, acf, var = ilf.window_acf(d, ta, nf)
        runs.append((rgs, lag, acf))
    rgs, lag = runs[0][0], runs[0][1]
    g = (rgs > 90) & (rgs < R_MAX)
    S = [10 * n.log10(n.maximum(spectra(a[g], lag, f) * (rgs[g, None] / 300.0) ** 2, 1e-3)) for _, _, a in runs]
    vmax = n.nanpercentile(S[0], 99.5)
    fig, ax = plt.subplots(2, 2, figsize=(12, 10), sharey=True)
    for col, (_, _, a) in enumerate(runs):
        with n.errstate(invalid="ignore", divide="ignore"):
            norm = (a / zero_lag(a)[:, None]).real
        m = ax[0, col].pcolormesh(lag * 1e6, rgs[g], norm[g], vmin=-0.3, vmax=1.0, cmap="viridis", shading="nearest")
        fig.colorbar(m, ax=ax[0, col], label="ACF / zero-lag power (real part)")
        ax[0, col].set_title("range-lag, %s" % labels[col])
        ax[0, col].set_xlabel("lag ($\\mu$s)")
        m = ax[1, col].pcolormesh(f / 1e3, rgs[g], S[col], vmin=vmax - 25, vmax=vmax, cmap="inferno", shading="nearest")
        fig.colorbar(m, ax=ax[1, col], label="power spectral density $\\times (r/300\\,\\mathrm{km})^2$ (dB)")
        ax[1, col].set_title("range-Doppler, %s" % labels[col])
        ax[1, col].set_xlabel("Doppler frequency (kHz)")
    for a in ax[:, 0]:
        a.set_ylabel("range (km)")
    for a in ax.flat:
        for r in (535, 590):
            a.axhline(r, color="w", lw=0.6, ls="--")
    ts = datetime.datetime.fromtimestamp(ta, datetime.timezone.utc)
    fig.suptitle("Zenith, %s UTC (%.0f s); dashed: the Starlink shells, 535-590 km" % (ts.strftime("%Y-%m-%d %H:%M:%S"), tb - ta))
    fig.tight_layout()
    fig.savefig(out, dpi=110)


if __name__ == "__main__":
    main()
