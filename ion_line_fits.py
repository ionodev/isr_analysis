#!/usr/bin/env python3
"""
The ion-line fit at a few altitudes for one fit window, in two runs (for
example no satellite handling, and the tracked satellite columns): the
measured ACF the fit saw, with its error bars, the fitted model, and the
ion-line spectrum of that model.

The measured ACF is rebuilt from the run's LPI files the way
fit_lpi.fit_lpifiles builds it (1/variance weighted mean over the window,
ground-clutter-cancelled ACF below 300 km, receiver gain, range averaging);
the model is the stored fit (pp-*.h5: Te, Ti, vi and the zero-lag power
ne/ne_const). As a check the gate is refitted with fit_lpi.fit_acf on the
rebuilt ACF and both parameter sets are printed.

Usage:
    python3 ion_line_fits.py <before_lpi_dir> <after_lpi_dir> <pp t0> <out.png> <range_km> [<range_km> ...] [--labels before after]

The two directories are the lpi_<postfix>/<channel> folders of the runs.
"""

import datetime
import glob
import os
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as n
import scipy.constants as c

import fit_lpi as fl

MAX_DT = 60
RANGE_LIMITS = n.array([0, 300, 700, 1500])
RANGE_AVG = n.array([1, 3, 5])


def window_acf(d, t0, n_files):
    """The averaged, gain-corrected, range-averaged ACF and its variance."""
    fs = []
    for f in sorted(glob.glob("%s/lpi-*.h5" % d)):
        with h5py.File(f, "r") as h:
            if t0 - 1e-3 <= h["i0"][()] < t0 + MAX_DT:
                fs.append(f)
    if len(fs) != n_files:
        raise SystemExit("%s: %d LPI files in the window, the fit used %d" % (d, len(fs), n_files))
    sa, sw, gains = 0.0, 0.0, []
    for f in fs:
        with h5py.File(f, "r") as h:
            rgs, lag = h["rgs_km"][()], h["lags"][()]
            a = h["acfs_e"][()]
            cut = n.where(rgs > 300)[0][0]
            a[:cut, :] = h["acfs_g"][()][:cut, :] / 2.0
            v = h["acfs_var"][()]
            gains.append(h["alpha"][()])
        sa = sa + n.nan_to_num(a / v)
        sw = sw + n.nan_to_num(1 / v)
    with n.errstate(invalid="ignore", divide="ignore"):
        acf, var = sa / sw, 1 / sw
    g = n.nanmedian(gains)
    acf, var = acf / g, var / g ** 2
    a0, v0 = n.copy(acf), n.copy(var)
    idx = [n.argmin(n.abs(rgs - x)) for x in RANGE_LIMITS]
    for i, ra in enumerate(RANGE_AVG):
        aa, vv = fl.range_average(a0, v0, rgs, ra)
        acf[idx[i]:idx[i + 1]] = aa[idx[i]:idx[i + 1]]
        var[idx[i]:idx[i + 1]] = vv[idx[i]:idx[i + 1]]
    return rgs, lag, acf, var


def spectrum(te, ti, mol, vi, p0, f):
    """Ion-line power spectrum of the model ACF (both signs of lag), on f in Hz."""
    tau = fl.ilf.lag
    r = p0 * fl.model_acf(te, ti, mol, vi, tau)
    k = n.exp(-2j * n.pi * n.outer(f, tau))
    s = 2 * n.real(k @ r) - n.real(r[0])
    return s * (tau[1] - tau[0])


def main():
    args = sys.argv[1:]
    labels = ["before", "after"]
    if "--labels" in args:
        i = args.index("--labels")
        labels = args[i + 1:i + 3]
        args = args[:i]
    dirs, t0, out = args[0:2], int(args[2]), args[3]
    ranges = [float(x) for x in args[4:]]
    fl._init_tables(440.2e6, table_dir=os.path.join(os.path.dirname(os.path.abspath(fl.__file__)), "data"))

    runs = []
    for d in dirs:
        pp = h5py.File("%s/pp-%d.h5" % (d, t0), "r")
        P = {k: pp[k][()] for k in ["rgs", "Te", "Ti", "vi", "ne", "ne_const", "dTe", "dTi", "dvi", "dne", "n_lpi_files", "t0", "t1"]}
        rgs, lag, acf, var = window_acf(d, P["t0"], P["n_lpi_files"])
        runs.append((P, rgs, lag, acf, var))

    f = n.linspace(-12e3, 12e3, 481)
    mlag = n.linspace(0, lag[-1], 200)
    nr = len(ranges)
    fig, ax = plt.subplots(nr, 3, figsize=(14, 3.6 * nr), squeeze=False)
    for row, rk in enumerate(ranges):
        ri = n.argmin(n.abs(runs[0][1] - rk))
        ylim = 0
        for col, (P, rgs, lag, acf, var) in enumerate(runs):
            h = rgs[ri]
            mol = fl.mh_molecular_ion_fraction(n.array([h]))[0]
            te, ti, vi = P["Te"][ri], P["Ti"][ri], P["vi"][ri]
            p0 = P["ne"][ri] / P["ne_const"][ri]
            m = p0 * fl.model_acf(te, ti, mol, vi, mlag)
            sd = n.sqrt(var[ri].real / 2)            # 1 sigma of the real or imaginary part
            a = ax[row, col]
            a.errorbar(lag * 1e6, acf[ri].real, yerr=sd, fmt="o", ms=3, color="C0", alpha=0.7, label="measured, real")
            a.errorbar(lag * 1e6, acf[ri].imag, yerr=sd, fmt="o", ms=3, color="C1", alpha=0.7, label="measured, imaginary")
            a.plot(mlag * 1e6, m.real, "-", color="C0", lw=2, label="fit, real")
            a.plot(mlag * 1e6, m.imag, "-", color="C1", lw=2, label="fit, imaginary")
            a.axhline(0, color="k", lw=0.5)
            a.set_title("%s, %.0f km\n$T_e$=%.0f$\\pm$%.0f K, $T_i$=%.0f$\\pm$%.0f K, $v_i$=%.0f$\\pm$%.0f m/s"
                        % (labels[col], h, te, P["dTe"][ri], ti, P["dTi"][ri], vi, P["dvi"][ri], ), fontsize=9)
            a.set_xlabel("lag ($\\mu$s)")
            a.set_ylabel("ACF (receiver units)")
            ylim = max(ylim, n.nanmax(n.abs(acf[ri].real[:5])), abs(p0))
            ax[row, 2].plot(f / 1e3, spectrum(te, ti, mol, vi, p0, f), lw=2, color=["0.4", "C3"][col],
                            label="%s: zero-lag power %.0f" % (labels[col], p0))
            # the check: refit the rebuilt ACF
            x, _, sg, _ = fl.fit_acf(n.copy(acf[ri]), lag, h, n.copy(var[ri]))
            print("%s %.0f km: stored Te %.0f Ti %.0f vi %.0f P0 %.4g | refit Te %.0f Ti %.0f vi %.0f P0 %.4g"
                  % (labels[col], h, te, ti, vi, p0, x[0] * x[1], x[1], x[2], x[3]))
        for col in range(2):
            ax[row, col].set_ylim(-0.6 * ylim, 1.3 * ylim)
            ax[row, col].set_xlim(0, lag[-1] * 1e6 + 10)
        ax[row, 2].set_title("fitted ion-line spectrum, %.0f km" % h, fontsize=9)
        ax[row, 2].set_xlabel("Doppler frequency (kHz)")
        ax[row, 2].set_ylabel("power spectral density (receiver units)")
        ax[row, 2].legend(fontsize=8)
    ax[0, 0].legend(fontsize=7, loc="upper right")
    ts = datetime.datetime.fromtimestamp(runs[0][0]["t0"], datetime.timezone.utc)
    fig.suptitle("Ion-line fits, zenith, %s UTC (%.0f s)" % (ts.strftime("%Y-%m-%d %H:%M:%S"), runs[0][0]["t1"] - runs[0][0]["t0"]))
    fig.tight_layout()
    fig.savefig(out, dpi=110)


if __name__ == "__main__":
    main()
