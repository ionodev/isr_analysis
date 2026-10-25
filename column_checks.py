#!/usr/bin/env python3
"""
Per-column checks of the satellite columns (Memo 46, option C): is there an
echo in the column, and is it a coherent point target?

A satellite is a coherent point target: within a pulse its echo is a delayed,
Doppler-shifted copy of the transmitted pulse, so the column's fitted
amplitude is s_j = A exp(-i 2 pi df tau_j) at every lag tau_j (memo 14): the
same size at every lag, only its phase turns with the small Doppler error df.
A column on plasma (a false detection where the plasma is) shrinks with lag
as the plasma's ACF does (to near zero beyond 200 us in the F region); a
column on noise fits zero at every lag. The test uses the 45 amplitudes and
their variances (sat_amp_e, sat_amp_var):

    fit     df on a grid of 25 Hz over +-3 kHz; for each, the inverse-variance
            weighted mean of s_j exp(+i 2 pi df tau_j) over the lags, A(df);
            the df of largest |A|. sigma_A^2 = 1 / sum(1/v_j)
    T1      an echo is there: S = |A| / sigma_A above s1. For a column on
            noise with independent lags |A|^2/sigma_A^2 is exponential, so
            P(S > x) = exp(-x^2) per trial of df; with about 3 independent df
            trials (6 kHz times 450 us) and a significance of 1e-3 that is
            x = 2.83. The lags are not independent (the same samples enter
            every lag), so the threshold is taken from the pure-plasma
            control's random columns instead (their 99.9th percentile of S),
            when given (--null)
    T2      the echo is coherent: with the fitted df, the de-rotated weighted
            mean over the long lags 31-45 (310-450 us), projected on that
            over the short lags 1-15 (10-150 us), rho = Re(A_long A_short*) /
            |A_short|^2. A point target gives 1, plasma about 0 (its ACF
            has decayed). The column fails T2 when rho is below 0.5 (the
            midpoint) by more than z = 3.09 (one-sided 1e-3) of its error
            sigma_rho = sqrt(var_long/2 + rho^2 var_short/2) / |A_short|

Classes (classify): "coherent" (T1 and T2 pass), "decaying" (T1 passes, T2
fails: not a point target's shape, e.g. plasma or interference), "no echo"
(T1 fails and no other column of the pulse overlaps it: a false detection,
or a tracked pulse beyond the echo) and "shared" (T1 fails, but another
column of the pulse overlaps it, so the echo may be split between them).

Limit of T1: the variance weights include the echo's own power (a running
mean of |m|^2 over ten pulses), so a bright echo's amplitude is uncertain in
proportion to itself and S stays of order ten even for a bright echo.

The residual (sat_resid, when the LPI files have it): the mean residual
power over the noise in the column's own rows after the solve, lags 1-10.

Detected and tracked columns are told apart with the catalogue: a column
whose pulse and delay are a catalogue detection's is "detected", the rest
were added by the tracker (or are random control columns).

Usage:
    python3 column_checks.py <lpi_dir> <out.npz> [--catalogue DIR] [--null NULL.npz] [--jobs 16]
"""

import argparse
import glob
import os
from multiprocessing import Pool

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import numpy as n

KM_PER_US = 299792458.0 / 2 / 1e9
DF = n.arange(-3000.0, 3000.0 + 1e-9, 25.0)
SHORT = slice(0, 15)
LONG = slice(30, 45)
S1_THEORY = float(n.sqrt(n.log(3 / 1e-3)))
Z2 = 3.09
RHO0 = 0.5
CHANNEL_DELAY_US = 10.62


def coherence(s, v, tau):
    """Per column: S (T1 statistic), rho, sigma_rho, df, |A|. s, v: (n_col, n_lag)."""
    w = n.where(n.isfinite(v) & (v > 0) & n.isfinite(s), 1.0 / n.where(v > 0, v, 1.0), 0.0)
    s0 = n.where(w > 0, s, 0.0)
    ph = n.exp(2j * n.pi * DF[:, None] * tau[None, :])           # (n_df, n_lag)
    num = (s0 * w) @ ph.T                                         # (n_col, n_df)
    W = w.sum(axis=1)
    with n.errstate(invalid="ignore", divide="ignore"):
        A = num / W[:, None]
        k = n.nanargmax(n.abs(n.nan_to_num(A)), axis=1)
        Abest = A[n.arange(len(k)), k]
        sigA = n.sqrt(1.0 / W)
        S = n.abs(Abest) / sigA
        rot = n.exp(2j * n.pi * DF[k][:, None] * tau[None, :])
        d = s0 * rot
        ws, wl = w[:, SHORT], w[:, LONG]
        As = (d[:, SHORT] * ws).sum(1) / ws.sum(1)
        Al = (d[:, LONG] * wl).sum(1) / wl.sum(1)
        vs, vl = 1.0 / ws.sum(1), 1.0 / wl.sum(1)
        rho = n.real(Al * n.conj(As)) / n.abs(As) ** 2
        srho = n.sqrt(vl / 2 + rho ** 2 * vs / 2) / n.abs(As)
    return S, rho, srho, DF[k], n.abs(Abest)


def one_file(args):
    f, cat_dir = args
    with h5py.File(f, "r") as h:
        if "sat_amp_e" not in h or h["sat_amp_e"].shape[0] == 0:
            return None
        s = h["sat_amp_e"][()].astype(n.complex128)
        v = h["sat_amp_var"][()].astype(n.float64)
        tau = h["lags"][()]
        keys = h["sat_keys"][()]
        delay = h["sat_delay_samples"][()]
        dop = h["sat_doppler_hz"][()]
        resid = n.nanmean(h["sat_resid"][()][:, 0:10], axis=1) if "sat_resid" in h else n.full(len(keys), n.nan)
        t0 = float(h["t0"][()])
    S, rho, srho, df, amp = coherence(s, v, tau)
    # the nearest other column of the same pulse, in delay
    nearest = n.full(len(keys), n.inf)
    order = n.lexsort((delay, keys))
    ks, ds = keys[order], delay[order]
    same = ks[1:] == ks[:-1]
    gap = n.where(same, n.abs(n.diff(ds)), n.inf)
    near_sorted = n.minimum(n.concatenate([[n.inf], gap]), n.concatenate([gap, [n.inf]]))
    nearest[order] = near_sorted
    detected = n.zeros(len(keys), dtype=bool)
    snr = n.full(len(keys), n.nan)
    if cat_dir:
        import satellite_columns as satcol
        cat = satcol.CatalogueDetections(cat_dir, "zenith-l", CHANNEL_DELAY_US)
        rows = cat.rows(int(t0 * 1e6), int((t0 + 10) * 1e6) + 40000)
        lut = {}
        for r in rows:
            lut.setdefault(int(r[0]), []).append((r[1], r[3]))
        for i, (k, d) in enumerate(zip(keys, delay)):
            for dd, sn in lut.get(int(k), []):
                if abs(dd - d) < 1e-3:
                    detected[i] = True
                    snr[i] = sn
    return dict(t0=n.full(len(keys), t0), key=keys, delay=delay, doppler=dop, S=S, rho=rho, srho=srho,
                df=df, amp=amp, detected=detected, snr=snr, resid=resid, nearest=nearest)


def classify(S, rho, srho, nearest, s1, overlap=480.0):
    """
    coherent   T1 and T2 pass
    decaying   T1 passes, T2 fails: an echo, but not a point target's shape
    no echo    T1 fails and no other column of the pulse lies within the
               pulse length (overlap samples): nothing in the column
    shared     T1 fails, but another column of the pulse overlaps it: the
               two can share one echo (their amplitudes are then anti-
               correlated and each uncertain), so the column alone cannot
               tell whether it holds one
    """
    cls = n.full(len(S), "coherent", dtype=object)
    low = ~(S > s1)
    decays = (rho + Z2 * srho < RHO0) & ~low
    cls[low & (nearest >= overlap)] = "no echo"
    cls[low & (nearest < overlap)] = "shared"
    cls[decays] = "decaying"
    return cls


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("lpi_dir")
    p.add_argument("out")
    p.add_argument("--catalogue")
    p.add_argument("--jobs", type=int, default=16)
    p.add_argument("--glob", default="lpi-*.h5")
    args = p.parse_args()
    files = sorted(glob.glob(os.path.join(args.lpi_dir, args.glob)))
    with Pool(args.jobs) as pool:
        res = [r for r in pool.map(one_file, [(f, args.catalogue) for f in files], chunksize=4) if r is not None]
    out = {k: n.concatenate([r[k] for r in res]) for k in res[0]}
    out["range_km"] = out["delay"] * KM_PER_US
    n.savez(args.out, **out)
    print("%d files, %d columns, %d detected" % (len(res), len(out["S"]), int(out["detected"].sum())))


if __name__ == "__main__":
    main()
