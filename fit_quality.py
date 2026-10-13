#!/usr/bin/env python3
"""
Fitted plasma parameters with quality flags, for the comparison of the LPI
variants and the eclipse plots.

load_fits() reads the fit_lpi output (pp-*.h5) of a directory into arrays
(time, range) and scales n_e by the density constant (Memo 43). flags()
marks each fitted gate:

    FAILED    a parameter or its uncertainty is not finite
    BOUND     the fit ends at a limit of its search: Te/Ti at 0.99 or 3,
              Ti at 150 or 3000 K, |vi| at 1500 m/s (fit_lpi's bounds)
    NOISY     dne/ne or dTe/Te above 0.5
    SATELLITE a catalogued echo of at least 20 dB within one gate of it
              during the fit window (for information: with the satellite
              columns it is modelled, without them it may bias the gate)

A gate is good when it carries none of the first three.
"""

import glob
import os

import h5py
import numpy as n

DENSITY_CONSTANT = 7.15e8          # Memo 43: lpi_60, 100 kHz, filter length 20
T0_RECORDING = 1712484000          # the first period of the eclipse2024 recording
PERIOD_S = 10

FAILED, BOUND, NOISY, SATELLITE = 1, 2, 4, 8


def load_fits(d, t_range=None, constant=DENSITY_CONSTANT):
    files = sorted(glob.glob(os.path.join(d, "lpi_60", "zenith-l", "pp-*.h5")))
    rows = []
    for f in files:
        t = int(os.path.basename(f)[3:-3])
        if t_range is not None and not (t_range[0] <= t < t_range[1]):
            continue
        with h5py.File(f, "r") as h:
            rows.append(dict(t0=float(h["t0"][()]), t1=float(h["t1"][()]), r=h["rgs"][()],
                             ne=h["ne"][()] * constant, dne=h["dne"][()] * constant,
                             Te=h["Te"][()], dTe=h["dTe"][()], Ti=h["Ti"][()], dTi=h["dTi"][()],
                             vi=h["vi"][()], dvi=h["dvi"][()]))
    if not rows:
        return None
    out = {k: n.array([x[k] for x in rows]) for k in rows[0] if k != "r"}
    out["r"] = rows[0]["r"]
    return out


def flags(F, echoes=None, min_snr_db=20.0):
    """Flag bits per (window, gate). echoes: rows (period index, range km, SNR dB)."""
    f = n.zeros(F["ne"].shape, dtype=n.int16)
    vals = [F[k] for k in ("ne", "dne", "Te", "dTe", "Ti", "dTi", "vi", "dvi")]
    bad = n.zeros(f.shape, bool)
    for v in vals:
        bad |= ~n.isfinite(v)
    f[bad] |= FAILED
    with n.errstate(invalid="ignore", divide="ignore"):
        ratio = F["Te"] / F["Ti"]
        bound = (ratio <= 0.995) | (ratio >= 2.98) | (F["Ti"] <= 151) | (F["Ti"] >= 2990) \
            | (n.abs(F["vi"]) >= 1490)
        noisy = (n.abs(F["dne"] / F["ne"]) > 0.5) | (n.abs(F["dTe"] / F["Te"]) > 0.5)
    f[bound & ~bad] |= BOUND
    f[noisy & ~bad] |= NOISY
    if echoes is not None:
        ai, rk, snr = echoes[:, 0], echoes[:, 1], echoes[:, 2]
        gate = n.median(n.diff(F["r"]))
        strong = snr >= min_snr_db
        for w in range(len(F["t0"])):
            a0 = (F["t0"][w] - T0_RECORDING) / PERIOD_S
            a1 = (F["t1"][w] - T0_RECORDING) / PERIOD_S
            m = strong & (ai >= a0) & (ai < a1)
            for r in n.unique(rk[m]):
                f[w, n.abs(F["r"] - r) <= gate] |= SATELLITE
    return f


def good(f):
    return (f & (FAILED | BOUND | NOISY)) == 0
