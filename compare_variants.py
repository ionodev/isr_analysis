#!/usr/bin/env python3
"""
Comparison measures for the LPI variants (weeks 7-9 of the project plan).

The five variants (equal weights; variance weights; with the background
column; with the satellite columns; with both) are run on the same test
intervals of zenith-l, each into its own output directory. This script
reads their LPI files and fits and computes, per variant and interval:

    retention   fraction of range-lag cells with an estimate, and of fitted
                gates (150-700 km) with a finite n_e
    sigma       median reported standard deviation of the complex lag-1 ACF at
                200-500 km (precision, meaningful only if the error bars are
                honest)
    honesty     the scatter between neighbouring periods of the same zenith
                cycle divided by the reported error: the robust standard
                deviation of (x_i - x_j)/sqrt(v_i + v_j) for the real part of
                lags 1-5 at 150-600 km. About 1 if the error bars are right;
                above 1 if they understate the error. The ionosphere changes
                a little in 10 s, so it is an upper bound
    residual    satellite power left at the catalogued echoes: for each echo
                above min_snr_db, the lag-1 ACF at the nearest gate minus
                the median of the gates 4-8 away on both sides, divided by
                its reported standard deviation; the median of |z| and the
                fraction above 3. Near 0 if the echo is removed

Bias against a known truth needs synthetic echoes in a quiet interval (the
injection test, to come); until then the quiet reference interval is the
baseline every variant is compared with.

Usage:
    python3 compare_variants.py <variants_dir> <intervals.json> <echoes.npy> <out_dir>

variants_dir holds one directory per variant (lpi_60/zenith-l/ inside);
intervals.json maps an interval name to its list of period indices;
echoes.npy holds the catalogued zenith-l echoes as rows (period index,
range km, SNR dB).
"""

import glob
import json
import os
import sys

import h5py
import numpy as n

# period index -> the unix second in the LPI file name: 10 s periods from the
# first pulse of the eclipse2024 recording
T0_RECORDING = 1712484000
PERIOD_S = 10


def lpi_file(vdir, ai):
    return os.path.join(vdir, "lpi_60", "zenith-l", "lpi-%d.h5" % (T0_RECORDING + PERIOD_S * ai))


def read_lpi(vdir, periods):
    out = {}
    for ai in periods:
        f = lpi_file(vdir, ai)
        if not os.path.exists(f):
            continue
        with h5py.File(f, "r") as h:
            out[ai] = dict(acf=h["acfs_e"][()], var=h["acfs_var"][()], r=h["rgs_km"][()])
    return out


def robust_std(x):
    x = x[n.isfinite(x)]
    if len(x) < 10:
        return n.nan
    return 1.4826 * n.median(n.abs(x - n.median(x)))


def cycles(periods):
    """Runs of consecutive periods: the zenith cycles."""
    p = sorted(periods)
    runs, cur = [], [p[0]]
    for a in p[1:]:
        if a == cur[-1] + 1:
            cur.append(a)
        else:
            runs.append(cur)
            cur = [a]
    runs.append(cur)
    return runs


def measures(L, periods, echoes, min_snr_db=20.0):
    if not L:
        return {}
    r = next(iter(L.values()))["r"]
    acf = n.array([L[a]["acf"] for a in sorted(L)])
    var = n.array([L[a]["var"] for a in sorted(L)])
    res = {}
    res["periods"] = len(L)
    res["retention_cells"] = float(n.isfinite(acf).mean())
    band = (r >= 200) & (r <= 500)
    res["sigma_lag1_200_500km"] = float(n.nanmedian(n.sqrt(var[:, band, 1])))

    # honesty: neighbouring periods of the same cycle
    hb = (r >= 150) & (r <= 600)
    z = []
    for run in cycles(list(L)):
        for a, b in zip(run[:-1], run[1:]):
            # acfs_var is the variance of the complex ACF; its real part
            # carries half of it
            d = (L[a]["acf"][hb, 1:6] - L[b]["acf"][hb, 1:6]).real
            s = n.sqrt((L[a]["var"][hb, 1:6] + L[b]["var"][hb, 1:6]) / 2.0)
            z.append((d / s).ravel())
    res["honesty"] = float(robust_std(n.concatenate(z))) if z else n.nan

    # residual satellite power at the catalogued echoes
    zz = []
    ai_e, rk, snr = echoes[:, 0].astype(int), echoes[:, 1], echoes[:, 2]
    for a in L:
        m = (ai_e == a) & (snr >= min_snr_db)
        x = L[a]["acf"][:, 1].real
        s = n.sqrt(L[a]["var"][:, 1] / 2.0)       # the real part's standard deviation
        for rr in n.unique(n.round(rk[m])):
            g = int(n.argmin(n.abs(r - rr)))
            nb = [g + k for k in list(range(-8, -3)) + list(range(4, 9)) if 0 <= g + k < len(r)]
            if not n.isfinite(x[g]) or not n.isfinite(s[g]) or len(nb) < 4:
                continue
            zz.append((x[g] - n.nanmedian(x[nb])) / s[g])
    zz = n.abs(n.array(zz))
    res["echoes_tested"] = int(len(zz))
    res["residual_median_absz"] = float(n.median(zz)) if len(zz) else n.nan
    res["residual_frac_above_3"] = float((zz > 3).mean()) if len(zz) else n.nan
    return res


def fit_measures(vdir, periods):
    """Fit success at 150-700 km for the fit windows inside the interval."""
    t_lo = T0_RECORDING + PERIOD_S * min(periods)
    t_hi = T0_RECORDING + PERIOD_S * (max(periods) + 1)
    ok, tot, rel = 0, 0, []
    for f in glob.glob(os.path.join(vdir, "lpi_60", "zenith-l", "pp-*.h5")):
        t = int(os.path.basename(f)[3:-3])
        if not (t_lo <= t < t_hi):
            continue
        with h5py.File(f, "r") as h:
            r, ne, dne = h["rgs"][()], h["ne"][()], h["dne"][()]
        b = (r >= 150) & (r <= 700)
        ok += int(n.isfinite(ne[b]).sum())
        tot += int(b.sum())
        rel.append(n.abs(dne[b] / ne[b]))
    out = {"fit_windows": len(rel), "retention_fits": ok / tot if tot else n.nan}
    out["median_dne_over_ne"] = float(n.nanmedian(n.concatenate(rel))) if rel else n.nan
    return out


def main():
    vroot, ivf, ecf, outd = sys.argv[1:5]
    intervals = json.load(open(ivf))
    echoes = n.load(ecf)
    os.makedirs(outd, exist_ok=True)
    rows = []
    for vdir in sorted(glob.glob(os.path.join(vroot, "v*"))):
        v = os.path.basename(vdir)
        for name, per in intervals.items():
            L = read_lpi(vdir, per)
            m = measures(L, per, echoes)
            m.update(fit_measures(vdir, per))
            m.update(variant=v, interval=name)
            rows.append(m)
            print(v, name, {k: (round(x, 3) if isinstance(x, float) else x) for k, x in m.items()
                            if k not in ("variant", "interval")}, flush=True)
    json.dump(rows, open(os.path.join(outd, "variant_measures.json"), "w"), indent=1)
    keys = ["periods", "retention_cells", "retention_fits", "sigma_lag1_200_500km", "honesty",
            "echoes_tested", "residual_median_absz", "residual_frac_above_3", "median_dne_over_ne"]
    with open(os.path.join(outd, "variant_measures.md"), "w") as f:
        f.write("| interval | variant | " + " | ".join(keys) + " |\n")
        f.write("|" + "---|" * (len(keys) + 2) + "\n")
        for m in sorted(rows, key=lambda m: (m["interval"], m["variant"])):
            cells = [("%.3g" % m[k]) if isinstance(m.get(k), float) else str(m.get(k)) for k in keys]
            f.write("| %s | %s | %s |\n" % (m["interval"], m["variant"], " | ".join(cells)))


if __name__ == "__main__":
    main()
