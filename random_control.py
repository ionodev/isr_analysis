#!/usr/bin/env python3
"""
The pure-plasma control (Memo 46, option F): what satellite columns where no
echo is do to the ACF and the fits.

Runs of the quiet reference interval with tracked satellite columns (the
base) and with, in addition, random columns (RandomControlColumns) at
several rates per period. The random columns of a run are those in its LPI
files (sat_keys, sat_delay_samples, sat_doppler_hz) that the base does not
have. Per rate:

    columns   random columns per period
    shift     at each gate (period, gate) holding a random column: the ACF
              with them minus the base, real part, in the base's standard
              deviation (z), mean over lags 1-5 and over all lags; and in
              the standard deviation of the difference itself, sqrt(var with
              - var base), as both estimates use the same data and the base
              is the more precise one (z_diff): mean about 0 if a column on
              plasma costs no bias, robust spread about 1 if the reported
              growth of the error bar is the real one
    sigma     median ratio of the reported standard deviations (lag 1) at
              the gates with columns and at all gates 150-800 km; the price
              of a column: the growth of the variance at its gate per column
    honesty   compare_variants' measure of the run with the columns
    fits      60 s fits at 150-800 km: good gates, median shift of n_e, T_e,
              T_i, v_i in the base fit's standard deviation, error ratio

Usage:
    python3 random_control.py <base_dir> <out_dir> <fig_dir> <rate>=<run_dir> ...
"""

import json
import os
import sys

import h5py
import numpy as n

import compare_variants as cv
import fit_quality as fq

KM_PER_US = 299792458.0 / 2 / 1e9
PERIODS = list(range(3906, 3924))
BANDS = [(150, 300), (300, 450), (450, 600), (600, 800)]


def read(d, periods=PERIODS):
    out = {}
    for ai in periods:
        f = cv.lpi_file(d, ai)
        if not os.path.exists(f):
            continue
        with h5py.File(f, "r") as h:
            cols = set()
            if "sat_keys" in h:
                cols = set(zip(h["sat_keys"][()].tolist(), n.round(h["sat_delay_samples"][()], 6).tolist(),
                               n.round(h["sat_doppler_hz"][()], 6).tolist()))
            out[ai] = dict(acf=h["acfs_e"][()], var=h["acfs_var"][()], r=h["rgs_km"][()], cols=cols)
    return out


def measure(B, R):
    common = sorted(set(B) & set(R))
    r = B[common[0]]["r"]
    rows = []
    zs, zd, sr_col, n_cols, dvar = [], [], [], [], []
    hit_all = []
    for a in common:
        extra = R[a]["cols"] - B[a]["cols"]
        n_cols.append(len(extra))
        gates = {}
        for k, d, f in extra:
            g = int(n.argmin(n.abs(r - d * KM_PER_US)))
            gates[g] = gates.get(g, 0) + 1
        hit = n.zeros(len(r), dtype=int)
        for g, c in gates.items():
            hit[g] = c
        hit_all.append(hit)
        gi = n.array(sorted(gates), dtype=int)
        if len(gi) == 0:
            continue
        dA = (R[a]["acf"][gi] - B[a]["acf"][gi]).real
        vb, vr = B[a]["var"][gi] / 2.0, R[a]["var"][gi] / 2.0
        zs.append(dA / n.sqrt(vb))
        with n.errstate(invalid="ignore"):
            zd.append(dA / n.sqrt(vr - vb))
        sr_col.append(n.sqrt(vr[:, 1] / vb[:, 1]))
        dvar.append(((vr[:, 1] - vb[:, 1]) / vb[:, 1]) / n.array([gates[g] for g in gi]))
    zs, zd = n.concatenate(zs), n.concatenate(zd)
    hit_all = n.array(hit_all)
    band = (r >= 150) & (r <= 800)
    ratio_all = n.sqrt(n.array([R[a]["var"][band, 1] / B[a]["var"][band, 1] for a in common]))
    res = dict(periods=len(common), columns_per_period=float(n.mean(n_cols)),
               gates_with_columns=int(sum((h > 0).sum() for h in hit_all)),
               shift_lags1_5_mean=float(n.nanmean(zs[:, 0:5])),
               shift_lags1_5_err=float(n.nanstd(n.nanmean(zs[:, 0:5], axis=1)) / n.sqrt(len(zs))),
               shift_all_mean=float(n.nanmean(zs)),
               shift_per_lag=n.nanmean(zs, axis=0).tolist(),
               shift_err_per_lag=(n.nanstd(zs, axis=0) / n.sqrt(n.sum(n.isfinite(zs), axis=0))).tolist(),
               zdiff_mean=float(n.nanmean(zd[:, 0:5])),
               zdiff_robust_std=float(cv.robust_std(zd[:, 0:5].ravel())),
               zdiff_std=float(n.nanstd(zd[:, 0:5])),
               sigma_ratio_at_columns=float(n.nanmedian(n.concatenate(sr_col))),
               sigma_ratio_all_150_800=float(n.nanmedian(ratio_all)),
               var_growth_per_column=float(n.nanmedian(n.concatenate(dvar))),
               var_growth_per_column_mean=float(n.nanmean(n.concatenate(dvar))))
    for name, L in (("base", B), ("run", R)):
        m = cv.measures({a: dict(acf=L[a]["acf"], var=L[a]["var"], r=L[a]["r"]) for a in common}, common,
                        n.zeros((0, 3)))
        res["honesty_" + name] = m["honesty"]
        res["sigma_lag1_200_500_" + name] = m["sigma_lag1_200_500km"]
    # sigma ratio against range, for the figure
    res["r_km"] = r.tolist()
    res["sigma_ratio_by_gate"] = n.nanmedian(n.sqrt(n.array([R[a]["var"][:, 1] / B[a]["var"][:, 1]
                                                             for a in common])), axis=0).tolist()
    res["columns_by_gate"] = hit_all.sum(axis=0).tolist()
    res["zs_lags1_5"] = n.nanmean(zs[:, 0:5], axis=1).tolist()
    return res


def fit_change(d_base, d_run, periods=PERIODS):
    tr = (fq.T0_RECORDING + fq.PERIOD_S * min(periods), fq.T0_RECORDING + fq.PERIOD_S * (max(periods) + 1))
    F0, F1 = fq.load_fits(d_base, tr), fq.load_fits(d_run, tr)
    if F0 is None or F1 is None:
        return {}
    t = sorted(set(F0["t0"]) & set(F1["t0"]))
    i0 = [list(F0["t0"]).index(x) for x in t]
    i1 = [list(F1["t0"]).index(x) for x in t]
    r = F0["r"]
    g0, g1 = fq.good(fq.flags(F0))[i0], fq.good(fq.flags(F1))[i1]
    b = (r >= 150) & (r <= 800)
    out = dict(windows=len(t), good_base=float(g0[:, b].mean()), good_run=float(g1[:, b].mean()))
    both = g0 & g1
    for lo, hi in BANDS:
        m = (r >= lo) & (r < hi)
        row = {}
        for p in ("ne", "Te", "Ti", "vi"):
            dz = (F1[p][i1][:, m] - F0[p][i0][:, m]) / F0["d" + p][i0][:, m]
            ds = F1["d" + p][i1][:, m] / F0["d" + p][i0][:, m]
            ok = both[:, m]
            row["shift_" + p] = float(n.nanmedian(dz[ok])) if ok.any() else n.nan
            row["mean_shift_" + p] = float(n.nanmean(dz[ok])) if ok.any() else n.nan
            row["err_ratio_" + p] = float(n.nanmedian(ds[ok])) if ok.any() else n.nan
        out["%d-%d" % (lo, hi)] = row
    return out


def figure(results, fig_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))
    rates = sorted(results, key=lambda x: int(x))
    for rate in rates:
        m = results[rate]["lpi"]
        z = n.array(m["zs_lags1_5"])
        ax[0].hist(z, bins=n.linspace(-3, 3, 41), histtype="step", density=True,
                   label="%s per period: mean %+.3f" % (rate, m["shift_lags1_5_mean"]))
        ax[1].plot(m["r_km"], m["sigma_ratio_by_gate"], label="%s per period" % rate)
    xx = n.linspace(-3, 3, 200)
    ax[0].set_xlabel(r"ACF shift at a column's gate, lags 1-5 ($\sigma$)")
    ax[0].set_ylabel("density")
    ax[0].set_title("Bias: shift at gates with random columns")
    ax[0].legend(fontsize=7)
    ax[1].axhline(1, color="k", lw=0.5)
    ax[1].set_xlim(100, 900)
    ax[1].set_xlabel("range (km)")
    ax[1].set_ylabel(r"$\sigma$ with / $\sigma$ without (lag 1)")
    ax[1].set_title("Precision: error bar growth")
    ax[1].legend(fontsize=7)
    x = [results[r]["lpi"]["columns_per_period"] for r in rates]
    y = [results[r]["lpi"]["sigma_ratio_at_columns"] for r in rates]
    y2 = [results[r]["lpi"]["sigma_ratio_all_150_800"] for r in rates]
    ax[2].plot(x, y, "o-", label="at the gates with columns")
    ax[2].plot(x, y2, "s-", label="all gates 150-800 km")
    ax[2].set_xscale("log")
    ax[2].axhline(1, color="k", lw=0.5)
    ax[2].set_xlabel("random columns per 10 s period")
    ax[2].set_ylabel(r"median $\sigma$ ratio (lag 1)")
    ax[2].set_title("Price against the number of columns")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "random_control.png"), dpi=150)
    plt.close(fig)


def main():
    base, out, fig_dir = sys.argv[1:4]
    runs = dict(a.split("=", 1) for a in sys.argv[4:])
    os.makedirs(out, exist_ok=True)
    B = read(base)
    results = {}
    for rate, d in runs.items():
        R = read(d)
        results[rate] = dict(lpi=measure(B, R), fits=fit_change(base, d))
        m, f = results[rate]["lpi"], results[rate]["fits"]
        print("rate %s: %.0f columns/period, %d gates; shift lags1-5 %+.3f +- %.3f, all lags %+.3f; "
              "z_diff mean %+.3f robust std %.2f (std %.2f); sigma ratio at columns %.3f, all %.4f; "
              "var growth per column %.4f (mean %.4f); honesty %.3f -> %.3f" % (
                  rate, m["columns_per_period"], m["gates_with_columns"], m["shift_lags1_5_mean"],
                  m["shift_lags1_5_err"], m["shift_all_mean"], m["zdiff_mean"], m["zdiff_robust_std"],
                  m["zdiff_std"], m["sigma_ratio_at_columns"], m["sigma_ratio_all_150_800"],
                  m["var_growth_per_column"], m["var_growth_per_column_mean"], m["honesty_base"],
                  m["honesty_run"]))
        if f:
            print("   fits %d windows: good %.3f -> %.3f" % (f["windows"], f["good_base"], f["good_run"]))
            for lo, hi in BANDS:
                b = f["%d-%d" % (lo, hi)]
                print("   %d-%d km: shift ne %+.2f Te %+.2f Ti %+.2f vi %+.2f (means %+.2f %+.2f %+.2f %+.2f);"
                      " error ratio ne %.3f Te %.3f" % (lo, hi, b["shift_ne"], b["shift_Te"], b["shift_Ti"],
                                                        b["shift_vi"], b["mean_shift_ne"], b["mean_shift_Te"],
                                                        b["mean_shift_Ti"], b["mean_shift_vi"],
                                                        b["err_ratio_ne"], b["err_ratio_Te"]))
    json.dump(results, open(os.path.join(out, "random_control.json"), "w"), indent=1)
    figure(results, fig_dir)


if __name__ == "__main__":
    main()
