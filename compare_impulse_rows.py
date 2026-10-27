#!/usr/bin/env python3
"""
The RFI mitigation by removing the impulse-hit rows (impulse_rows) measured
on minutes with power-line sparks and on the quiet reference interval.

Two runs of the same LPI variant (variant 5: variance weights, background
and satellite columns), without and with "impulse_rows", on the same
periods. Per interval:

    removed     the fraction of rows removed per lag (data cost) and the
                fraction of echo-window samples found to be impulses
    shift       the ACF of the default run minus that with the rows removed,
                in the latter's standard deviation of the real part, mean
                over 150-600 km per lag: the impulse power the default lets
                into the ACF (its bias), and on the quiet interval nothing
    sigma       the median ratio of the reported standard deviations
                (with/without removal) at 200-500 km, per lag: the cost (or
                gain) in precision
    honesty     the scatter between neighbouring periods over the reported
                error (compare_variants), with and without
    background  the background column's estimate (noise_e) per lag, which
                absorbs any range-independent power such as spread impulses
    fits        the 60 s fits: the fraction of good gates at 150-700 km, and
                the median change of n_e, T_e, T_i, v_i in the default fit's
                standard deviation, by altitude

Usage:
    python3 compare_impulse_rows.py <run_off> <run_on> <intervals.json> <out_dir> [<figure_dir>]
"""

import json
import os
import sys

import h5py
import numpy as n

import compare_variants as cv
import fit_quality as fq

BANDS = [(150, 250), (250, 350), (350, 450), (450, 600), (600, 700)]


def read(d, periods):
    out = {}
    for ai in periods:
        f = cv.lpi_file(d, ai)
        if not os.path.exists(f):
            continue
        with h5py.File(f, "r") as h:
            x = dict(acf=h["acfs_e"][()], var=h["acfs_var"][()], r=h["rgs_km"][()],
                     noise=h["noise_e"][()], lags=h["lags"][()])
            if "impulse_rows_removed_fraction" in h:
                x["removed"] = h["impulse_rows_removed_fraction"][()]
                x["imp_samples"] = float(h["impulse_sample_fraction"][()])
            out[ai] = x
    return out


def interval_measures(off, on):
    common = sorted(set(off) & set(on))
    r = off[common[0]]["r"]
    lags = off[common[0]]["lags"] * 1e6
    A0 = n.array([off[a]["acf"] for a in common])
    V0 = n.array([off[a]["var"] for a in common])
    A1 = n.array([on[a]["acf"] for a in common])
    V1 = n.array([on[a]["var"] for a in common])
    b = (r >= 150) & (r <= 600)
    s1 = n.sqrt(V1[:, b, :] / 2.0)
    z = (A0[:, b, :] - A1[:, b, :]).real / s1
    sig = (r >= 200) & (r <= 500)
    ratio = n.sqrt(V1[:, sig, :] / V0[:, sig, :])
    res = dict(periods=len(common), lags_us=lags.tolist(),
               removed_per_lag=n.mean([on[a]["removed"] for a in common], axis=0).tolist(),
               impulse_sample_fraction=float(n.mean([on[a]["imp_samples"] for a in common])),
               shift_per_lag=n.nanmean(z, axis=(0, 1)).tolist(),
               shift_err_per_lag=(n.nanstd(n.nanmean(z, axis=1), axis=0) / n.sqrt(len(common))).tolist(),
               sigma_ratio_per_lag=n.nanmedian(ratio, axis=(0, 1)).tolist(),
               noise_off=n.mean([off[a]["noise"] for a in common], axis=0).real.tolist(),
               noise_on=n.mean([on[a]["noise"] for a in common], axis=0).real.tolist())
    res["shift_lags1_5"] = float(n.nanmean(z[:, :, 0:5]))
    res["shift_lags6_45"] = float(n.nanmean(z[:, :, 5:]))
    res["sigma_ratio_lag1"] = float(res["sigma_ratio_per_lag"][1])
    res["sigma_ratio_all"] = float(n.nanmedian(ratio))
    res["removed_mean"] = float(n.mean(res["removed_per_lag"]))
    for name, L in (("off", off), ("on", on)):
        m = cv.measures({a: dict(acf=L[a]["acf"], var=L[a]["var"], r=L[a]["r"]) for a in common}, common,
                        n.zeros((0, 3)))
        res["honesty_" + name] = m["honesty"]
        res["sigma_lag1_" + name] = m["sigma_lag1_200_500km"]
    return res


def fit_change(d_off, d_on, periods):
    tr = (fq.T0_RECORDING + fq.PERIOD_S * min(periods), fq.T0_RECORDING + fq.PERIOD_S * (max(periods) + 1))
    F0, F1 = fq.load_fits(d_off, tr), fq.load_fits(d_on, tr)
    if F0 is None or F1 is None:
        return {}
    t = sorted(set(F0["t0"]) & set(F1["t0"]))
    i0 = [list(F0["t0"]).index(x) for x in t]
    i1 = [list(F1["t0"]).index(x) for x in t]
    r = F0["r"]
    g0, g1 = fq.good(fq.flags(F0))[i0], fq.good(fq.flags(F1))[i1]
    b = (r >= 150) & (r <= 700)
    out = dict(windows=len(t), good_off=float(g0[:, b].mean()), good_on=float(g1[:, b].mean()))
    both = g0 & g1
    for lo, hi in BANDS:
        m = (r >= lo) & (r < hi)
        row = {}
        for p in ("ne", "Te", "Ti", "vi"):
            dz = (F1[p][i1][:, m] - F0[p][i0][:, m]) / F0["d" + p][i0][:, m]
            ds = F1["d" + p][i1][:, m] / F0["d" + p][i0][:, m]
            ok = both[:, m]
            row["shift_" + p] = float(n.nanmedian(dz[ok])) if ok.any() else n.nan
            row["abs_shift_" + p] = float(n.nanmedian(n.abs(dz[ok]))) if ok.any() else n.nan
            row["err_ratio_" + p] = float(n.nanmedian(ds[ok])) if ok.any() else n.nan
        out["%d-%d" % (lo, hi)] = row
    return out


def figures(results, off_dir, on_dir, intervals, fig_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = list(results)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    for name in names:
        r = results[name]["lpi"]
        lab = name.replace("_", " ")
        ax[0].plot(r["lags_us"], 100 * n.array(r["removed_per_lag"]), label=lab)
        ax[1].errorbar(r["lags_us"], r["shift_per_lag"], yerr=r["shift_err_per_lag"], label=lab, capsize=0)
        ax[2].plot(r["lags_us"], r["sigma_ratio_per_lag"], label=lab)
    ax[0].set_ylabel("rows removed (%)")
    ax[0].set_title("Data cost")
    ax[1].axhline(0, color="k", lw=0.5)
    ax[1].set_ylabel(r"ACF without $-$ with removal ($\sigma$)")
    ax[1].set_title("Impulse power left in by default\n(mean over 150-600 km)")
    ax[2].axhline(1, color="k", lw=0.5)
    ax[2].set_ylabel(r"$\sigma$ with / $\sigma$ without removal")
    ax[2].set_title("Error bars (200-500 km)")
    for a in ax:
        a.set_xlabel(r"lag ($\mu$s)")
    ax[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "impulse_rows_measures.png"), dpi=150)
    plt.close(fig)


def example_figure(fig_dir, period, data_dir="/mnt/data/juha/millstone_hill/isr/eclipse2024/"
                   "usrp-rx0-r_20240407T100000_20240409T110000"):
    """One pulse of a spark period: the raw power with the impulses found, and
    the rows each lag loses."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from digital_rf import DigitalMetadataReader
    import impulse_rows as ir
    from outlier_lpi import simple_decimator
    from radar_timing import TMM
    from raw_reader import RawReader
    i0 = int((fq.T0_RECORDING + fq.PERIOD_S * period) * 1e6)
    sid = DigitalMetadataReader("%s/metadata/id_metadata" % data_dir).read(i0, i0 + 2000000, "sweepid")
    rd = RawReader("%s/rf_data/" % data_dir)
    best = None
    for k in sorted(sid)[3:200]:
        if int(sid[k]) not in TMM or TMM[int(sid[k])]["read_length"] > 10000:
            continue
        z = rd.read_vector_1d(k, 10000, "zenith-l").astype(n.complex64)
        t = TMM[int(sid[k])]
        hot = ir.detect(z, t["gc"], t["last_echo"])
        if best is None or 3 <= hot.sum() // 5 <= 12 and hot.sum() > best[2].sum():
            best = (k, z, hot, t)
    k, z, hot, t = best
    decim = simple_decimator(L=10000, dec=60)
    touched = ir.widen(hot, 20)
    lags = n.arange(1, 46) * 10
    H = n.array([ir.hit_rows(touched, lag, decim, 13, 137) for lag in lags])
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.8), gridspec_kw=dict(width_ratios=[1.4, 1]))
    p = n.abs(z) ** 2
    s = n.arange(len(z))
    ax[0].semilogy(s[t["gc"]:t["last_echo"]], p[t["gc"]:t["last_echo"]], lw=0.4, color="0.4")
    ax[0].semilogy(s[hot], p[hot], "r.", ms=4, label="impulse samples found (with guard)")
    ax[0].set_xlabel(r"sample in the pulse ($\mu$s after the pulse start)")
    ax[0].set_ylabel("raw power (arb.)")
    ax[0].set_title("One zenith-l pulse, %s UTC" % _utc(k))
    ax[0].legend(fontsize=8, loc="upper right")
    rows_us = (n.arange(13, 137)) * 60
    ax[1].pcolormesh(rows_us, lags, H, cmap="Greys", shading="auto")
    ax[1].set_xlabel(r"row start ($\mu$s)")
    ax[1].set_ylabel(r"lag ($\mu$s)")
    ax[1].set_title("Rows removed for this pulse (black): %.0f%%" % (100 * H.mean()))
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "impulse_rows_example.png"), dpi=150)
    plt.close(fig)
    return dict(pulse=int(k), impulses=int(n.sum(n.diff(n.concatenate([[0], hot.astype(int)])) == 1)),
                rows_removed=float(H.mean()))


def _utc(key):
    import datetime as dt
    return dt.datetime.fromtimestamp(key / 1e6, dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")


def fit_profile_figure(off_dir, on_dir, periods, name, fig_dir):
    """Median profiles of the 60 s fits of one interval, without and with removal."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tr = (fq.T0_RECORDING + fq.PERIOD_S * min(periods), fq.T0_RECORDING + fq.PERIOD_S * (max(periods) + 1))
    F0, F1 = fq.load_fits(off_dir, tr), fq.load_fits(on_dir, tr)
    r = F0["r"]
    fig, ax = plt.subplots(1, 4, figsize=(13, 4.2), sharey=True)
    for F, lab, c in ((F0, "default", "C0"), (F1, "impulse rows removed", "C3")):
        g = fq.good(fq.flags(F))
        for a, p, sc in zip(ax, ("ne", "Te", "Ti", "vi"), (1, 1, 1, 1)):
            v = n.where(g, F[p], n.nan)
            dv = n.where(g, F["d" + p], n.nan)
            med = n.nanmedian(v, axis=0)
            err = n.nanmedian(dv, axis=0) / n.sqrt(n.maximum(n.sum(n.isfinite(v), axis=0), 1))
            a.errorbar(med, r, xerr=err, color=c, label=lab, lw=1, capsize=0)
    ax[0].set_xscale("log")
    for a, lab in zip(ax, (r"$n_e$ (m$^{-3}$)", r"$T_e$ (K)", r"$T_i$ (K)", r"$v_i$ (m/s)")):
        a.set_xlabel(lab)
    ax[0].set_ylabel("range (km)")
    ax[0].set_ylim(120, 700)
    ax[0].legend(fontsize=8)
    fig.suptitle("Median of the good 60 s fits, %s" % name.replace("_", " "))
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "impulse_rows_fits_%s.png" % name), dpi=150)
    plt.close(fig)


def main():
    off_dir, on_dir, ivf, out = sys.argv[1:5]
    fig_dir = sys.argv[5] if len(sys.argv) > 5 else None
    intervals = json.load(open(ivf))
    os.makedirs(out, exist_ok=True)
    results = {}
    for name, per in intervals.items():
        off, on = read(off_dir, per), read(on_dir, per)
        if not off or not on:
            print(name, "missing")
            continue
        results[name] = dict(lpi=interval_measures(off, on), fits=fit_change(off_dir, on_dir, per))
        m = results[name]["lpi"]
        print("%-24s periods %d  impulse samples %.4f  rows removed %.3f (lag 1 %.3f)  shift lags 1-5 %+.3f  "
              "6-45 %+.3f  sigma ratio lag1 %.3f all %.3f  honesty %.2f -> %.2f" % (
                  name, m["periods"], m["impulse_sample_fraction"], m["removed_mean"], m["removed_per_lag"][0],
                  m["shift_lags1_5"], m["shift_lags6_45"], m["sigma_ratio_lag1"], m["sigma_ratio_all"],
                  m["honesty_off"], m["honesty_on"]))
        f = results[name]["fits"]
        if f:
            print("    fits: windows %d, good %.3f -> %.3f" % (f["windows"], f["good_off"], f["good_on"]))
            for lo, hi in BANDS:
                b = f["%d-%d" % (lo, hi)]
                print("    %d-%d km: median shift ne %+.2f Te %+.2f Ti %+.2f vi %+.2f sigma;"
                      " |shift| ne %.2f Te %.2f; error ratio ne %.3f Te %.3f" % (
                          lo, hi, b["shift_ne"], b["shift_Te"], b["shift_Ti"], b["shift_vi"],
                          b["abs_shift_ne"], b["abs_shift_Te"], b["err_ratio_ne"], b["err_ratio_Te"]))
    json.dump(results, open(os.path.join(out, "impulse_rows_measures.json"), "w"), indent=1)
    if fig_dir:
        figures(results, off_dir, on_dir, intervals, fig_dir)
        for name, per in intervals.items():
            if name.startswith("spark"):
                fit_profile_figure(off_dir, on_dir, per, name, fig_dir)
        ex = example_figure(fig_dir, intervals[[k for k in intervals if k.startswith("spark")][-1]][8])
        print("example pulse:", ex)
        json.dump(ex, open(os.path.join(out, "example_pulse.json"), "w"))


if __name__ == "__main__":
    main()
