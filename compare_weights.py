#!/usr/bin/env python3
"""
Variance weights against leave-one-out weights in the final product, by the
project plan's four measures (Henrik, 1 November: "compare leave one out
weights with variance weights").

The variance weights (variants 2-5) weight each lagged-product row by the
inverse of a running mean of |m|^2 over ten pulses. When that mean includes
the row itself, the weighted estimate is pulled towards zero by about
1/((1+SNR) N) (weight_bias_mc.py, Memo 51); the leave-one-out weights
(weights_exclude_self) leave the row out. On the quiet cycle the
leave-one-out ACF agrees with the unbiased equal-weight estimate to 0.4 %
(Memo 51), so here the leave-one-out run is the reference, and a difference
is the variance weights' bias.

Both runs are the final product (satellite removal v2, impulse-row removal)
on 7-8 April: the same voltages and fit windows. Each run's n_e is put on an
absolute scale with its own density constant from the plasma lines
(plasma_line_density_scale.py), so that the comparison shows the altitude
dependence the constant cannot absorb. Per altitude band and day:

  bias         median of (leave-one-out - variance) for T_e, T_i, v_i, and
               median of the ratio for n_e, over the (window, gate) pairs
               where both fits are good (fit_quality.good); also in units of
               the variance run's stated 1-sigma error
  variance     stated: median 1-sigma error of each run (relative for n_e);
               empirical: the scatter between consecutive fit windows on the
               zenith beam (which does not move, so the difference holds the
               noise plus the ionosphere's change in about a minute), 1.4826
               MAD of the difference / sqrt(2), over pairs at most PAIR_MAX_S
               apart
  retention    share of the gates whose fit is good
  residual     at the Starlink shells (535-590 km): share of the fitted
               gates with T_e < 1000 K, an impossible value that a satellite
               echo left in the data produces (Memos 45, 48)

Usage:
    python3 compare_weights.py <variance dir> <loo dir> <C variance> <C loo> <out dir>
The dirs contain "{day}", replaced by 0407 and 0408.
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as n

import fit_quality as fq
from change_table import shells

BANDS = [(150, 250), (250, 350), (350, 450), (450, 520), (520, 600), (600, 700)]
PROFILE_KM = n.arange(150, 701, 25)
PAIR_MAX_S = 90.0


def load_pair(dv, dl, Cv, Cl):
    V, L = fq.load_fits(dv, constant=Cv), fq.load_fits(dl, constant=Cl)
    if not n.array_equal(V["r"], L["r"]):
        raise SystemExit("the two runs do not have the same range gates")
    common, iv, il = n.intersect1d(V["t0"], L["t0"], return_indices=True)
    V = {k: (v[iv] if k != "r" else v) for k, v in V.items()}
    L = {k: (v[il] if k != "r" else v) for k, v in L.items()}
    return V, L


def robust_sd(x):
    x = x[n.isfinite(x)]
    return 1.4826 * n.median(n.abs(x - n.median(x))) if len(x) > 5 else n.nan


def successive_scatter(F, ok, p, s):
    """Empirical 1-sigma of p (relative for n_e) from consecutive windows."""
    pair = n.where(n.diff(F["t0"]) <= PAIR_MAX_S)[0]
    x, g = F[p][:, s], ok[:, s]
    d = []
    for k in pair:
        m = g[k] & g[k + 1]
        if p == "ne":
            d.append((x[k + 1][m] - x[k][m]) / (0.5 * (x[k + 1][m] + x[k][m])))
        else:
            d.append(x[k + 1][m] - x[k][m])
    d = n.concatenate(d) if d else n.zeros(0)
    return float(robust_sd(d) / n.sqrt(2)), len(d)


def stated(F, ok, p, s):
    e = n.abs(F["d" + p][:, s][ok[:, s]])
    if p == "ne":
        e = e / F["ne"][:, s][ok[:, s]]
    return float(n.median(e)) if len(e) else n.nan


def band_rows(V, L, okv, okl, bands):
    r = V["r"]
    both = okv & okl
    rows = []
    for lo, hi in bands:
        s = (r >= lo) & (r < hi)
        b = both[:, s]
        row = dict(band="%d-%d" % (lo, hi), gates_both=int(b.sum()))
        with n.errstate(invalid="ignore", divide="ignore"):
            q = (L["ne"][:, s] / V["ne"][:, s])[b]
            row["ne_ratio"] = float(n.median(q)) if len(q) else n.nan
            row["ne_ratio_16_84"] = [float(x) for x in n.percentile(q, [16, 84])] if len(q) else [n.nan, n.nan]
            for p in ("Te", "Ti", "vi"):
                d = (L[p] - V[p])[:, s][b]
                row["d" + p] = float(n.median(d)) if len(d) else n.nan
                row["d" + p + "_sigma"] = float(n.median(d / V["d" + p][:, s][b])) if len(d) else n.nan
        for name, F, ok in (("var", V, okv), ("loo", L, okl)):
            row["retention_" + name] = float(ok[:, s].mean())
            for p in ("ne", "Te"):
                row["stated_%s_%s" % (p, name)] = stated(F, ok, p, s)
                row["empirical_%s_%s" % (p, name)], row["pairs_" + name] = successive_scatter(F, ok, p, s)
        rows.append(row)
    return rows


def main():
    dv, dl, Cv, Cl, out = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
    os.makedirs(out, exist_ok=True)
    res = dict(variance=dv, loo=dl, C_variance=Cv, C_loo=Cl, days={})
    prof = {}
    for day in ("0407", "0408"):
        V, L = load_pair(dv.replace("{day}", day), dl.replace("{day}", day), Cv, Cl)
        okv, okl = fq.good(fq.flags(V)), fq.good(fq.flags(L))
        rows = band_rows(V, L, okv, okl, BANDS)
        res["days"][day] = dict(windows=len(V["t0"]), bands=rows, shells_var=shells(V), shells_loo=shells(L))
        prof[day] = band_rows(V, L, okv, okl, [(a, a + 25) for a in PROFILE_KM[:-1]])
        print("%s April: %d common 60-s fit windows" % (day[2:], len(V["t0"])))
        print("  band     ne L/V (16-84)      dTe K (sig)    dTi K (sig)   dvi m/s  good V/L    "
              "dne/ne stated V/L  empir. V/L    dTe stated V/L  empir. V/L")
        for x in rows:
            print("  %-8s %.3f (%.2f-%.2f) %+5.0f (%+.2f) %+5.0f (%+.2f) %+6.1f  %.2f/%.2f   "
                  "%.3f/%.3f   %.3f/%.3f   %4.0f/%4.0f   %4.0f/%4.0f" % (
                      x["band"], x["ne_ratio"], x["ne_ratio_16_84"][0], x["ne_ratio_16_84"][1],
                      x["dTe"], x["dTe_sigma"], x["dTi"], x["dTi_sigma"], x["dvi"],
                      x["retention_var"], x["retention_loo"],
                      x["stated_ne_var"], x["stated_ne_loo"], x["empirical_ne_var"], x["empirical_ne_loo"],
                      x["stated_Te_var"], x["stated_Te_loo"], x["empirical_Te_var"], x["empirical_Te_loo"]))
        for k in ("var", "loo"):
            s = res["days"][day]["shells_" + k]
            print("  shells %-3s: T_e < 1000 K in %.3f of the fitted gates, median T_e (good) %.0f K"
                  % (k, s["shell_Te_below_1000"], s["shell_median_Te_good"]))
    json.dump(res, open(os.path.join(out, "compare_weights.json"), "w"), indent=1)

    zc = PROFILE_KM[:-1] + 12.5
    fig, ax = plt.subplots(1, 4, figsize=(13, 5.2), sharey=True)
    for day, c in (("0407", "C0"), ("0408", "C3")):
        P = prof[day]
        lab = "%s April" % day[2:].lstrip("0")
        g = lambda k: n.array([x[k] for x in P])
        lo, hi = n.array([x["ne_ratio_16_84"] for x in P]).T
        ax[0].plot(g("ne_ratio"), zc, color=c, label=lab)
        ax[0].fill_betweenx(zc, lo, hi, color=c, alpha=0.15)
        ax[1].plot(g("dTe"), zc, color=c, label="$T_e$, " + lab)
        ax[1].plot(g("dTi"), zc, color=c, ls="--", label="$T_i$, " + lab)
        ax[2].plot(g("empirical_ne_loo") / g("empirical_ne_var"), zc, color=c, label="$n_e$, " + lab)
        ax[2].plot(g("empirical_Te_loo") / g("empirical_Te_var"), zc, color=c, ls="--", label="$T_e$, " + lab)
        ax[3].plot(g("retention_var"), zc, color=c, ls="--", label="variance, " + lab)
        ax[3].plot(g("retention_loo"), zc, color=c, label="leave-one-out, " + lab)
    ax[0].axvline(1, color="0.4", lw=0.8)
    ax[0].set_xlabel("$n_e$ leave-one-out / variance\n(each on its own plasma-line scale)")
    ax[0].set_ylabel("range (km)")
    ax[1].axvline(0, color="0.4", lw=0.8)
    ax[1].set_xlabel("leave-one-out $-$ variance (K)\nmedian over gates good in both")
    ax[2].axvline(1, color="0.4", lw=0.8)
    ax[2].set_xlabel("empirical scatter, leave-one-out / variance\n(consecutive 60-s windows)")
    ax[3].set_xlabel("share of gates with a good fit")
    for a in ax:
        a.legend(fontsize=7)
        a.grid(alpha=0.3)
    fig.suptitle("Final product, 7-8 April 2024, zenith: leave-one-out against variance weights", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "compare_weights.png"), dpi=130)


if __name__ == "__main__":
    main()
