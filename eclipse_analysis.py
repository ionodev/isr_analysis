#!/usr/bin/env python3
"""
The eclipse analysis of Memo 52: 8 April 2024 (the eclipse) against 7 April
(the control day), 11-23 UTC, zenith, from the 60 s fits of one run, using
only the good gates (fit_quality.good: not failed, not at a fit limit, not
noisy).

Per fit window:
    NmF2, hmF2   the F-region peak: a parabola in altitude fitted to ln n_e
                 at the five good gates around the largest n_e at 200-450 km,
                 weighted by (n_e / dn_e)^2; NmF2 = exp(peak), hmF2 its
                 altitude, with errors from the fit's covariance
    Te300, Ti300 the error-weighted mean of the good gates at 270-330 km
    Te560        the same for T_e at 530-590 km, the Starlink shells, where
                 the satellite handling matters
Both days are then binned in 10-minute bins of UTC (median of the windows in
the bin; error from the windows' scatter, or the fits' errors for a single
window), and compared bin by bin: the ratio NmF2(8 April) / NmF2(7 April),
the difference of hmF2, Te300 and Ti300. The pre-eclipse hours (11-18 UTC)
give the day-to-day baseline: the mean and scatter of the same comparison
when there is no eclipse.

    python3 eclipse_analysis.py <out.json> <name>=<fits_0407>,<fits_0408> [...]
(a run covering both days is given twice). With --plot <prefix>, figures for
the first run: the time series with the obscuration (eclipse_timing.py) and
the profiles before, at and after maximum.
"""
import datetime
import json
import sys

import numpy as n

import fit_quality as fq

DAYS = ("2024-04-07", "2024-04-08")
BIN_S = 600
PRE = (11.0, 18.0)          # hours UTC: the baseline, before the eclipse
MAX_GROUND_H = 19 + 29.5 / 60   # eclipse_timing.py, ground
MAX_300_H = 19 + 34 / 60 + 10 / 3600


def day0(day):
    return datetime.datetime.fromisoformat(day).replace(tzinfo=datetime.timezone.utc).timestamp()


def per_window(F):
    r = F["r"]
    good = fq.good(fq.flags(F))
    out = []
    for i in range(len(F["t0"])):
        ne, dne = F["ne"][i], F["dne"][i]
        ok = good[i] & n.isfinite(ne) & (ne > 0) & n.isfinite(dne) & (dne > 0)
        rec = dict(t=0.5 * (F["t0"][i] + F["t1"][i]), NmF2=n.nan, dNmF2=n.nan, hmF2=n.nan, dhmF2=n.nan)
        cand = n.where(ok & (r >= 200) & (r <= 450))[0]
        if len(cand) >= 5:
            k = cand[n.argmax(ne[cand])]
            idx = [j for j in range(k - 2, k + 3) if 0 <= j < len(r) and ok[j]]
            if len(idx) >= 4:
                x = r[idx] - r[k]
                y = n.log(ne[idx])
                w = (ne[idx] / dne[idx]) ** 2
                A = n.vstack([x ** 2, x, n.ones_like(x)]).T
                try:
                    C = n.linalg.inv(A.T @ (A * w[:, None]))
                    a, b, c = C @ (A.T @ (w * y))
                    if a < 0:
                        x0 = -b / (2 * a)
                        if abs(x0) <= 2.5 * (r[1] - r[0]):
                            peak = c - b * b / (4 * a)
                            # gradients of x0 and peak with respect to (a, b, c)
                            gx = n.array([b / (2 * a * a), -1 / (2 * a), 0.0])
                            gp = n.array([b * b / (4 * a * a), -b / (2 * a), 1.0])
                            rec.update(NmF2=float(n.exp(peak)), dNmF2=float(n.exp(peak) * n.sqrt(gp @ C @ gp)),
                                       hmF2=float(r[k] + x0), dhmF2=float(n.sqrt(gx @ C @ gx)))
                except n.linalg.LinAlgError:
                    pass
        for key, tag, lo, hi in (("Te", "300", 270, 330), ("Ti", "300", 270, 330), ("Te", "560", 530, 590)):
            m = ok & (r >= lo) & (r <= hi) & n.isfinite(F[key][i]) & (F["d" + key][i] > 0)
            if m.any():
                wt = 1 / F["d" + key][i][m] ** 2
                rec[key + tag] = float(n.sum(wt * F[key][i][m]) / n.sum(wt))
                rec["d" + key + tag] = float(1 / n.sqrt(n.sum(wt)))
            else:
                rec[key + tag] = rec["d" + key + tag] = n.nan
        out.append(rec)
    return out


def binned(rows, d0, key):
    t = n.array([x["t"] for x in rows])
    v = n.array([x[key] for x in rows])
    e = n.array([x["d" + key] for x in rows])
    edges = n.arange(d0 + 11 * 3600, d0 + 23 * 3600 + 1, BIN_S)
    hv, he, hh = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (t >= a) & (t < b) & n.isfinite(v)
        hh.append((0.5 * (a + b) - d0) / 3600)
        if m.sum() == 0:
            hv.append(n.nan); he.append(n.nan)
        elif m.sum() == 1:
            hv.append(v[m][0]); he.append(e[m][0])
        else:
            hv.append(n.median(v[m]))
            he.append(max(n.std(v[m], ddof=1) / n.sqrt(m.sum()), n.sqrt(n.mean(e[m] ** 2) / m.sum())))
    return n.array(hh), n.array(hv), n.array(he)


def compare(runs_rows):
    """runs_rows: {day: rows}. Bin-by-bin comparison and the key numbers."""
    res = {}
    h, n7, e7 = binned(runs_rows[DAYS[0]], day0(DAYS[0]), "NmF2")
    _, n8, e8 = binned(runs_rows[DAYS[1]], day0(DAYS[1]), "NmF2")
    ratio = n8 / n7
    eratio = ratio * n.sqrt((e8 / n8) ** 2 + (e7 / n7) ** 2)
    pre = (h >= PRE[0]) & (h < PRE[1]) & n.isfinite(ratio)
    base, base_sd = n.mean(ratio[pre]), n.std(ratio[pre], ddof=1)
    rel = ratio / base
    ok = n.isfinite(rel) & (h >= 18) & (h <= 22)
    i = n.where(ok)[0][n.argmin(rel[ok])]
    res["NmF2"] = dict(hours=h.tolist(), ratio=ratio.tolist(), ratio_err=eratio.tolist(),
                       baseline=float(base), baseline_sd=float(base_sd),
                       min_rel=float(rel[i]), min_rel_err=float(n.hypot(eratio[i] / base, base_sd / base)),
                       min_hour=float(h[i]), lag_after_max_ground_min=float((h[i] - MAX_GROUND_H) * 60))
    for key in ("hmF2", "Te300", "Ti300", "Te560"):
        _, a7, s7 = binned(runs_rows[DAYS[0]], day0(DAYS[0]), key)
        _, a8, s8 = binned(runs_rows[DAYS[1]], day0(DAYS[1]), key)
        d = a8 - a7
        ed = n.hypot(s7, s8)
        pre = (h >= PRE[0]) & (h < PRE[1]) & n.isfinite(d)
        base, base_sd = n.mean(d[pre]), n.std(d[pre], ddof=1)
        rel = d - base
        ok = n.isfinite(rel) & (h >= 18) & (h <= 22)
        j = n.where(ok)[0]
        res[key] = dict(hours=h.tolist(), day7=a7.tolist(), day8=a8.tolist(), diff=d.tolist(), diff_err=ed.tolist(),
                        baseline=float(base), baseline_sd=float(base_sd))
        # the lowest and the highest change against the baseline, 18-22 UTC
        for tag, k in (("min", j[n.argmin(rel[j])]), ("max", j[n.argmax(rel[j])])):
            res[key].update({tag: float(rel[k]), tag + "_err": float(n.hypot(ed[k], base_sd)), tag + "_hour": float(h[k]),
                             tag + "_lag_min": float((h[k] - MAX_GROUND_H) * 60)})
    return res


def main():
    args = sys.argv[1:]
    plot = None
    if "--plot" in args:
        i = args.index("--plot"); plot = args[i + 1]; args = args[:i] + args[i + 2:]
    out = args[0]
    allres = {}
    first = None
    for spec in args[1:]:
        name, dirs = spec.split("=", 1)
        d7, d8 = dirs.split(",")
        rows = {}
        for day, d in zip(DAYS, (d7, d8)):
            F = fq.load_fits(d, (day0(day) + 11 * 3600, day0(day) + 23 * 3600))
            rows[day] = per_window(F)
        res = compare(rows)
        res["windows"] = rows
        allres[name] = res
        first = first or (name, rows, res)
        x = res["NmF2"]
        print("%-12s NmF2(8)/NmF2(7): baseline %.3f +- %.3f (11-18 UTC); minimum %.3f +- %.3f of baseline at %.2f h UTC "
              "(%+.0f min after maximum at the ground)" % (name, x["baseline"], x["baseline_sd"], x["min_rel"],
                                                          x["min_rel_err"], x["min_hour"], x["lag_after_max_ground_min"]))
        for key, unit in (("hmF2", "km"), ("Te300", "K"), ("Ti300", "K"), ("Te560", "K")):
            x = res[key]
            print("%-12s %s(8)-(7): baseline %+.0f +- %.0f %s; lowest %+.0f +- %.0f at %.2f h (%+.0f min), highest %+.0f +- %.0f at %.2f h (%+.0f min)" % (
                name, key, x["baseline"], x["baseline_sd"], unit, x["min"], x["min_err"], x["min_hour"], x["min_lag_min"],
                x["max"], x["max_err"], x["max_hour"], x["max_lag_min"]))
    json.dump(allres, open(out, "w"), indent=1, default=float)
    if plot:
        make_plots(plot, *first)


def make_plots(prefix, name, rows, res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        import eclipse_timing as et
        from astropy.time import Time
        import astropy.units as u
        tt = Time("2024-04-08T17:00:00", scale="utc") + n.arange(0, 5 * 3600, 60) * u.s
        obs, _ = et.obscuration(tt, 300.0)
        th = 17 + n.arange(0, 5 * 3600, 60) / 3600
    except Exception:
        obs = None
    fig, ax = plt.subplots(4, 1, figsize=(10, 11), sharex=True)
    for day, c in zip(DAYS, ("tab:blue", "tab:red")):
        d0 = day0(day)
        t = (n.array([x["t"] for x in rows[day]]) - d0) / 3600
        for a, key, sc in ((ax[0], "NmF2", 1e-12), (ax[1], "hmF2", 1), (ax[2], "Te300", 1), (ax[3], "Ti300", 1)):
            v = n.array([x[key] for x in rows[day]]) * sc
            a.plot(t, v, ".", color=c, ms=2, alpha=0.4)
            h, bv, be = binned(rows[day], d0, key)
            a.errorbar(h, bv * sc, be * sc, fmt="-", color=c, lw=1.2, label="%s (10 min medians)" % day)
    for a in ax:
        a.axvspan(18 + 16 / 60, 20 + 39 / 60, color="0.9", zorder=0)
        a.axvline(MAX_GROUND_H, color="k", ls="--", lw=0.8)
    if obs is not None:
        ax0b = ax[0].twinx()
        ax0b.plot(th, obs, color="0.4", lw=1)
        ax0b.set_ylabel("obscuration at 300 km", color="0.4")
        ax0b.set_ylim(0, 1.05)
    ax[0].set_ylabel("NmF2 ($10^{12}$ m$^{-3}$)")
    ax[1].set_ylabel("hmF2 (km)")
    ax[2].set_ylabel("$T_e$ at 270-330 km (K)")
    ax[3].set_ylabel("$T_i$ at 270-330 km (K)")
    ax[3].set_xlabel("hours (UTC)")
    ax[0].legend(fontsize=8, loc="lower left")
    ax[0].set_title("Zenith, %s: grey band = partial eclipse at Millstone Hill, dashed = maximum" % name)
    fig.tight_layout()
    fig.savefig(prefix + "_timeseries.png", dpi=110)


if __name__ == "__main__":
    main()
