#!/usr/bin/env python3
"""
The LPI variants compared at the level of the fitted plasma parameters.

1. Injection test (inject_satellites.py): for each variant, the fits of the
   injected run against the fits of the same variant without the injection,
   at the gate nearest the synthetic satellite (450 km): the shift of n_e,
   T_e, T_i and v_i in the clean fit's standard deviation, the mean over the
   fit windows.
2. Real intervals: for each variant and test interval, the fraction of good
   gates at 150-700 km (fit_quality.good), the median relative error of n_e,
   and the shift of n_e, T_e and T_i from variant 5 (variance weights,
   background and satellite columns) in variant 5's standard deviation, at
   the Starlink shells (535-590 km) and at 200-500 km.

Usage:
    python3 compare_fits.py <injection_dir> <variants_dir> <intervals.json> <echoes.npy> <out_dir>
"""

import glob
import json
import os
import sys

import numpy as n

import fit_quality as fq

PARAMS = ("ne", "Te", "Ti", "vi")


def injection(inj_root, clean_root):
    rows = []
    for d in sorted(glob.glob(os.path.join(inj_root, "*", "v*"))):
        meta = json.load(open(os.path.join(d, "injection.json")))
        per = [p["period"] for p in meta["periods"]]
        tr = (fq.T0_RECORDING + fq.PERIOD_S * min(per), fq.T0_RECORDING + fq.PERIOD_S * (max(per) + 1))
        Fi = fq.load_fits(d, tr)
        Fc = fq.load_fits(os.path.join(clean_root, meta["variant"]), tr)
        if Fi is None or Fc is None:
            continue
        g = int(n.argmin(n.abs(Fc["r"] - meta["range_km"])))
        row = dict(scenario=meta["scenario"], variant=meta["variant"], windows=len(Fi["t0"]),
                   gate_km=float(Fc["r"][g]))
        for p in PARAMS:
            z = (Fi[p][:, g] - Fc[p][:, g]) / Fc["d" + p][:, g]
            row["shift_" + p] = float(n.nanmean(z))
        rows.append(row)
    return rows


def real(var_root, intervals, echoes):
    rows = []
    ref = os.path.join(var_root, "v5_variance_background_satellites")
    for name, per in intervals.items():
        tr = (fq.T0_RECORDING + fq.PERIOD_S * min(per), fq.T0_RECORDING + fq.PERIOD_S * (max(per) + 1))
        R = fq.load_fits(ref, tr)
        for d in sorted(glob.glob(os.path.join(var_root, "v*"))):
            F = fq.load_fits(d, tr)
            if F is None:
                continue
            f = fq.flags(F, echoes)
            band = (F["r"] >= 150) & (F["r"] <= 700)
            ok = fq.good(f)
            row = dict(interval=name, variant=os.path.basename(d), windows=len(F["t0"]),
                       good_frac=float(ok[:, band].mean()),
                       median_rel_dne=float(n.nanmedian((F["dne"] / F["ne"])[:, band][ok[:, band]])))
            for lab, (a, b) in (("shells", (535, 590)), ("200_500", (200, 500))):
                sel = (F["r"] >= a) & (F["r"] <= b)
                both = ok[:, sel] & fq.good(fq.flags(R))[:, sel] if R is not None and len(R["t0"]) == len(F["t0"]) else None
                for p in ("ne", "Te", "Ti"):
                    if both is None or not both.any():
                        row["%s_shift_vs_v5_%s" % (p, lab)] = n.nan
                        continue
                    z = ((F[p][:, sel] - R[p][:, sel]) / R["d" + p][:, sel])[both]
                    row["%s_shift_vs_v5_%s" % (p, lab)] = float(n.nanmedian(z))
            rows.append(row)
    return rows


def main():
    inj_root, var_root, ivf, ecf, outd = sys.argv[1:6]
    os.makedirs(outd, exist_ok=True)
    out = dict(injection=injection(inj_root, var_root),
               real=real(var_root, json.load(open(ivf)), n.load(ecf)))
    json.dump(out, open(os.path.join(outd, "fit_measures.json"), "w"), indent=1)
    for r in out["injection"]:
        print("inj  %-6s %-34s %s" % (r["scenario"], r["variant"],
              "  ".join("%s %+.1f" % (p, r["shift_" + p]) for p in PARAMS)))
    for r in out["real"]:
        print("real %-18s %-34s good %.2f  rel dne %.2f  ne@shells %+.2f  Te@shells %+.2f" % (
            r["interval"], r["variant"], r["good_frac"], r["median_rel_dne"],
            r["ne_shift_vs_v5_shells"], r["Te_shift_vs_v5_shells"]))


if __name__ == "__main__":
    main()
