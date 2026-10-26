#!/usr/bin/env python3
"""
The measures of Memos 42 and 45 for several runs on the test intervals:
precision, honesty, echo residual at the catalogued echoes of at least
20 dB (from the default catalogue, P_FA 1e-10, the same echoes for every
run), the fraction of good fitted gates, and at the Starlink shells
(535-590 km) the fraction of fitted gates with T_e below 1000 K.

    python3 interval_measures.py <intervals.json> <out.json> name=dir [name=dir ...]
"""
import json
import sys

import numpy as n
from digital_rf import DigitalMetadataReader

import compare_variants as cv
import fit_quality as fq

CAT = "/mnt/data/henrik/satellite_detections_2026-10-02"


def echoes(periods):
    rd = DigitalMetadataReader(CAT)
    rows = []
    for a in periods:
        i0 = int((cv.T0_RECORDING + 10 * a) * 1e6)
        for k, rec in rd.read(i0, i0 + 10_000_000).items():
            ch = rec["channel"]
            ch = ch.decode() if isinstance(ch, bytes) else str(ch)
            if ch != "zenith-l":
                continue
            for r, s in zip(rec["range_km"], rec["snr_db"]):
                rows.append((a, float(r), float(s)))
    return n.array(rows) if rows else n.zeros((0, 3))


def shells(d, periods):
    tr = (cv.T0_RECORDING + 10 * min(periods), cv.T0_RECORDING + 10 * (max(periods) + 1))
    F = fq.load_fits(d, tr)
    if F is None:
        return {}
    m = (F["r"] >= 535) & (F["r"] <= 590)
    te = F["Te"][:, m]
    g = fq.good(fq.flags(F))
    b = (F["r"] >= 150) & (F["r"] <= 700)
    ok = n.isfinite(te)
    return dict(shell_Te_below_1000=float(n.mean(te[ok] < 1000)) if ok.any() else n.nan,
                shell_median_Te=float(n.nanmedian(te)), good_150_700=float(g[:, b].mean()),
                fit_windows=int(len(F["t0"])))


def main():
    iv = json.load(open(sys.argv[1]))
    out = sys.argv[2]
    runs = dict(a.split("=", 1) for a in sys.argv[3:])
    res = []
    for name, per in iv.items():
        E = echoes(per)
        for rn, d in runs.items():
            L = cv.read_lpi(d, per)
            if not L:
                continue
            m = cv.measures(L, per, E)
            m.update(shells(d, per))
            m.update(run=rn, interval=name)
            res.append(m)
            print("%-20s %-12s periods %2d sigma %.0f honesty %.3f echoes %d residual |z| %.2f (>3: %.3f) "
                  "good %.3f shell Te<1000 %.3f median Te %.0f" % (
                      name, rn, m["periods"], m["sigma_lag1_200_500km"], m["honesty"], m["echoes_tested"],
                      m["residual_median_absz"], m["residual_frac_above_3"], m.get("good_150_700", n.nan),
                      m.get("shell_Te_below_1000", n.nan), m.get("shell_median_Te", n.nan)), flush=True)
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
