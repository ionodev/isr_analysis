#!/usr/bin/env python3
"""
The change table of Memos 45 and 48: how the fitted plasma parameters of one
run ("after") differ from another ("before") on the same fit windows and
range gates, by altitude band, 11-23 UTC of one day.

Per band, over all (window, gate) pairs where both runs have a value:

    ne_ratio    median of n_e after / n_e before
    dTe, dTi    median of after - before (K)
    abs_dvi     median of |v_i after - v_i before| (m/s)
    moved_2sig  share of the gates whose n_e or T_e moved by more than
                2 sigma, sigma being the "before" run's stated 1-sigma error
                at that gate (Memo 45's definition: it reproduces that
                table). Both runs use the same voltages, so their errors are
                not independent and are not added in quadrature.

and at the Starlink shells (535-590 km), for each run, the share of fitted
gates with T_e below 1000 K and the median T_e of the good gates
(fit_quality.good), as in Memo 45's table.

    python3 change_table.py <before_dir> <after_dir> <YYYY-MM-DD> [out.json]
"""
import datetime
import json
import sys

import numpy as n

import fit_quality as fq

BANDS = [(150, 250), (250, 350), (350, 450), (450, 520), (520, 600), (600, 700)]


def aligned(before, after, day):
    day0 = datetime.datetime.fromisoformat(day).replace(tzinfo=datetime.timezone.utc).timestamp()
    tr = (day0 + 11 * 3600, day0 + 23 * 3600)
    A, B = fq.load_fits(before, tr), fq.load_fits(after, tr)
    if not n.array_equal(A["r"], B["r"]):
        raise SystemExit("the two runs do not have the same range gates")
    # the fit windows both runs have
    common, ia, ib = n.intersect1d(A["t0"], B["t0"], return_indices=True)
    A = {k: (v[ia] if k != "r" else v) for k, v in A.items()}
    B = {k: (v[ib] if k != "r" else v) for k, v in B.items()}
    return A, B, len(common), len(set(A["t0"])), len(set(B["t0"]))


def table(A, B):
    r = A["r"]
    with n.errstate(invalid="ignore", divide="ignore"):
        ratio = B["ne"] / A["ne"]
        dTe = B["Te"] - A["Te"]
        dTi = B["Ti"] - A["Ti"]
        dvi = n.abs(B["vi"] - A["vi"])
        sne = A["dne"]
        sTe = A["dTe"]
        moved = (n.abs(B["ne"] - A["ne"]) > 2 * sne) | (n.abs(dTe) > 2 * sTe)
    both = n.isfinite(ratio) & n.isfinite(dTe)
    rows = []
    for lo, hi in BANDS:
        g = both & ((r >= lo) & (r < hi))[None, :]
        rows.append(dict(band="%d-%d" % (lo, hi), gates=int(g.sum()),
                         ne_ratio=float(n.nanmedian(ratio[g])), dTe=float(n.nanmedian(dTe[g])),
                         dTi=float(n.nanmedian(dTi[g])), abs_dvi=float(n.nanmedian(dvi[g])),
                         moved_2sig=float(n.mean(moved[g]))))
    return rows


def shells(F):
    m = (F["r"] >= 535) & (F["r"] <= 590)
    te = F["Te"][:, m]
    ok = n.isfinite(te)
    good = fq.good(fq.flags(F))[:, m]
    return dict(shell_Te_below_1000=float(n.mean(te[ok] < 1000)), shell_median_Te_good=float(n.nanmedian(te[good])))


def main():
    before, after, day = sys.argv[1:4]
    A, B, nc, na, nb = aligned(before, after, day)
    res = dict(day=day, before=before, after=after, windows_common=nc, windows_before=na, windows_after=nb,
               bands=table(A, B), shells_before=shells(A), shells_after=shells(B))
    print("%s: %d common fit windows (before %d, after %d)" % (day, nc, na, nb))
    for x in res["bands"]:
        print("  %-8s ne %.2f  dTe %+6.0f K  dTi %+5.0f K  |dvi| %4.0f m/s  >2sigma %4.1f %%  (%d gates)" % (
            x["band"], x["ne_ratio"], x["dTe"], x["dTi"], x["abs_dvi"], 100 * x["moved_2sig"], x["gates"]))
    for k in ("before", "after"):
        s = res["shells_" + k]
        print("  shells %-6s Te<1000 K %.3f, median Te (good gates) %.0f K" % (k, s["shell_Te_below_1000"], s["shell_median_Te_good"]))
    if len(sys.argv) > 4:
        json.dump(res, open(sys.argv[4], "w"), indent=1)


if __name__ == "__main__":
    main()
