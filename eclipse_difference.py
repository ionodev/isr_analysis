#!/usr/bin/env python3
"""
Where the satellite handling changed the plasma parameters: an altitude-time
plot of one run minus another (for example the tracked satellite columns
minus no satellite handling), in physical units, drawn like eclipse_plots.py
--raw: every fitted value as it is, no flags, no gaps.

    n_e   the ratio, as log10(n_e after / n_e before): negative where the
          satellite's power had raised the density
    T_e, T_i   after - before, in K: positive where the satellite had made
          the plasma look colder
    v_i   after - before, in m/s

The two runs must have the same fit windows and range gates (the same
periods and fit settings).

Usage:
    python3 eclipse_difference.py <before_fits_dir> <after_fits_dir> <day YYYY-MM-DD> <out.png> [<label_before> <label_after>]
"""

import datetime
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as n

import eclipse_plots as ep
import fit_quality as fq


def main():
    before, after, day, out = sys.argv[1:5]
    lb = sys.argv[5] if len(sys.argv) > 5 else "before"
    la = sys.argv[6] if len(sys.argv) > 6 else "after"
    ep.MAX_COLUMN_S = n.inf
    day0 = datetime.datetime.fromisoformat(day).replace(tzinfo=datetime.timezone.utc).timestamp()
    tr = (day0 + 11 * 3600, day0 + 23 * 3600)
    A, B = fq.load_fits(before, tr), fq.load_fits(after, tr)
    ia, ib = n.argsort(A["t0"]), n.argsort(B["t0"])
    if not (n.array_equal(A["t0"][ia], B["t0"][ib]) and n.array_equal(A["r"], B["r"])):
        raise SystemExit("the two runs do not have the same fit windows and gates")
    A = {k: (v[ia] if k != "r" else v) for k, v in A.items()}
    B = {k: (v[ib] if k != "r" else v) for k, v in B.items()}
    r = A["r"]
    redges = n.concatenate([r - 0.5 * (r[1] - r[0]), [r[-1] + 0.5 * (r[1] - r[0])]])
    with n.errstate(invalid="ignore", divide="ignore"):
        d = {"ne": n.log10(B["ne"] / A["ne"]), "Te": B["Te"] - A["Te"], "Ti": B["Ti"] - A["Ti"], "vi": B["vi"] - A["vi"]}
    G = {k: ep.grid(A, v) for k, v in d.items()}
    tg = G["ne"][0]
    fig, ax = plt.subplots(4, 1, figsize=(11, 12), sharex=True)
    ep.panel(ax[0], tg, redges, G["ne"][1], day0, "RdBu_r", -0.3, 0.3, "log$_{10}$($n_e$ ratio)")
    ep.panel(ax[1], tg, redges, G["Te"][1], day0, "RdBu_r", -2000, 2000, "$\\Delta T_e$ (K)")
    ep.panel(ax[2], tg, redges, G["Ti"][1], day0, "RdBu_r", -1000, 1000, "$\\Delta T_i$ (K)")
    ep.panel(ax[3], tg, redges, G["vi"][1], day0, "RdBu_r", -150, 150, "$\\Delta v_i$ (m/s)")
    ax[0].set_title("%s, zenith: %s minus %s; dashed = eclipse" % (day, la, lb))
    ax[-1].set_xlabel("hours (UTC)")
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    s = (r >= 535) & (r <= 590)
    o = (r >= 200) & (r <= 450)
    print("535-590 km: median dTe %+.0f K, dTi %+.0f K, ne ratio %.2f; 200-450 km: median dTe %+.0f K, |dTe| median %.0f K" % (
        n.nanmedian(d["Te"][:, s]), n.nanmedian(d["Ti"][:, s]), 10 ** n.nanmedian(d["ne"][:, s]),
        n.nanmedian(d["Te"][:, o]), n.nanmedian(n.abs(d["Te"][:, o]))))


if __name__ == "__main__":
    main()
