#!/usr/bin/env python3
"""
What the outlier rejection removes when it runs together with tracked
satellite columns: per interval, the fraction of rows each test removes (the
ratio test, the power cut, either; outlier_*_fraction in the LPI files, mean
over periods and lags), and how the rows it removes are spread over the
satellite columns by the class of column_checks (coherent, decaying, no
echo, shared), the classes taken from the run without rejection (the same
columns: tracked columns do not depend on the rejection).

    python3 rejection_measures.py <intervals.json> <run_with_rejection> <run_without> <out.json>
"""
import json
import os
import sys

import h5py
import numpy as n

import column_checks as cc
import compare_variants as cv


def main():
    iv = json.load(open(sys.argv[1]))
    d_on, d_off, out = sys.argv[2], sys.argv[3], sys.argv[4]
    res = {}
    for name, per in iv.items():
        fr = {"ratio": [], "power": [], "any": []}
        by_class = {c: [] for c in ("coherent", "decaying", "no echo", "shared")}
        n_cols = 0
        for a in per:
            f_on, f_off = cv.lpi_file(d_on, a), cv.lpi_file(d_off, a)
            if not (os.path.exists(f_on) and os.path.exists(f_off)):
                continue
            with h5py.File(f_on, "r") as h:
                if "outlier_any_fraction" not in h:
                    continue
                for k in fr:
                    fr[k].append(h["outlier_%s_fraction" % k][()])
                rej = n.nanmean(h["sat_rows_rejected"][()], axis=1) if h["sat_keys"].shape[0] else n.zeros(0)
                keys_on = h["sat_keys"][()]
                del_on = h["sat_delay_samples"][()]
            r = cc.one_file((f_off, None))
            if r is None or len(rej) == 0:
                continue
            cls = cc.classify(r["S"], r["rho"], r["srho"], r["nearest"], cc.S1_THEORY)
            lut = {(int(k), round(float(d), 4)): c for k, d, c in zip(r["key"], r["delay"], cls)}
            for k, d, q in zip(keys_on, del_on, rej):
                c = lut.get((int(k), round(float(d), 4)))
                if c is not None and n.isfinite(q):
                    by_class[c].append(q)
                    n_cols += 1
        if not fr["any"]:
            continue
        res[name] = dict(periods=len(fr["any"]), columns=n_cols,
                         **{"rows_removed_" + k: float(n.nanmean(v)) for k, v in fr.items()},
                         **{"col_rows_rejected_" + c: (float(n.mean(v)) if v else n.nan) for c, v in by_class.items()},
                         **{"cols_" + c: len(v) for c, v in by_class.items()})
        m = res[name]
        print("%-20s periods %d rows removed: ratio %.4f power %.4f either %.4f | columns %d; mean fraction of a "
              "column's rows removed: coherent %.3f (%d) decaying %.3f (%d) no echo %.3f (%d) shared %.3f (%d)" % (
                  name, m["periods"], m["rows_removed_ratio"], m["rows_removed_power"], m["rows_removed_any"],
                  n_cols, m["col_rows_rejected_coherent"], m["cols_coherent"], m["col_rows_rejected_decaying"],
                  m["cols_decaying"], m["col_rows_rejected_no echo"], m["cols_no echo"],
                  m["col_rows_rejected_shared"], m["cols_shared"]), flush=True)
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
