#!/usr/bin/env python3
"""
How much did the faulty corner weights of il_interp.ilint.getspec move the
fitted plasma parameters?

Before the fix, two corner weights of the 16-corner interpolation reused the
names of the Te/Ti axis weights (w10, w11) and overwrote them, and corner 0111
read the upper n_e table row. fit_lpi.py evaluates the model at n_e=1e12, a
table node, so all the weight sits on the upper n_e corners, the ones built
from the overwritten weights.

For a few integration periods of the archived lpi_30 products this fits every
range gate twice with fit_lpi.fit_lpifiles, unchanged: once with the fixed
interpolation and once with getspec from OLD_COMMIT. The topside O+ fraction
(fit_acf_ts, above 700 km) is recorded here by watching the fit's
minimizations, also where the fit's covariance fails; fit_lpifiles writes it
to the pp files as O_frac since fit-ts-composition. Every fit is its own
process, all at once.

    python3 validate_il_interp.py --out /some/dir
    python3 validate_il_interp.py --out /some/dir --report
"""

import os
for _v in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ[_v] = "1"
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import argparse
import glob
import subprocess
import sys
import types

import numpy as n
import h5py

DATA = "/mnt/data/juha/millstone_hill/isr/eclipse2024/usrp-rx0-r_20240407T100000_20240409T110000"
CHANNEL = "zenith-l"
POSTFIX = "_30"
# the commit before the fix
OLD_COMMIT = "b35df55"
# first file of a zenith-l block: 2024-04-08 16:01:50 UT (local noon),
# 04:12:50 UT (night) and 19:28:40 UT (eclipse maximum at Millstone Hill)
PERIODS = [1712592110, 1712549570, 1712604520]
FIT = dict(max_dt=300, first_lag=0, range_avg=n.array([1, 3, 5]), reanalyze=True, plot=False)
ARMS = ["old", "new"]


def lpi_files():
    fl = glob.glob("%s/lpi%s/%s/lpi-*.h5" % (DATA, POSTFIX, CHANNEL))
    return sorted(fl, key=lambda f: int(os.path.basename(f)[4:-3]))


def stage(out, t0, arm):
    """
    A directory holding the period's files only, plus the first file after it:
    fit_lpifiles closes a period when it meets the next one's first file.
    """
    base = "%s/%s_%d" % (out, arm, t0)
    d = "%s/lpi%s/%s" % (base, POSTFIX, CHANNEL)
    os.makedirs(d, exist_ok=True)
    fl = lpi_files()
    t = n.array([int(os.path.basename(f)[4:-3]) for f in fl])
    sel = n.where((t >= t0) & (t < t0 + FIT["max_dt"]))[0]
    sel = list(sel) + [n.where(t >= t0 + FIT["max_dt"])[0][0]]
    for i in sel:
        dst = "%s/%s" % (d, os.path.basename(fl[i]))
        if not os.path.exists(dst):
            os.symlink(fl[i], dst)
    return base


def old_ilint():
    """The ilint class as it was at OLD_COMMIT."""
    src = subprocess.run(["git", "show", "%s:il_interp.py" % OLD_COMMIT], check=True, capture_output=True,
                         text=True, cwd=os.path.dirname(os.path.abspath(__file__))).stdout
    mod = types.ModuleType("il_interp_old")
    # its table directory is found next to __file__
    mod.__file__ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "il_interp.py")
    exec(compile(src, "il_interp_old.py", "exec"), mod.__dict__)
    return mod.ilint


def run_one(out, t0, arm):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import fit_lpi
    base = stage(out, t0, arm)
    fit_lpi._init_tables(440.2e6)
    if arm == "old":
        old = old_ilint()
        fit_lpi.ilf.getspec = types.MethodType(old.getspec, fit_lpi.ilf)
        fit_lpi.ilf_ho.getspec = types.MethodType(old.getspec, fit_lpi.ilf_ho)

    # fit_acf_ts keeps the best of three minimizations and drops the fitted
    # composition, x[4]; repeat its choice to recover it
    runs = []
    so = fit_lpi.so

    def minimize(fun, x0, **kw):
        r = so.minimize(fun, x0, **kw)
        if len(x0) == 5:
            runs.append((r.x, fun(r.x)))
        return r
    fit_lpi.so = types.SimpleNamespace(minimize=minimize)
    fit_acf_ts = fit_lpi.fit_acf_ts
    comp = {}

    def fit_acf_ts_recorded(acf, lags, rgs, var, **kw):
        runs.clear()
        ok = False
        try:
            res = fit_acf_ts(acf, lags, rgs, var, **kw)
            ok = True
            return res
        finally:
            # the minimum is recorded also when the covariance that follows
            # fails, so that the failures can be looked at
            if len(runs) == 3:
                bx, sb = runs[0]
                for x, f in runs[1:]:
                    if f < sb:
                        bx = x
                comp[float(rgs)] = (bx, ok)
    fit_lpi.fit_acf_ts = fit_acf_ts_recorded

    fit_lpi.fit_lpifiles(dirn=base, channel=CHANNEL, postfix=POSTFIX, **FIT)
    k = sorted(comp)
    x = n.array([comp[h][0] for h in k]).reshape(-1, 5)
    n.savez("%s/composition.npz" % base, hgt=n.array(k), o_frac=1 - 1 / x[:, 4], x=x,
            fit_ok=n.array([comp[h][1] for h in k], dtype=bool))


def load(out, t0, arm):
    base = "%s/%s_%d" % (out, arm, t0)
    f = glob.glob("%s/lpi%s/%s/pp-*.h5" % (base, POSTFIX, CHANNEL))
    if not f:
        return None
    with h5py.File(f[0], "r") as h:
        r = {k: h[k][()] for k in ("Te", "Ti", "ne", "dTe_Ti", "dTi", "dne", "rgs")}
    r["Te_Ti"] = r["Te"] / r["Ti"]
    c = n.load("%s/composition.npz" % base)
    # only the gates whose fit completed, as for the other parameters
    r["o_hgt"], r["o_frac"] = c["hgt"][c["fit_ok"]], c["o_frac"][c["fit_ok"]]
    r["ts_hgt"], r["ts_x"], r["ts_ok"] = c["hgt"], c["x"], c["fit_ok"]
    return r


BANDS = [(0, 200), (200, 700), (700, 1500), (0, 1500)]


def report(out):
    lines = ["Fit with the fixed interpolation (new) against the faulty one (old), zenith-l, lpi_30.",
             "Per height band: median and largest |new-old|, and median |new-old|/sigma(new),",
             "over gates where both fits are finite and sigma(Ti) < 500 K.", ""]
    for t0 in PERIODS:
        o, w = load(out, t0, "old"), load(out, t0, "new")
        if o is None or w is None:
            lines.append("%d missing" % t0)
            continue
        lines.append("period starting %d (%s UT)" % (t0, n.datetime64(t0, "s")))
        h = w["rgs"]
        ok = n.isfinite(o["Ti"]) & n.isfinite(w["Ti"]) & (w["dTi"] < 500)
        for lo, hi in BANDS:
            m = ok & (h >= lo) & (h < hi)
            if not n.any(m):
                continue
            txt = "  %4d-%4d km %3d gates" % (lo, hi, n.sum(m))
            for k, dk, fmt in (("Te_Ti", "dTe_Ti", "%.3f"), ("Ti", "dTi", "%.0f K"), ("ne", "dne", "%.1f%%")):
                d = w[k][m] - o[k][m]
                z = n.abs(d) / w[dk][m]
                if k == "ne":
                    d = 100 * d / w[k][m]
                txt += (" | %s med " + fmt + " max " + fmt + " med %.2f sig") % (
                    k, n.median(n.abs(d)), n.max(n.abs(d)), n.nanmedian(z))
            lines.append(txt)
        # gates whose fit failed (NaN), and fitted Te/Ti sitting on a table
        # node: the faulty weights made the model flat in Te/Ti between nodes
        for arm, r in (("old", o), ("new", w)):
            fin = n.isfinite(r["Ti"]) & (h < 1500)
            node = n.min(n.abs(r["Te_Ti"][fin, None] - n.linspace(1, 3, 6)[None, 1:-1]), axis=1) < 2e-3
            txt = "  %s: %d gates below 1500 km fitted, %d failed (%s); Te/Ti within 0.002 of an inner table node: %d (%.0f%%)" % (
                arm, n.sum(fin), n.sum(~fin & (h < 1500) & (h > 80)),
                ", ".join("%d-%d km %d" % (lo, hi, n.sum(~n.isfinite(r["Ti"]) & (h >= max(lo, 80)) & (h < hi)))
                          for lo, hi in BANDS[:3]),
                n.sum(node), 100.0 * n.mean(node))
            lines.append(txt)
        # topside fits that ended on the Te/Ti bound of 3, the table's edge.
        # there the index is clamped and the finite difference step in Te/Ti
        # changes nothing, so the Jacobian loses a column and its inverse fails
        for arm, r in (("old", o), ("new", w)):
            at3 = r["ts_x"][:, 0] > 2.999
            lines.append("  %s topside: %d fits, %d at Te/Ti=3, of which %d failed; %d failed otherwise" % (
                arm, len(at3), n.sum(at3), n.sum(at3 & ~r["ts_ok"]), n.sum(~at3 & ~r["ts_ok"])))
        # the O+ fraction fitted above 700 km
        common = n.intersect1d(o["o_hgt"], w["o_hgt"])
        if len(common):
            fo = o["o_frac"][n.isin(o["o_hgt"], common)]
            fw = w["o_frac"][n.isin(w["o_hgt"], common)]
            lines.append("  O+ fraction above 700 km, %d gates: old median %.3f, new median %.3f, "
                         "median |new-old| %.3f, max %.3f" % (len(common), n.median(fo), n.median(fw),
                                                               n.median(n.abs(fw - fo)), n.max(n.abs(fw - fo))))
        lines.append("")
    txt = "\n".join(lines)
    print(txt)
    with open("%s/report.txt" % out, "w") as f:
        f.write(txt + "\n")


def plot(out, t0, fname):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    o, w = load(out, t0, "old"), load(out, t0, "new")
    h = w["rgs"]
    ok = n.isfinite(o["Ti"]) & n.isfinite(w["Ti"]) & (w["dTi"] < 500) & (h < 1500)
    fig, ax = plt.subplots(1, 4, figsize=(13, 4.5), sharey=True, layout="constrained")
    for a, k, lab in ((ax[0], "Te_Ti", r"$T_e/T_i$"), (ax[1], "Ti", r"$T_i$ (K)"), (ax[2], "Te", r"$T_e$ (K)")):
        a.plot(o[k][ok], h[ok], "o", mfc="none", ms=6, color="tab:red", label="old")
        a.plot(w[k][ok], h[ok], ".", color="tab:blue", label="fixed")
        a.set_xlabel(lab)
    ax[3].plot(o["o_frac"], o["o_hgt"], "o", mfc="none", ms=6, color="tab:red", label="old")
    ax[3].plot(w["o_frac"], w["o_hgt"], ".", color="tab:blue", label="fixed")
    ax[3].set_xlabel(r"O$^+$ fraction")
    ax[0].set_ylabel("Height (km)")
    ax[0].legend()
    fig.suptitle("zenith-l, period from %s UT" % n.datetime64(t0, "s"))
    fig.savefig(fname, dpi=150)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    p.add_argument("--report", action="store_true")
    p.add_argument("--plot", help="figure of the first period")
    p.add_argument("--one", nargs=2)
    args = p.parse_args()
    if args.one:
        run_one(args.out, int(args.one[0]), args.one[1])
        return
    if args.plot:
        plot(args.out, PERIODS[0], args.plot)
        return
    if not args.report:
        os.makedirs(args.out, exist_ok=True)
        running = []
        for t0 in PERIODS:
            for arm in ARMS:
                log = open("%s/%s_%d.log" % (args.out, arm, t0), "w")
                running.append(subprocess.Popen([sys.executable, os.path.abspath(__file__), "--out", args.out,
                                                 "--one", str(t0), arm], stdout=log, stderr=subprocess.STDOUT))
        for r in running:
            r.wait()
    report(args.out)


if __name__ == "__main__":
    main()
