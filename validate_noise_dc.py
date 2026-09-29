#!/usr/bin/env python3
"""
Test of the DC offset estimate in the noise injection calibration of
outlier_lpi.py (noise_dc="period" against the earlier "window").

Memo 7 found that subtracting each pulse's own background window mean biases
T_sys low by a factor that grows as the passband narrows and as T_sys itself
grows, and derived the bias with no free parameter (compare_tsys_sky.py,
lpi_bias_factor). A flat spectrum source must give the same T_sys at every
bandwidth, so the test is:

    1. with noise_dc="period", T_sys agrees across bandwidths;
    2. with noise_dc="window", it spreads as memo 7 found;
    3. the "window" values, corrected by lpi_bias_factor, equal the "period" ones;
    4. the ACFs are bit identical between the two, as only the calibration changes.

Two periods on 2024-04-08, zenith-l: the Cygnus A transit through the beam
(11:36:20 UTC), where T_sys is ~1600 K and the bias largest, and a baseline
1.9 h earlier. Each inversion is its own process, all run at once.

    python3 validate_noise_dc.py [--out plots/noise_dc_test]
"""

import os
for _v in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ[_v] = "1"
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import argparse
import glob
import json
import subprocess
import sys

import numpy as n
import h5py

DATA = "/mnt/data/juha/millstone_hill/isr/eclipse2024/usrp-rx0-r_20240407T100000_20240409T110000"
CHANNEL = "zenith-l"
PERIODS = {"transit": 9218, "baseline": 8519}
PASS_BANDS = {"18kHz": 0.018e6, "50kHz": 0.05e6, "full": 0.4e6}
MODES = ["window", "period"]
# as the archived lpi_30 products, apart from the gate, which T_sys does not depend on
LPI = dict(avg_dur=10, rg=60, min_tx_frac=0.2, filter_len=100, maximum_range_delay=7200,
           save_acf_images=False, lag_avg=1, reanalyze=True)


def run_one(out, period, band, mode):
    import outlier_lpi as olpi
    olpi.lpi_files(dirname=DATA, channel=CHANNEL, output_base="%s/%s_%s_%s" % (out, period, band, mode),
                   periods=[PERIODS[period]], pass_band=PASS_BANDS[band], noise_dc=mode, **LPI)


def load(out, period, band, mode):
    f = glob.glob("%s/%s_%s_%s/lpi_%d/%s/*.h5" % (out, period, band, mode, LPI["rg"], CHANNEL))[0]
    with h5py.File(f, "r") as h:
        return {k: h[k][()] for k in h.keys()}


def report(out):
    from compare_tsys_sky import lpi_bias_factor
    lines = ["T_sys (K), zenith-l 2024-04-08; 'corrected' is the window value through memo 7's lpi_bias_factor", ""]
    res = {}
    for period in PERIODS:
        lines.append("%s (period %d)" % (period, PERIODS[period]))
        lines.append("  %-6s %10s %10s %12s %14s" % ("band", "window", "period", "corrected", "period/corr."))
        tp = []
        for band, pb in PASS_BANDS.items():
            w, p = load(out, period, band, "window"), load(out, period, band, "period")
            corr = float(lpi_bias_factor(float(w["T_sys"]), pb))
            lines.append("  %-6s %10.1f %10.1f %12.1f %14.4f" % (band, w["T_sys"], p["T_sys"], corr, p["T_sys"] / corr))
            same = all(n.array_equal(w[k], p[k], equal_nan=True) for k in ["acfs_e", "acfs_g", "acfs_var", "noise_e", "noise_g"])
            res["%s_%s" % (period, band)] = dict(window=float(w["T_sys"]), period=float(p["T_sys"]), corrected=corr,
                                                  alpha_ratio=float(p["alpha"] / w["alpha"]), acfs_identical=bool(same))
            tp.append(float(p["T_sys"]))
        tw = [res["%s_%s" % (period, b)]["window"] for b in PASS_BANDS]
        lines.append("  spread across bandwidths: window %.1f %%, period %.1f %%"
                     % (100 * (max(tw) - min(tw)) / n.mean(tw), 100 * (max(tp) - min(tp)) / n.mean(tp)))
        lines.append("")
    lines.append("ACFs bit identical between the two modes in every run: %s"
                 % all(v["acfs_identical"] for v in res.values()))
    lines.append("alpha, period over window: %s" % ", ".join("%s %.4f" % (k, v["alpha_ratio"]) for k, v in res.items()))
    txt = "\n".join(lines)
    print(txt)
    open("%s/report.txt" % out, "w").write(txt + "\n")
    json.dump(res, open("%s/results.json" % out, "w"), indent=1)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default="plots/noise_dc_test")
    p.add_argument("--one", nargs=3, metavar=("PERIOD", "BAND", "MODE"), help=argparse.SUPPRESS)
    p.add_argument("--report", action="store_true", help="only report on finished runs")
    args = p.parse_args()
    if args.one:
        run_one(args.out, *args.one)
        return
    if not args.report:
        os.makedirs(args.out, exist_ok=True)
        procs = []
        for period in PERIODS:
            for band in PASS_BANDS:
                for mode in MODES:
                    log = open("%s/log_%s_%s_%s.txt" % (args.out, period, band, mode), "w")
                    procs.append(subprocess.Popen([sys.executable, __file__, "--out", args.out, "--one", period, band, mode],
                                                  stdout=log, stderr=subprocess.STDOUT))
        print("%d inversions running in parallel" % len(procs), flush=True)
        failed = [pr.args for pr in procs if pr.wait() != 0]
        if failed:
            print("failed:", failed)
    report(args.out)


if __name__ == "__main__":
    main()
