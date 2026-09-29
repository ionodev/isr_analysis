#!/usr/bin/env python3
"""
Does the line interference of memo 20 reach the lag profile inversion's ACFs?

The interference survey (documents/analysis: survey_misa_rfi_psd.py,
analyse_misa_rfi.py) classifies the range-spread pulses; class N are narrow
lines from an outside transmitter.  For the integration periods where those
lines add the most power after the inversion's filter, and for periods where
they add a moderate amount, each period is inverted

    default   as the pipeline runs, with its outlier rejection
    norej     without outlier rejection
    masked    without the class N pulses (and without outlier rejection),
              the reference for norej
    masked_rej  without the class N pulses, with outlier rejection, the
              reference for default: the rejection also acts on other
              signals (satellites), which this cancels

at the eclipse configuration's filter (pass_band 100 kHz, filter_len 20) and
the archived products' (18 kHz, 100).  The difference of default and norej
from masked, in units of masked's standard deviation, is the interference's
effect; T_sys and alpha too.  Every inversion is its own process, all at once.
Periods holding only MISA's low elevation mode 800 are left out: the inversion
does not handle that mode (the range-Doppler path does).

    python3 validate_rfi_lpi.py --classes .../results/misa_rfi_classes.npz
    python3 validate_rfi_lpi.py --report
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
CLASSES = os.path.expanduser("~/isr_project/documents/analysis/results/misa_rfi_classes.npz")
T0 = 1712484000
CONFIGS = {"f100": dict(pass_band=0.1e6, filter_len=20), "f18": dict(pass_band=0.018e6, filter_len=100)}
ARMS = ["default", "norej", "masked", "masked_rej"]
LPI = dict(avg_dur=10, rg=60, min_tx_frac=0.2, maximum_range_delay=7200, save_acf_images=False, lag_avg=1,
           reanalyze=True)
OWN = {"misa-l": 1.0, "zenith-l": -1.0}


def lpi_periods(npr):
    """Periods holding pulses of a code the inversion handles (not only mode 800)."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from raw_reader import pulse_index
    ix = pulse_index(DATA)
    per = ((ix["key"] / 1e6 - T0) // 10).astype(int)
    ok = n.bincount(per[ix["sweepid"] != 800], minlength=npr)[:npr] > 100
    return ok


def choose(classes, n_worst, n_moderate):
    z = n.load(classes)
    out = []
    for ch in ("misa", "zenith"):
        b = z["%s_p_lpi100" % ch]
        ok = n.isfinite(b) & lpi_periods(len(b))
        worst = [int(i) for i in n.argsort(n.where(ok, b, -1))[::-1][:n_worst[ch]]]
        mid = n.where(ok & (b > 0.01) & (b < 0.2))[0]
        mid = [int(i) for i in mid[n.linspace(0, len(mid) - 1, n_moderate[ch]).astype(int)]]
        out += [(ch + "-l", p, "worst", float(b[p])) for p in worst] + [(ch + "-l", p, "moderate", float(b[p])) for p in mid]
    return out


def masked_keys(classes, channel, period):
    z = n.load(classes)
    ch = channel.split("-")[0]
    k = z["keys"][z["%s_class" % ch] == "N"]
    t = k / 1e6
    return k[(t >= T0 + 10 * period - 1) & (t < T0 + 10 * period + 11)]


def run_one(out, channel, period, config, arm, classes):
    import outlier_lpi as olpi
    import millstone_radar_state as mrs
    kw = dict(LPI, **CONFIGS[config])
    base = "%s/%s_p%d_%s_%s" % (out, channel, period, config, arm)
    if arm.startswith("masked"):
        mask = set(int(k) for k in masked_keys(classes, channel, period))
        print("masking %d pulses" % len(mask))
        original = mrs.get_antenna_select
        away = -OWN[channel]            # report the pulse as sent on the other antenna, so it is skipped

        def masked_select(path):
            tx_ant, rx_ant = original(path)
            return (lambda k: away if int(k) in mask else tx_ant(k)), rx_ant
        olpi.mrs.get_antenna_select = masked_select
    olpi.lpi_files(dirname=DATA, channel=channel, output_base=base, periods=[period],
                   outlier_rejection=arm in ("default", "masked_rej"), **kw)


def load(out, channel, period, config, arm):
    f = glob.glob("%s/%s_p%d_%s_%s/lpi_%d/%s/*.h5" % (out, channel, period, config, arm, LPI["rg"], channel))
    if not f:
        return None
    with h5py.File(f[0], "r") as h:
        return {k: h[k][()] for k in ("acfs_e", "acfs_var", "T_sys", "alpha", "rgs_km", "lags")}


def report(out, chosen):
    lines = ["Effect of the class N line interference on the inversion; z = (arm - reference) / std(reference),",
             "reference: masked_rej for default, masked for norej; over all gates and lags with finite values;",
             "|z|_max, median |z|, mean z of the real part", ""]
    res = []
    for ch, p, kind, extra in chosen:
        for cfg in CONFIGS:
            refs = {"default": load(out, ch, p, cfg, "masked_rej"), "norej": load(out, ch, p, cfg, "masked")}
            if refs["norej"] is None:
                continue
            row = dict(channel=ch, period=p, kind=kind, extra_power_f100=extra, config=cfg)
            txt = "%-8s p%-5d %-8s extra %7.2f  %-4s" % (ch, p, kind, extra, cfg)
            for arm in ("default", "norej"):
                a = load(out, ch, p, cfg, arm)
                m = refs[arm]
                if m is None:
                    a = None
                if a is None:
                    txt += " | %s: missing" % arm
                    continue
                zz = (a["acfs_e"] - m["acfs_e"]) / n.sqrt(m["acfs_var"])
                zr = n.real(zz)[n.isfinite(zz)]
                row[arm] = dict(zmax=float(n.max(n.abs(zz[n.isfinite(zz)]))), zmed=float(n.median(n.abs(zr))),
                                zmean=float(n.mean(zr)), tsys_ratio=float(a["T_sys"] / m["T_sys"]),
                                alpha_ratio=float(a["alpha"] / m["alpha"]))
                txt += " | %s: |z|max %6.1f med %.2f mean %+.2f  Tsys x%.3f" % (
                    arm, row[arm]["zmax"], row[arm]["zmed"], row[arm]["zmean"], row[arm]["tsys_ratio"])
            lines.append(txt)
            res.append(row)
    txt = "\n".join(lines)
    print(txt)
    open("%s/report.txt" % out, "w").write(txt + "\n")
    json.dump(res, open("%s/results.json" % out, "w"), indent=1)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default="plots/rfi_lpi_test")
    p.add_argument("--classes", default=CLASSES)
    p.add_argument("--report", action="store_true")
    p.add_argument("--one", nargs=4, help=argparse.SUPPRESS)
    args = p.parse_args()
    if args.one:
        ch, per, cfg, arm = args.one
        run_one(args.out, ch, int(per), cfg, arm, args.classes)
        return
    os.makedirs(args.out, exist_ok=True)
    chosen = choose(args.classes, {"misa": 12, "zenith": 8}, {"misa": 12, "zenith": 8})
    json.dump(chosen, open("%s/chosen.json" % args.out, "w"), indent=1)
    if not args.report:
        jobs = [(ch, per, cfg, arm) for ch, per, _, _ in chosen for cfg in CONFIGS for arm in ARMS
                if load(args.out, ch, per, cfg, arm) is None]
        ncpu = os.cpu_count()
        print("%d inversions, %d at a time" % (len(jobs), ncpu), flush=True)
        running = []
        for j in jobs:
            while len(running) >= ncpu:
                running = [r for r in running if r.poll() is None]
                if len(running) >= ncpu:
                    running[0].wait()
            log = open("%s/log_%s_p%d_%s_%s.txt" % ((args.out,) + j), "w")
            running.append(subprocess.Popen([sys.executable, __file__, "--out", args.out, "--classes", args.classes,
                                             "--one"] + [str(x) for x in j], stdout=log, stderr=subprocess.STDOUT))
        for r in running:
            r.wait()
    report(args.out, chosen)


if __name__ == "__main__":
    main()
