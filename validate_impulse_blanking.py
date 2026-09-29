#!/usr/bin/env python3
"""
Effect of the power-line impulses (memo 21) on the lag profile inversion, and
of blanking them (impulse_blanking.BlankingReader, lpi_files blank_impulses).

Periods are taken from the 60 s windows of the power-line survey
(documents/analysis/results/powerline_survey.npz) in which the receiver's own
antenna transmits, at 400 kW or more, a mode the inversion handles: for each
receiver the 8 windows with the most impulses and the 4 with the fewest, two
periods from each.  Every period is inverted as the pipeline runs (default) and
with the impulses blanked, at the eclipse configuration's filter (100 kHz,
filter_len 20) and the archived products' (18 kHz, 100).  The difference
default - blanked is given in units of the blanked inversion's standard
deviation, over all gates and lags and as the mean over gates per lag, with the
fraction of samples blanked.  Every inversion is its own process, all at once.

    python3 validate_impulse_blanking.py [--out plots/impulse_blanking_test]
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
SURVEY = os.path.expanduser("~/isr_project/documents/analysis/results/powerline_survey.npz")
T0 = 1712484000
CONFIGS = {"f100": dict(pass_band=0.1e6, filter_len=20), "f18": dict(pass_band=0.018e6, filter_len=100)}
ARMS = ("default", "blanked")
LPI = dict(avg_dur=10, rg=60, min_tx_frac=0.2, maximum_range_delay=7200, save_acf_images=False, lag_avg=1,
           reanalyze=True)
OWN = {"misa-l": 1.0, "zenith-l": -1.0}


def choose(n_most=8, n_least=4):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import millstone_radar_state as mrs
    from raw_reader import pulse_index
    z = n.load(SURVEY, allow_pickle=True)
    ix = pulse_index(DATA)
    tx_ant, rx_ant = mrs.get_antenna_select("%s/metadata/antenna_control_metadata" % DATA)
    zpm, mpm = mrs.get_tx_power_model("%s/metadata/powermeter" % DATA)
    out = []
    for ch in OWN:
        m = n.where(z["channel"] == ch)[0]
        ok = []
        for i in m:
            t = z["t0"][i] + 30
            k = t * 1e6
            pw = mpm(t) if OWN[ch] > 0 else zpm(t)
            a, b = n.searchsorted(ix["key"], [z["t0"][i] * 1e6, (z["t0"][i] + 60) * 1e6])
            if tx_ant(k) * OWN[ch] > 0.99 and rx_ant(k) * OWN[ch] > 0.99 and pw >= 400e3 and \
                    n.mean(ix["sweepid"][a:b] == 800) < 0.5:
                ok.append(i)
        ok = n.array(ok)
        rate = z["n_used"][ok] / z["n_pulses"][ok]
        o = n.argsort(rate)
        for kind, sel in (("impulsive", o[::-1][:n_most]), ("quiet", o[:n_least])):
            for i, r in zip(ok[sel], rate[sel]):
                for d in (10, 30):
                    out.append((ch, int((z["t0"][i] + d - T0) // 10), kind, float(r)))
    return out


def run_one(out, ch, period, cfg, arm):
    import outlier_lpi as olpi
    readers = []

    class Counting(olpi.BlankingReader):
        def __init__(self, *a, **k):
            super().__init__(*a, **k)
            readers.append(self)
    olpi.BlankingReader = Counting
    olpi.lpi_files(dirname=DATA, channel=ch, output_base="%s/%s_p%d_%s_%s" % (out, ch, period, cfg, arm),
                   periods=[period], blank_impulses=(arm == "blanked"), **dict(LPI, **CONFIGS[cfg]))
    for r in readers:
        print("blanked fraction %.6f of %d samples" % (r.n_blanked / max(r.n_samples, 1), r.n_samples))


def load(out, ch, period, cfg, arm):
    f = glob.glob("%s/%s_p%d_%s_%s/lpi_%d/%s/*.h5" % (out, ch, period, cfg, arm, LPI["rg"], ch))
    if not f:
        return None
    with h5py.File(f[0], "r") as h:
        return {k: h[k][()] for k in ("acfs_e", "acfs_var", "noise_e", "T_sys", "alpha", "lags", "rgs_km")}


def report(out, chosen):
    lines = ["default - blanked, in units of the blanked inversion's standard deviation", ""]
    res = []
    per_lag = {}
    for ch, p, kind, rate in chosen:
        for cfg in CONFIGS:
            a, b = load(out, ch, p, cfg, "default"), load(out, ch, p, cfg, "blanked")
            if a is None or b is None:
                lines.append("%-8s p%-5d %-9s %-4s missing" % (ch, p, kind, cfg))
                continue
            z = (a["acfs_e"] - b["acfs_e"]) / n.sqrt(b["acfs_var"])
            zr = n.real(z)
            ok = n.isfinite(zr)
            lag_mean = n.nanmean(zr, axis=0)
            per_lag.setdefault((ch, kind, cfg), []).append(lag_mean)
            blk = open("%s/log_%s_p%d_%s_blanked.txt" % (out, ch, p, cfg)).read()
            row = dict(channel=ch, period=p, kind=kind, rate=rate, config=cfg,
                       zmax=float(n.nanmax(n.abs(z))), zmed=float(n.nanmedian(n.abs(zr[ok]))),
                       zmean=float(n.mean(zr[ok])), zmean_lag1=float(lag_mean[0]),
                       tsys_ratio=float(a["T_sys"] / b["T_sys"]),
                       blanked=float(blk.split("blanked fraction")[-1].split()[0]) if "blanked fraction" in blk else n.nan)
            res.append(row)
            lines.append("%-8s p%-5d %-9s impulse rate %.2f  %-4s | |z|max %6.1f  med %.2f  mean %+.3f  first lag mean %+.2f"
                         "  Tsys x%.4f  blanked %.5f" % (ch, p, kind, rate, cfg, row["zmax"], row["zmed"], row["zmean"],
                                                         row["zmean_lag1"], row["tsys_ratio"], row["blanked"]))
    lines.append("")
    lines.append("mean over gates and periods of the real part of z, per lag (us):")
    for (ch, kind, cfg), v in sorted(per_lag.items()):
        v = n.nanmean(n.array(v), axis=0)
        lines.append("  %-8s %-9s %-4s %s" % (ch, kind, cfg, " ".join("%+.2f" % x for x in v[:12])))
    txt = "\n".join(lines)
    print(txt)
    open("%s/report.txt" % out, "w").write(txt + "\n")
    json.dump(res, open("%s/results.json" % out, "w"), indent=1)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default="plots/impulse_blanking_test")
    p.add_argument("--report", action="store_true")
    p.add_argument("--one", nargs=4, help=argparse.SUPPRESS)
    args = p.parse_args()
    if args.one:
        ch, per, cfg, arm = args.one
        run_one(args.out, ch, int(per), cfg, arm)
        return
    os.makedirs(args.out, exist_ok=True)
    chosen = choose()
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
            running.append(subprocess.Popen([sys.executable, __file__, "--out", args.out, "--one"] + [str(x) for x in j],
                                            stdout=log, stderr=subprocess.STDOUT))
        for r in running:
            r.wait()
    report(args.out, chosen)


if __name__ == "__main__":
    main()
