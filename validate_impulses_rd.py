#!/usr/bin/env python3
"""
Do the power-line impulses of memo 21 reach the range-Doppler spectra of the
long pulse (mode 300, avg_range_doppler_spec.py)?

The periods of validate_impulse_blanking.py (for each receiver the most and
the least impulsive, with the receiver's own antenna transmitting) are
processed with the reader as recorded and with impulse_blanking.BlankingReader,
both averaged with the outlier test (the pipeline's default) and with a plain
mean.  The blanked and the unblanked calibration windows differ (the
injection window hides impulses the background shows), so the spectra are
compared as signal-to-noise ratios, (RDS - noise) / noise with the noise from
the range the path itself uses (700-1200 km for mode 300), which is what
fit_lp.py fits and does not depend on alpha.  The difference is given in units
of the blanked SNR's standard deviation.

    python3 validate_impulses_rd.py [--out plots/impulses_rd_test]
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
CHOSEN = os.path.expanduser("~/isr_project/isr_analysis_powerline/plots/impulse_blanking_test/chosen.json")
ARMS = {"outlier": ("outlier_removal", False), "outlier_blanked": ("outlier_removal", True),
        "mean": ("mean", False), "mean_blanked": ("mean", True)}
REF = {"outlier": "outlier_blanked", "mean": "mean_blanked"}
NOISE_KM = (700, 1200)


def run_one(out, ch, period, arm, tx_delay_us):
    import avg_range_doppler_spec as ards
    from raw_reader import RawReader
    from impulse_blanking import BlankingReader
    avg_type, blank = ARMS[arm]
    ards.DigitalRFReader = (lambda path: BlankingReader(RawReader(path), DATA)) if blank else RawReader
    ards.avg_range_doppler_spectra(dirname=DATA, channel=ch, mode=300, avg_dur=10, step=10, avg_type=avg_type,
                                   postfix="_" + arm, reanalyze=True, output_base="%s/%s_p%d" % (out, ch, period),
                                   tx_delay_us=tx_delay_us, periods=[period])


def load(out, ch, period, arm):
    f = glob.glob("%s/%s_p%d/range_doppler_300_%s/%s/il_*.h5" % (out, ch, period, arm, ch))
    if not f:
        return None
    with h5py.File(f[0], "r") as h:
        return {k: h[k][()] for k in ("RDS_LP", "RDS_LP_var", "alpha", "T_sys", "n_pulses", "rgs_km")}


def snr(d):
    r = (d["rgs_km"] > NOISE_KM[0]) & (d["rgs_km"] < NOISE_KM[1])
    noise = n.median(d["RDS_LP"][r])
    return (d["RDS_LP"] - noise) / noise, noise


def report(out, chosen):
    lines = ["mode 300 range-Doppler spectra, as recorded against impulses blanked; z = dSNR / sigma(SNR, blanked)", ""]
    res = []
    for ch, p, kind, rate in chosen:
        txt = "%-8s p%-5d %-9s impulse rate %.2f" % (ch, p, kind, rate)
        row = dict(channel=ch, period=p, kind=kind, rate=rate)
        for arm, rarm in REF.items():
            a, b = load(out, ch, p, arm), load(out, ch, p, rarm)
            if a is None or b is None:
                txt += " | %s: missing" % arm
                continue
            sa, na = snr(a)
            sb, nb = snr(b)
            sig = n.sqrt(b["RDS_LP_var"] / b["alpha"]) / nb
            z = (sa - sb) / sig
            ok = n.isfinite(z)
            low = b["rgs_km"] < 600
            row[arm] = dict(zmax=float(n.max(n.abs(z[ok]))), zmed=float(n.median(n.abs(z[ok]))), zmean=float(n.mean(z[ok])),
                            zmean_low=float(n.nanmean(z[low])), noise_ratio=float((na * a["alpha"]) / (nb * b["alpha"])),
                            tsys_ratio=float(a["T_sys"] / b["T_sys"]))
            r = row[arm]
            txt += " | %s: |z|max %6.1f med %.2f mean %+.3f (below 600 km %+.3f), raw noise x%.3f, T_sys x%.3f" % (
                arm, r["zmax"], r["zmed"], r["zmean"], r["zmean_low"], r["noise_ratio"], r["tsys_ratio"])
        lines.append(txt)
        res.append(row)
    txt = "\n".join(lines)
    print(txt)
    open("%s/report.txt" % out, "w").write(txt + "\n")
    json.dump(res, open("%s/results.json" % out, "w"), indent=1)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default="plots/impulses_rd_test")
    p.add_argument("--report", action="store_true")
    p.add_argument("--one", nargs=4, help=argparse.SUPPRESS)
    args = p.parse_args()
    if args.one:
        run_one(args.out, args.one[0], int(args.one[1]), args.one[2], float(args.one[3]))
        return
    os.makedirs(args.out, exist_ok=True)
    chosen = [tuple(c) for c in json.load(open(CHOSEN))]
    json.dump(chosen, open("%s/chosen.json" % args.out, "w"), indent=1)
    if not args.report:
        import millstone_radar_state as mrs
        import tx_delay as txd
        zpm, mpm = mrs.get_tx_power_model("%s/metadata/powermeter" % DATA)
        tx_ant, rx_ant = mrs.get_antenna_select("%s/metadata/antenna_control_metadata" % DATA)
        delay = {}
        for ch in sorted(set(c[0] for c in chosen)):
            d = txd.estimate_channel_delay(DATA, ch, zpm=zpm, mpm=mpm, tx_ant=tx_ant, rx_ant=rx_ant, min_tx_pwr=400e3)[0]
            delay[ch] = txd.DEFAULT_DELAY_US if d is None else d
        print("tx delays (us):", delay, flush=True)
        jobs = [(ch, per, arm) for ch, per, _, _ in chosen for arm in ARMS if load(args.out, ch, per, arm) is None]
        print("%d runs" % len(jobs), flush=True)
        running = []
        for ch, per, arm in jobs:
            while len(running) >= 64:
                running = [r for r in running if r.poll() is None]
                if len(running) >= 64:
                    running[0].wait()
            running.append(subprocess.Popen(
                ["systemd-run", "--user", "--scope", "--quiet", "-p", "MemoryMax=8G", sys.executable, __file__,
                 "--out", args.out, "--one", ch, str(per), arm, str(delay[ch])],
                stdout=open("%s/log_%s_p%d_%s.txt" % (args.out, ch, per, arm), "w"), stderr=subprocess.STDOUT))
        for r in running:
            r.wait()
    report(args.out, chosen)


if __name__ == "__main__":
    main()
