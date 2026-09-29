#!/usr/bin/env python3
"""
Does the line interference of memo 20 reach the range-Doppler spectra of
MISA's low elevation mode 800?

Mode 800 is processed by avg_range_doppler_spec.py, not the lag profile
inversion: a 60 kHz low-pass filter, then per range gate the spectrum of the
echo times the transmitted pulse, keeping |f| < 50 kHz, averaged over the 10 s
period with or without an outlier test.  The worst periods of the class N
interference in mode 800, ranked by its power within 55 kHz of the radar
frequency, and moderate ones, are each processed

    outlier       avg_type "outlier_removal", as the pipeline runs
    mean          avg_type "mean", no outlier test
    outlier_masked, mean_masked
                  the same without the class N pulses: the references

The difference from the reference is given in units of the reference's
standard deviation of the period mean, sigma = sqrt(RDS_LP_var / alpha) (the
file stores the pulse variance / alpha / n_pulses), and as the change of the
signal-to-noise ratio that fit_lp.py fits, (RDS - noise) / noise, noise from
3600-4500 km.  Every run is its own process, all at once.

    python3 validate_rfi_rd800.py [--out plots/rfi_rd800_test]
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
RESULTS = os.path.expanduser("~/isr_project/documents/analysis/results")
CACHE = os.path.join(os.path.expanduser("~/.cache/isr_analysis"), os.path.basename(DATA))
T0 = 1712484000
N_PER = 17640
ARMS = {"outlier": ("outlier_removal", False), "mean": ("mean", False),
        "outlier_masked": ("outlier_removal", True), "mean_masked": ("mean", True)}
REF = {"outlier": "outlier_masked", "mean": "mean_masked"}
NOISE_KM = (3600, 4500)


def class_n_mode800():
    """Keys of the class N misa-l pulses sent in mode 800 on MISA, and their power within 55 kHz."""
    s = n.load(os.path.join(RESULTS, "misa_rfi_psd_summary.npz"))
    z = n.load(os.path.join(RESULTS, "misa_rfi_classes.npz"))
    P = n.load(os.path.join(CACHE, "misa_rfi_psd.npz"))["psd"][:, 0].astype(n.float32)
    f = s["freq_hz"]
    m800 = (s["sweepid"] == 800) & (s["tx_ant"] > 0.99)
    band = n.abs(f) < 55e3
    ref = 10**(n.nanmedian(P[s["control"] & m800], axis=0) / 10)
    excess = (10**(P[:, band] / 10)).sum(1) / ref[band].sum() - 1
    N = (z["misa_class"] == "N") & s["flag_misa"] & m800
    return s["keys"][N], excess[N]


def choose(n_worst=12, n_moderate=8):
    keys, ex = class_n_mode800()
    per = ((keys / 1e6 - T0) // 10).astype(int)
    tot = n.bincount(per, weights=n.nan_to_num(ex), minlength=N_PER) / 289.0
    worst = [int(p) for p in n.argsort(tot)[::-1][:n_worst]]
    mid = n.where((tot > 0.01) & (tot < 0.5))[0]
    mid = [int(p) for p in mid[n.linspace(0, len(mid) - 1, n_moderate).astype(int)]]
    return [(p, "worst", float(tot[p])) for p in worst] + [(p, "moderate", float(tot[p])) for p in mid]


def run_one(out, period, arm, tx_delay_us):
    import avg_range_doppler_spec as ards
    import millstone_radar_state as mrs
    avg_type, masked = ARMS[arm]
    if masked:
        keys, _ = class_n_mode800()
        t = keys / 1e6
        mask = set(int(k) for k in keys[(t >= T0 + 10 * period - 1) & (t < T0 + 10 * period + 12)])
        print("masking %d pulses" % len(mask))
        original = mrs.get_antenna_select

        def masked_select(path):
            tx_ant, rx_ant = original(path)
            return (lambda k: -1.0 if int(k) in mask else tx_ant(k)), rx_ant   # as if sent on zenith: skipped
        ards.mrs.get_antenna_select = masked_select
    ards.avg_range_doppler_spectra(dirname=DATA, channel="misa-l", mode=800, avg_dur=10, step=10,
                                   avg_type=avg_type, postfix="_" + arm, reanalyze=True,
                                   output_base="%s/p%d" % (out, period), tx_delay_us=tx_delay_us,
                                   periods=[period])


def load(out, period, arm):
    f = glob.glob("%s/p%d/range_doppler_800_%s/misa-l/il_*.h5" % (out, period, arm))
    if not f:
        return None
    with h5py.File(f[0], "r") as h:
        return {k: h[k][()] for k in ("RDS_LP", "RDS_LP_var", "alpha", "T_sys", "n_pulses", "rgs_km", "dop_hz")}


def snr(d):
    r = (d["rgs_km"] > NOISE_KM[0]) & (d["rgs_km"] < NOISE_KM[1])
    noise = n.median(d["RDS_LP"][r])
    return (d["RDS_LP"] - noise) / noise


def report(out, chosen):
    lines = ["Effect of the class N interference on mode 800 range-Doppler spectra (misa-l).",
             "z = (arm - reference) / sigma_reference over all ranges and Doppler bins;",
             "dSNR = change of (RDS - noise)/noise at ranges below 1500 km", ""]
    res = []
    for p, kind, extra in chosen:
        txt = "p%-5d %-8s in-band excess %7.2f" % (p, kind, extra)
        row = dict(period=p, kind=kind, inband_excess=extra)
        for arm, rarm in REF.items():
            a, m = load(out, p, arm), load(out, p, rarm)
            if a is None or m is None:
                txt += " | %s: missing" % arm
                continue
            sig = n.sqrt(m["RDS_LP_var"] / m["alpha"])
            zz = (a["RDS_LP"] - m["RDS_LP"]) / sig
            ok = n.isfinite(zz)
            low = m["rgs_km"] < 1500
            ds = (snr(a) - snr(m))[low]
            row[arm] = dict(zmax=float(n.max(n.abs(zz[ok]))), zmed=float(n.median(n.abs(zz[ok]))),
                            zmean=float(n.mean(zz[ok])), dsnr_max=float(n.nanmax(n.abs(ds))),
                            dsnr_med=float(n.nanmedian(n.abs(ds))), snr_max_ref=float(n.nanmax(snr(m)[low])),
                            tsys_ratio=float(a["T_sys"] / m["T_sys"]), n_pulses=int(a["n_pulses"]),
                            n_pulses_ref=int(m["n_pulses"]))
            r = row[arm]
            txt += " | %s: |z|max %6.1f med %.2f mean %+.2f  |dSNR|max %.3f  Tsys x%.3f  n %d/%d" % (
                arm, r["zmax"], r["zmed"], r["zmean"], r["dsnr_max"], r["tsys_ratio"], r["n_pulses"], r["n_pulses_ref"])
        lines.append(txt)
        res.append(row)
    txt = "\n".join(lines)
    print(txt)
    open("%s/report.txt" % out, "w").write(txt + "\n")
    json.dump(res, open("%s/results.json" % out, "w"), indent=1)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default="plots/rfi_rd800_test")
    p.add_argument("--report", action="store_true")
    p.add_argument("--one", nargs=3, help=argparse.SUPPRESS)
    args = p.parse_args()
    if args.one:
        run_one(args.out, int(args.one[0]), args.one[1], float(args.one[2]))
        return
    os.makedirs(args.out, exist_ok=True)
    chosen = choose()
    json.dump(chosen, open("%s/chosen.json" % args.out, "w"), indent=1)
    if not args.report:
        # one transmit delay for all runs, measured as the pipeline would
        import millstone_radar_state as mrs
        import tx_delay as txd
        zpm, mpm = mrs.get_tx_power_model("%s/metadata/powermeter" % DATA)
        tx_ant, rx_ant = mrs.get_antenna_select("%s/metadata/antenna_control_metadata" % DATA)
        d = txd.estimate_channel_delay(DATA, "misa-l", zpm=zpm, mpm=mpm, tx_ant=tx_ant, rx_ant=rx_ant,
                                       min_tx_pwr=400e3)[0]
        d = txd.DEFAULT_DELAY_US if d is None else d
        print("tx delay %.3f us" % d, flush=True)
        jobs = [(per, arm) for per, _, _ in chosen for arm in ARMS if load(args.out, per, arm) is None]
        print("%d runs" % len(jobs), flush=True)
        procs = [subprocess.Popen([sys.executable, __file__, "--out", args.out, "--one", str(per), arm, str(d)],
                                  stdout=open("%s/log_p%d_%s.txt" % (args.out, per, arm), "w"), stderr=subprocess.STDOUT)
                 for per, arm in jobs]
        for pr in procs:
            pr.wait()
    report(args.out, chosen)


if __name__ == "__main__":
    main()
