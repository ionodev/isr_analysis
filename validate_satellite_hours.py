#!/usr/bin/env python3
"""
The satellite columns of the lag profile inversion over whole hours, on
zenith-l and misa-l.

Memos 15 and 16 tested the columns on injected echoes and on two real passes.
This runs every 10-second period of one hour per channel three ways:

    nohandling  no satellite handling, no outlier rejection
    current     the pipeline's default, the outlier rejection
    fit         catalogued echoes as unknowns (satellite columns), no rejection

and reports, over the hour:
  - away from the echoes (gates more than 1200 us from every catalogued echo
    of the period), fit against nohandling in units of nohandling's standard
    deviation: the columns must not disturb the plasma elsewhere;
  - at the echoes (within 600 us), fit against current, and the ratio of their
    error bars, and how many lagged products each leaves unmeasured;
  - the satellite amplitudes: finite fraction, and how flat |amplitude| is
    with lag (it is the echo energy and should not change);
  - failures and run time.

The catalogue's delays are converted with the channel delay (tx_delay.py),
measured at the start of every block of consecutive periods, since it can
change within an hour (Memo 26).

    python3 validate_satellite_hours.py [--out plots/satellite_hours] [--hour 2024-04-08T16]

Only periods in which the channel's own antenna sends at least half of the
pulses in the inversion's modes (coded, 300) are run.
"""

import os
for _v in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ[_v] = "1"
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import argparse
import datetime
import glob
import json
import subprocess
import sys
import time

import numpy as n
import h5py

DATA = "/mnt/data/juha/millstone_hill/isr/eclipse2024/usrp-rx0-r_20240407T100000_20240409T110000"
SAT_MD = "%s/metadata/satellite_detections" % DATA
CHANNELS = ["zenith-l", "misa-l"]
LPI = dict(avg_dur=10, rg=60, min_tx_frac=0.5, pass_band=1e5, filter_len=20, maximum_range_delay=7000,
           save_acf_images=False, lag_avg=1, reanalyze=True)
ARMS = ["nohandling", "current", "fit"]
AWAY_US, FOOT_US = 1200, 600


def run_block(out, ch, arm, bi, periods, delay_us):
    """one arm over one block of periods: lpi_files reads its metadata once per call"""
    import outlier_lpi as olpi
    import satellite_columns as satcol
    t = time.time()
    kw = dict(LPI)
    if arm == "fit":
        kw["satellite_detections"] = satcol.CatalogueDetections(SAT_MD, ch, delay_us)
    olpi.lpi_files(dirname=DATA, channel=ch, output_base="%s/%s/%s" % (out, ch, arm), periods=periods,
                   outlier_rejection=(arm == "current"), **kw)
    json.dump(dict(seconds=time.time() - t, periods=periods), open("%s/%s/time_%s_%d.json" % (out, ch, arm, bi), "w"))


def load(out, ch, arm, t0):
    f = "%s/%s/%s/lpi_%d/%s/lpi-%d.h5" % (out, ch, arm, LPI["rg"], ch, t0)
    if not os.path.exists(f):
        return None
    with h5py.File(f, "r") as h:
        return {k: h[k][()] for k in h.keys()}


def report(out, plan):
    import satellite_columns as satcol
    lines = []
    res = {}
    for ch in CHANNELS:
        p = plan[ch]
        dets = [satcol.CatalogueDetections(SAT_MD, ch, b["delay_us"]) for b in p["blocks"]]
        block_of = {ai: bi for bi, b in enumerate(p["blocks"]) for ai in b["periods"]}
        files = sorted(glob.glob("%s/%s/fit/lpi_%d/%s/lpi-*.h5" % (out, ch, LPI["rg"], ch)))
        z_away, z_foot, eb_foot, nan_cur, nan_fit, n_foot, flat, fin, n_cols, t_fit, t_noh, n_bright = ([] for _ in range(12))
        n_per, n_sat_per, fails = 0, 0, []
        for f in files:
            t0 = int(os.path.basename(f)[4:-3])
            R = {a: load(out, ch, a, t0) for a in ARMS}
            if any(R[a] is None for a in ARMS):
                fails.append(t0)
                continue
            n_per += 1
            i0 = int(R["fit"]["i0"] * 1e6)
            det = dets[block_of[int(round((i0 - p["b0"]) / (LPI["avg_dur"] * 1e6)))]]
            sat = det.for_period(i0, i0 + int(LPI["avg_dur"] * 1e6) + 40000)
            delays = n.array([d for x in sat.values() for d, _ in x])
            gate = n.arange(R["fit"]["acfs_e"].shape[0]) * LPI["rg"]
            if len(delays):
                dist = n.min(n.abs(gate[:, None] - delays[None, :]), axis=1)
                n_sat_per += 1
            else:
                dist = n.full(len(gate), n.inf)
            away, foot = dist > AWAY_US, dist < FOOT_US
            a, b = R["fit"], R["nohandling"]
            z = (a["acfs_e"] - b["acfs_e"]) / n.sqrt(b["acfs_var"])
            z_away.append(z[away][n.isfinite(z[away])])
            if foot.any():
                c = R["current"]
                z = (a["acfs_e"] - c["acfs_e"]) / n.sqrt(c["acfs_var"])
                z_foot.append(z[foot][n.isfinite(z[foot])])
                eb = n.sqrt(a["acfs_var"][foot] / c["acfs_var"][foot])
                eb_foot.append(eb[n.isfinite(eb)])
                nan_cur.append(int(n.sum(~n.isfinite(c["acfs_e"][foot]))))
                nan_fit.append(int(n.sum(~n.isfinite(a["acfs_e"][foot]))))
                n_foot.append(int(foot.sum() * a["acfs_e"].shape[1]))
            if "sat_amp_e" in a and len(a["sat_amp_e"]):
                x = a["sat_amp_e"]
                n_cols.append(len(x))
                ok = n.all(n.isfinite(x[:, :36]), axis=1)
                fin.append(ok.mean())
                # flatness only where the echo stands 5 sigma above its own noise
                bright = ok & (n.abs(x[:, 0])**2 > 25 * a["sat_amp_var"][:, 0])
                n_bright.append(int(bright.sum()))
                if bright.any():
                    flat.append(n.abs(x[bright, :36]) / n.abs(x[bright, :1]))
        for bi in range(len(p["blocks"])):
            tt = {a: "%s/%s/time_%s_%d.json" % (out, ch, a, bi) for a in ("fit", "nohandling")}
            if all(os.path.exists(f) for f in tt.values()):
                t_fit.append(json.load(open(tt["fit"]))["seconds"])
                t_noh.append(json.load(open(tt["nohandling"]))["seconds"])
        dl = n.array([b["delay_us"] for b in p["blocks"]])
        za, zf, ebf = (n.concatenate(v) if v else n.zeros(0) for v in (z_away, z_foot, eb_foot))
        fl = n.median(n.concatenate(flat), axis=0) if flat else n.zeros(0)
        r = dict(periods=n_per, periods_with_echoes=n_sat_per, missing=len(fails), delay_us_min_max=[dl.min(), dl.max()],
                 blocks=len(dl),
                 away_rms=float(n.sqrt(n.mean(n.abs(za)**2))) if len(za) else None,
                 away_max=float(n.max(n.abs(za))) if len(za) else None,
                 away_frac_gt3=float(n.mean(n.abs(za) > 3)) if len(za) else None,
                 foot_rms=float(n.sqrt(n.mean(n.abs(zf)**2))) if len(zf) else None,
                 foot_errbar_median=float(n.median(ebf)) if len(ebf) else None,
                 foot_errbar_p16_p84=[float(x) for x in n.percentile(ebf, [16, 84])] if len(ebf) else None,
                 foot_unmeasured_current=float(n.sum(nan_cur) / max(1, n.sum(n_foot))),
                 foot_unmeasured_fit=float(n.sum(nan_fit) / max(1, n.sum(n_foot))),
                 columns_median=float(n.median(n_cols)) if n_cols else 0, columns_max=int(max(n_cols)) if n_cols else 0,
                 amp_finite=float(n.mean(fin)) if fin else None, bright_columns=int(sum(n_bright)),
                 delay_us_blocks=[float(x) for x in dl],
                 amp_flatness_min_max=[float(fl.min()), float(fl.max())] if len(fl) else None,
                 time_fit_over_nohandling=float(n.median(n.array(t_fit) / n.array(t_noh))) if t_fit else None,
                 time_fit_block_median_s=float(n.median(t_fit)) if t_fit else None)
        res[ch] = r
        lines += ["%s, hour from %s: %d periods (%d missing an arm), %d with catalogued echoes; channel delay %.2f-%.2f us "
                  "over %d blocks" % (ch, p["hour"], n_per, len(fails), n_sat_per, dl.min(), dl.max(), len(dl)),
                  "  away from echoes, fit - nohandling: rms %s sigma, max %s, fraction > 3 sigma %s"
                  % tuple("%.3g" % x if x is not None else "-" for x in (r["away_rms"], r["away_max"], r["away_frac_gt3"])),
                  "  at echoes, fit - current: rms %s sigma; error bar fit/current median %s (16-84%% %s);"
                  " unmeasured lagged products: current %.3f, fit %.3f"
                  % ("%.3g" % r["foot_rms"] if r["foot_rms"] is not None else "-",
                     "%.3f" % r["foot_errbar_median"] if r["foot_errbar_median"] is not None else "-",
                     "%.3f-%.3f" % tuple(r["foot_errbar_p16_p84"]) if r["foot_errbar_p16_p84"] else "-",
                     r["foot_unmeasured_current"], r["foot_unmeasured_fit"]),
                  "  satellite columns per period: median %s, max %s; amplitudes finite %s; %d columns 5 sigma bright, "
                  "median |amp(lag)|/|amp(first)| over them %s"
                  % (r["columns_median"], r["columns_max"], "%.3f" % r["amp_finite"] if r["amp_finite"] is not None else "-",
                     r["bright_columns"],
                     "%.3f-%.3f" % tuple(r["amp_flatness_min_max"]) if r["amp_flatness_min_max"] else "-"),
                  "  run time: fit %s s per block, %s times nohandling" %
                  ("%.0f" % r["time_fit_block_median_s"] if r["time_fit_block_median_s"] else "-",
                   "%.2f" % r["time_fit_over_nohandling"] if r["time_fit_over_nohandling"] else "-"),
                  "  channel delay per block (us): " + " ".join("%.2f" % x for x in dl), ""]
    txt = "\n".join(lines)
    print(txt)
    open("%s/report.txt" % out, "w").write(txt)
    json.dump(res, open("%s/results.json" % out, "w"), indent=1)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default="plots/satellite_hours")
    p.add_argument("--hour", default="2024-04-08T16")
    p.add_argument("--jobs", type=int, default=48)
    p.add_argument("--max-periods", type=int, default=360)
    p.add_argument("--report", action="store_true")
    p.add_argument("--block", type=int, default=8, help="periods per process")
    p.add_argument("--one", nargs=3, help=argparse.SUPPRESS)
    args = p.parse_args()
    pf = "%s/plan.json" % args.out
    if args.one:
        ch, arm, bi = args.one[0], args.one[1], int(args.one[2])
        b = json.load(open(pf))[ch]["blocks"][bi]
        run_block(args.out, ch, arm, bi, b["periods"], b["delay_us"])
        return
    os.makedirs(args.out, exist_ok=True)
    if os.path.exists(pf):
        plan = json.load(open(pf))
    else:
        plan = make_plan(args)
        json.dump(plan, open(pf, "w"), indent=1)
    if not args.report:
        jobs = [(ch, arm, bi) for ch in CHANNELS for bi in range(len(plan[ch]["blocks"])) for arm in ARMS
                if not os.path.exists("%s/%s/time_%s_%d.json" % (args.out, ch, arm, bi))]
        print("%d blocks to run" % len(jobs), flush=True)
        running = []
        for ch, arm, bi in jobs:
            while len(running) >= args.jobs:
                running = [r for r in running if r.poll() is None]
                if len(running) >= args.jobs:
                    time.sleep(1)
            os.makedirs("%s/%s" % (args.out, ch), exist_ok=True)
            running.append(subprocess.Popen(
                ["systemd-run", "--user", "--scope", "--quiet", "-p", "MemoryMax=6G", sys.executable, __file__,
                 "--out", args.out, "--one", ch, arm, str(bi)],
                stdout=open("%s/%s/log_%s_%d.txt" % (args.out, ch, arm, bi), "w"), stderr=subprocess.STDOUT))
        for r in running:
            r.wait()
    report(args.out, plan)


def make_plan(args):
    """The periods of the hour in which the channel's own antenna sends at least
    half of the pulses in the inversion's modes, in blocks of consecutive ones,
    each with the channel delay measured at its start."""
    from digital_rf import DigitalMetadataReader
    import millstone_radar_state as mrs
    import tx_delay
    from raw_reader import pulse_index
    b0 = DigitalMetadataReader("%s/metadata/id_metadata" % DATA).get_bounds()[0]
    t_hour = datetime.datetime.strptime(args.hour, "%Y-%m-%dT%H").replace(tzinfo=datetime.timezone.utc).timestamp()
    a0 = int(round((t_hour - b0 / 1e6) / LPI["avg_dur"]))
    zpm, mpm = mrs.get_tx_power_model("%s/metadata/powermeter" % DATA)
    tx_ant, rx_ant = mrs.get_antenna_select("%s/metadata/antenna_control_metadata" % DATA)
    ix = pulse_index(DATA)
    step = int(LPI["avg_dur"] * 1e6)
    plan = {}
    for ch in CHANNELS:
        own = -1 if ch.startswith("zenith") else 1
        good = []
        for ai in range(a0, a0 + args.max_periods):
            i0 = b0 + ai * step
            j0, j1 = n.searchsorted(ix["key"], [i0, i0 + step])
            sid = ix["sweepid"][j0:j1]
            inv = ((sid >= 1) & (sid <= 32)) | (sid == 300)
            mine = n.array([tx_ant(int(k)) * own > 0.99 for k in ix["key"][j0:j1]])
            if len(sid) and n.mean(mine & inv) >= 0.5:
                good.append(ai)
        runs = n.split(good, n.where(n.diff(good) > 1)[0] + 1) if good else []
        blocks = []
        for r in runs:
            for c in range(0, len(r), args.block):
                per = [int(x) for x in r[c:c + args.block]]
                d = tx_delay.estimate_channel_delay(DATA, ch, t0_unix=(b0 + per[0] * step) / 1e6, verbose=False,
                                                    zpm=zpm, mpm=mpm, tx_ant=tx_ant, rx_ant=rx_ant)[0]
                if d is None:
                    print("%s: no channel delay at period %d, block skipped" % (ch, per[0]), flush=True)
                    continue
                blocks.append(dict(periods=per, delay_us=float(d)))
        plan[ch] = dict(hour=args.hour, b0=int(b0), blocks=blocks)
        print("%s: %d periods in %d blocks; channel delay %s us" % (ch, sum(len(b["periods"]) for b in blocks),
              len(blocks), " ".join("%.2f" % b["delay_us"] for b in blocks)), flush=True)
    return plan


if __name__ == "__main__":
    main()
