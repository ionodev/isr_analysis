#!/usr/bin/env python3
"""
Test of raw_reader.RawReader in the lag profile inversion (fast_read=True
against DigitalRFReader).

    1. speed: two interleaved sets of integration periods, spread over the whole
       recording and not in the page cache, are inverted, one set with each
       reader, every period its own process and each set all at once, as a
       production run would be; the wall time of each set is compared;
    2. identity: the fast set is inverted again with DigitalRFReader, and every
       dataset of every output file must be bit identical.

Only periods in which zenith transmits and receives are used, and only ones
whose data are not in the page cache.  That is probed by timing one 4 kB read
from the middle of each of a period's files (fincore cannot be used: for files
another user owns the kernel reports every page as cached).  Each inversion
also reports the time it spent reading.

    python3 validate_raw_reader.py [--out plots/raw_reader_test] [--n 48]
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
CHANNEL = "zenith-l"
T0 = 1712484000                    # first second of the recording; period i starts at T0 + 10 i
N_PERIODS = 17640
LPI = dict(avg_dur=10, rg=60, min_tx_frac=0.2, filter_len=100, maximum_range_delay=7200,
           save_acf_images=False, lag_avg=1, reanalyze=True, pass_band=0.018e6)


def run_one(out, period, fast):
    import outlier_lpi as olpi
    import raw_reader
    from digital_rf import DigitalRFReader
    spent = [0.0]
    for cls in (raw_reader.RawReader, DigitalRFReader):
        f = cls.read_vector_1d

        def timed(self, *a, _f=f, **k):
            t = time.time()
            try:
                return _f(self, *a, **k)
            finally:
                spent[0] += time.time() - t
        cls.read_vector_1d = timed
    t = time.time()
    olpi.lpi_files(dirname=DATA, channel=CHANNEL, output_base="%s/%s/p%d" % (out, "fast" if fast else "drf", period),
                   periods=[period], fast_read=fast, **LPI)
    print("elapsed %.1f s reading %.1f s" % (time.time() - t, spent[0]))


def files_of(period):
    out = []
    for s in range(T0 + 10 * period, T0 + 10 * period + 10):
        sub = datetime.datetime.fromtimestamp(s - s % 3600, datetime.timezone.utc).strftime("%Y-%m-%dT%H-%M-%S")
        out.append("%s/rf_data/%s/%s/rf@%d.000.h5" % (DATA, CHANNEL, sub, s))
    return [f for f in out if os.path.exists(f)]


def probe_cached(period, limit_s=0.5e-3):
    """Fraction of the period's files whose middle page comes back faster than a disk could deliver it."""
    fs = files_of(period)
    hits = 0
    for f in fs:
        fd = os.open(f, os.O_RDONLY)
        try:
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_RANDOM)     # no read-ahead around the probe
            t = time.perf_counter()
            os.pread(fd, 4096, 2_000_000)
            hits += time.perf_counter() - t < limit_s
        finally:
            os.close(fd)
    return hits / max(len(fs), 1)


def zenith_periods(candidates):
    import millstone_radar_state as mrs
    tx_ant, rx_ant = mrs.get_antenna_select("%s/metadata/antenna_control_metadata" % DATA)
    zpm, _ = mrs.get_tx_power_model("%s/metadata/powermeter" % DATA)
    ok = []
    for q in candidates:
        ks = [(T0 + 10 * q + d) * 10**6 for d in (1, 5, 9)]
        if all(tx_ant(k) <= -0.99 and rx_ant(k) <= -0.99 and zpm(k / 1e6) >= 400e3 for k in ks):
            ok.append(int(q))
    return ok


def batch(out, periods, fast, log_tag):
    os.makedirs(out, exist_ok=True)
    t = time.time()
    procs = [subprocess.Popen([sys.executable, __file__, "--out", out, "--one", str(p), str(int(fast))],
                              stdout=open("%s/log_%s_p%d.txt" % (out, log_tag, p), "w"), stderr=subprocess.STDOUT)
             for p in periods]
    bad = [p for p, pr in zip(periods, procs) if pr.wait() != 0]
    return time.time() - t, bad


def outputs(out, which, period):
    return sorted(glob.glob("%s/%s/p%d/lpi_%d/%s/*.h5" % (out, which, period, LPI["rg"], CHANNEL)))


def identical(fa, fb):
    with h5py.File(fa, "r") as a, h5py.File(fb, "r") as b:
        if set(a.keys()) != set(b.keys()):
            return False
        for k in a.keys():
            x, y = a[k][()], b[k][()]
            if isinstance(x, n.ndarray) and x.dtype.kind in "fc":
                if not n.array_equal(x, y, equal_nan=True):
                    return False
            elif not n.array_equal(x, y):
                return False
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default="plots/raw_reader_test")
    p.add_argument("--n", type=int, default=48, help="periods per set")
    p.add_argument("--one", nargs=2, help=argparse.SUPPRESS)
    args = p.parse_args()
    if args.one:
        run_one(args.out, int(args.one[0]), bool(int(args.one[1])))
        return
    cand = zenith_periods(n.linspace(30, N_PERIODS - 30, 12 * args.n).astype(int))
    cache = {}
    cold = []
    for q in cand:
        cache[q] = probe_cached(q)
        if cache[q] == 0:
            cold.append(q)
    print("%d candidate periods with zenith transmitting, %d of them not in the page cache" % (len(cand), len(cold)), flush=True)
    pick = [cold[i] for i in n.linspace(0, len(cold) - 1, min(2 * args.n, len(cold))).astype(int)]
    set_drf, set_fast = pick[0::2], pick[1::2]
    res = dict(periods_drf=set_drf, periods_fast=set_fast, cached_probe=cache)
    res["wall_drf"], res["failed_drf"] = batch(args.out, set_drf, False, "drf")
    print("DigitalRFReader, %d periods at once: %.0f s" % (len(set_drf), res["wall_drf"]), flush=True)
    res["wall_fast"], res["failed_fast"] = batch(args.out, set_fast, True, "fast")
    print("RawReader,       %d periods at once: %.0f s" % (len(set_fast), res["wall_fast"]), flush=True)
    # identity: the fast periods again with digital_rf (now from the page cache)
    res["wall_drf_warm"], _ = batch(args.out, set_fast, False, "drf_warm")
    print("DigitalRFReader on the same periods, from the page cache: %.0f s" % res["wall_drf_warm"], flush=True)
    comp = {}
    for q in set_fast:
        fa, fb = outputs(args.out, "fast", q), outputs(args.out, "drf", q)
        comp[q] = dict(files=len(fa), same_names=[os.path.basename(f) for f in fa] == [os.path.basename(f) for f in fb],
                       identical=len(fa) == len(fb) and all(identical(a, b) for a, b in zip(fa, fb)))
    res["comparison"] = comp
    n_out = sum(c["files"] for c in comp.values())
    print("outputs compared: %d files from %d periods; all bit identical: %s"
          % (n_out, len(comp), all(c["identical"] and c["same_names"] for c in comp.values())), flush=True)
    for tag in ("drf", "fast", "drf_warm"):
        el, rd = [], []
        for q in (set_drf if tag == "drf" else set_fast):
            for line in open("%s/log_%s_p%d.txt" % (args.out, tag, q)):
                if line.startswith("elapsed"):
                    w = line.split()
                    el.append(float(w[1]))
                    rd.append(float(w[4]))
        res["time_%s" % tag] = dict(elapsed_mean=float(n.mean(el)), reading_mean=float(n.mean(rd)))
        print("%-8s per inversion: %.1f s, of which reading %.1f s" % (tag, n.mean(el), n.mean(rd)), flush=True)
    json.dump(res, open("%s/results.json" % args.out, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
