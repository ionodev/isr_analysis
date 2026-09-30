#!/usr/bin/env python3
"""
Regression benchmark for the review process (REVIEW_PROCESS.md).

Runs the pipeline's default stages on the fixed benchmark of
review/benchmark.json, once with the code of a base commit (default: main)
and once with the code of the commit under review, and compares every
output:

    lpi        outlier_lpi.lpi_files on single periods (the pipeline's defaults,
               outlier rejection on)
    fit_lpi    fit_lpi.fit_lpifiles on all LPI files of a channel, including a
               300 s stretch per channel
    rd300      avg_range_doppler_spec on single mode-300 periods

The outputs of a commit are cached under ~/isr_project/regression/runs/<sha>,
so the base usually runs only once. Each job runs in its own memory-limited
scope. The report goes to stdout and to ~/isr_project/regression/reports/.

    python3 review/regression.py <branch or commit> [--base main] [--jobs 16]

Exit status: 0 if every output is bit-identical to the base's, 1 if any
output differs, 2 if something failed or is missing. Gate A needs 0. Gate B
uses the report to quantify the change.
"""

import argparse
import datetime
import glob
import json
import os
import subprocess
import sys
import time

import numpy as n

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
ROOT = os.path.expanduser("~/isr_project/regression")
BENCH = json.load(open(os.path.join(HERE, "benchmark.json")))
DATA = BENCH["data"]
LPI = dict(rg=60, avg_dur=10, min_tx_frac=0.5, pass_band=1e5, filter_len=20, maximum_range_delay=7000,
           save_acf_images=False, lag_avg=1, reanalyze=True)
CHUNK = 8          # LPI periods per process: the metadata are read once per process


def git(*a, cwd=REPO):
    return subprocess.run(["git", *a], cwd=cwd, capture_output=True, text=True, check=True).stdout.strip()


# ---------------------------------------------------------------- workers

def worker(kind, code, out, channel, periods):
    """Run one job with the pipeline code of the worktree `code`."""
    sys.path.insert(0, code)
    os.chdir(code)
    if kind == "lpi":
        import outlier_lpi as olpi
        assert os.path.dirname(os.path.abspath(olpi.__file__)) == code
        olpi.lpi_files(dirname=DATA, channel=channel, output_base=out, periods=periods, **LPI)
    elif kind == "rd300":
        import avg_range_doppler_spec as ards
        ards.avg_range_doppler_spectra(dirname=DATA, channel=channel, mode=300, avg_dur=10, step=10,
                                       avg_type="outlier_removal", postfix="_outlier", min_tx_pulses=100,
                                       reanalyze=True, output_base=out, periods=periods)
    elif kind == "fit_lpi":
        import fit_lpi as flpi
        flpi.fit_lpifiles(dirn=DATA, channel=channel, postfix="_%d" % LPI["rg"], max_dt=300, plot=False,
                          first_lag=0, reanalyze=True, range_avg=n.array([1, 3, 5]), output_base=out,
                          radar_freq_hz=440.2e6, table_dir=BENCH["table_dir"])


# ---------------------------------------------------------------- running a commit

def jobs_for(stage):
    if stage == "lpi":
        out = []
        for ch in ("zenith-l", "misa-l"):
            per = sorted(int(p) for p in BENCH["lpi"][ch])
            a, b = BENCH["fit_lpi"][ch]
            per = sorted(set(per) | set(range(a, b + 1)))
            out += [(ch, per[i:i + CHUNK]) for i in range(0, len(per), CHUNK)]
        return [("lpi", ch, p) for ch, p in out]
    if stage == "rd300":
        return [("rd300", ch, [p]) for ch, ps in BENCH["range_doppler_300"].items() for p in ps]
    if stage == "fit_lpi":
        return [("fit_lpi", ch, []) for ch in BENCH["fit_lpi"]]


def run_jobs(jobs, code, out, logdir, n_jobs):
    """Run jobs in parallel, each in its own scope; returns the failed ones."""
    running, failed = [], []
    for i, (kind, ch, per) in enumerate(jobs):
        while len([r for r in running if r[0].poll() is None]) >= n_jobs:
            time.sleep(2)
        log = open(os.path.join(logdir, "%s_%s_%d.log" % (kind, ch, i)), "w")
        cmd = ["systemd-run", "--user", "--scope", "--quiet", "-p", "MemoryMax=10G", "-p", "MemorySwapMax=0",
               sys.executable, os.path.abspath(__file__), "--worker", kind, code, out, ch, ",".join(map(str, per))]
        env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                   HDF5_USE_FILE_LOCKING="FALSE", MPLBACKEND="Agg")
        running.append((subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env), (kind, ch, per)))
    for p, job in running:
        if p.wait() != 0:
            failed.append(job)
    return failed


def run_commit(sha, n_jobs):
    """Outputs of commit `sha` on the benchmark, cached."""
    out = os.path.join(ROOT, "runs", sha)
    if os.path.exists(os.path.join(out, "DONE")):
        print("%s: cached outputs in %s" % (sha[:7], out), flush=True)
        return out, []
    code = os.path.join(ROOT, "worktrees", sha)
    if not os.path.isdir(code):
        os.makedirs(os.path.dirname(code), exist_ok=True)
        git("worktree", "add", "--detach", code, sha)
    logdir = os.path.join(out, "logs")
    os.makedirs(logdir, exist_ok=True)
    t0 = time.time()
    print("%s: running the benchmark with %d parallel jobs" % (sha[:7], n_jobs), flush=True)
    failed = run_jobs(jobs_for("lpi") + jobs_for("rd300"), code, out, logdir, n_jobs)
    failed += run_jobs(jobs_for("fit_lpi"), code, out, logdir, n_jobs)
    if not failed:
        json.dump(dict(sha=sha, seconds=time.time() - t0, when=datetime.datetime.now(datetime.UTC).isoformat()),
                  open(os.path.join(out, "DONE"), "w"))
    print("%s: done in %.0f s, %d failed jobs" % (sha[:7], time.time() - t0, len(failed)), flush=True)
    return out, failed


# ---------------------------------------------------------------- comparing

def equal(a, b):
    if a.shape != b.shape or a.dtype.kind != b.dtype.kind:
        return False
    if a.dtype.kind in "fc":
        return bool(n.array_equal(a, b, equal_nan=True))
    return bool(n.array_equal(a, b))


def describe(name, a, b, h_base):
    """How a dataset changed, in the most meaningful unit available."""
    if a.shape != b.shape:
        return "shape %s -> %s" % (a.shape, b.shape)
    if a.dtype.kind not in "fciu":
        return "changed"
    a = a.astype(complex if a.dtype.kind == "c" else float)
    b = b.astype(complex if b.dtype.kind == "c" else float)
    fa, fb = n.isfinite(a), n.isfinite(b)
    txt = []
    if (fa & ~fb).any():
        txt.append("%d values lost (finite -> NaN)" % (fa & ~fb).sum())
    if (~fa & fb).any():
        txt.append("%d values recovered (NaN -> finite)" % (~fa & fb).sum())
    both = fa & fb
    var = {"acfs_e": "acfs_var", "acfs_g": "acfs_var"}.get(name)
    if var and var in h_base and both.any():
        s = n.sqrt(n.abs(h_base[var][()]))
        z = n.abs(b - a)[both] / s[both]
        z = z[n.isfinite(z)]
        if len(z):
            txt.append("|change| / sigma: median %.3g, max %.3g" % (n.median(z), n.max(z)))
    elif both.any():
        d = n.abs(b - a)[both]
        r = d / n.maximum(n.abs(a)[both], 1e-30)
        txt.append("relative change: median %.3g, max %.3g" % (n.median(r), n.max(r)))
    return "; ".join(txt) or "changed"


def compare(base_out, test_out):
    import h5py
    rel = lambda d: {os.path.relpath(f, d) for f in glob.glob(os.path.join(d, "**", "*.h5"), recursive=True)}
    fb, ft = rel(base_out), rel(test_out)
    rows, n_ident = [], 0
    for f in sorted(fb - ft):
        rows.append((f, "", "only in base"))
    for f in sorted(ft - fb):
        rows.append((f, "", "only in the commit under review"))
    for f in sorted(fb & ft):
        with h5py.File(os.path.join(base_out, f), "r") as hb, h5py.File(os.path.join(test_out, f), "r") as ht:
            kb, kt = set(hb.keys()), set(ht.keys())
            changed = False
            for k in sorted(kb ^ kt):
                rows.append((f, k, "dataset only in %s" % ("base" if k in kb else "the commit under review")))
                changed = True
            for k in sorted(kb & kt):
                a, b = n.asarray(hb[k][()]), n.asarray(ht[k][()])
                if not equal(a, b):
                    rows.append((f, k, describe(k, a, b, hb)))
                    changed = True
            n_ident += not changed
    return rows, n_ident, len(fb | ft)


# ---------------------------------------------------------------- main

def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        kind, code, out, ch, per = sys.argv[2:7]
        worker(kind, code, out, ch, [int(p) for p in per.split(",") if p])
        return
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("target", help="branch or commit under review")
    p.add_argument("--base", default="main")
    p.add_argument("--jobs", type=int, default=16)
    args = p.parse_args()
    base, test = git("rev-parse", args.base), git("rev-parse", args.target)
    print("base %s (%s), under review %s (%s)" % (base[:7], args.base, test[:7], args.target), flush=True)
    base_out, fail_b = run_commit(base, args.jobs)
    test_out, fail_t = run_commit(test, args.jobs)
    rows, n_ident, n_files = compare(base_out, test_out)
    lines = ["# Regression benchmark: %s (%s) against %s (%s)" % (args.target, test[:7], args.base, base[:7]), "",
             "Run %s UTC. %d output files; %d bit-identical." % (
                 datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%d %H:%M"), n_files, n_ident), ""]
    for label, fails in (("base", fail_b), ("under review", fail_t)):
        for kind, ch, per in fails:
            lines.append("- FAILED job (%s): %s %s periods %s; see the logs under %s" % (
                label, kind, ch, per, os.path.join(base_out if label == "base" else test_out, "logs")))
    if rows:
        lines += ["| file | dataset | change |", "|---|---|---|"] + ["| %s | %s | %s |" % r for r in rows]
    else:
        lines.append("Every output is bit-identical.")
    txt = "\n".join(lines) + "\n"
    print(txt)
    os.makedirs(os.path.join(ROOT, "reports"), exist_ok=True)
    name = "%s-%s_vs_%s-%s.md" % (args.target.replace("/", "_"), test[:7], args.base.replace("/", "_"), base[:7])
    open(os.path.join(ROOT, "reports", name), "w").write(txt)
    print("report: %s" % os.path.join(ROOT, "reports", name))
    sys.exit(2 if (fail_b or fail_t) else (1 if rows else 0))


if __name__ == "__main__":
    main()
