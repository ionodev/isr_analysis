#!/usr/bin/env python3
"""
Regression benchmark for the review process (REVIEW_PROCESS.md).

Runs pipeline stages on the fixed benchmark of review/benchmark.json, once
with the code of a base commit (default: main) and once with the code of the
commit under review, and compares every output byte for byte:

    lpi        outlier_lpi.lpi_files on single periods and on a 300 s stretch,
               with the arguments run_analysis.py passes for
               config/millstone_eclipse2024.json (outlier rejection on)
    fit_lpi    fit_lpi.fit_lpifiles on all LPI files of a channel
    rd300      avg_range_doppler_spec on single mode-300 periods
    fit_lp     fit_lp.fit_spectra on the rd300 output

Not covered: run_analysis.py itself (its config mapping), how the theory
tables are generated (the tables are read from table_dir as they are),
mode 800, options that are off by default, and data outside the benchmark.

A commit's outputs are cached under ~/isr_project/regression/runs/, keyed by
the commit and by a hash of everything else the outputs depend on: this
script, benchmark.json, the installed Python packages and the table files.
Each job runs in its own memory-limited scope. The report goes to stdout and
to ~/isr_project/regression/reports/.

    python3 review/regression.py <branch or commit> [--base main] [--jobs 16]

Run main's copy of this script (see REVIEW_PROCESS.md). The commit under
review must contain the base (merge or rebase main into the branch first).

Exit status: 0 every output is bit-identical to the base's; 1 some output
differs; 2 the run is not valid: a job failed, an expected output is missing
on either side, nothing was compared, the target does not contain the base,
or the tool itself failed. Gate A needs 0.
"""

import argparse
import datetime
import fcntl
import glob
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import traceback

import numpy as n

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
ROOT = os.path.expanduser("~/isr_project/regression")
BENCH = json.load(open(os.path.join(HERE, "benchmark.json")))
DATA = BENCH["data"]
LPI = dict(rg=60, avg_dur=10, min_tx_frac=0.5, pass_band=1e5, filter_len=20, maximum_range_delay=7000,
           save_acf_images=False, lag_avg=1, reanalyze=True)
CHUNK = 8                   # LPI periods per process: the metadata are read once per process
MIN_STRETCH_FILES = 10      # the fit stretch must give at least this many LPI files per channel


class Invalid(Exception):
    """The run cannot be judged: exit status 2."""


def git(*a, cwd=REPO):
    r = subprocess.run(["git", *a], cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        raise Invalid("git %s: %s" % (" ".join(a), r.stderr.strip()))
    return r.stdout.strip()


# ---------------------------------------------------------------- workers

def worker(kind, code, out, channel, periods):
    """Run one job with the pipeline code of the worktree `code`."""
    sys.path.insert(0, code)
    os.chdir(code)
    if kind == "lpi":
        import outlier_lpi as olpi
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
    elif kind == "fit_lp":
        import fit_lp as flp
        flp.fit_spectra(dirname=DATA, channel=channel, postfix="_300_outlier", avg_dur=600, ridx=[35, 230],
                        remove_space_objects=False, reanalyze=True, output_base=out, radar_freq_hz=440.2e6,
                        table_dir=BENCH["table_dir"], fit_bandwidth_hz=50e3, notch_bands_hz=None,
                        pulse_length_us=None)
    # every module of the repository that was loaded must come from this worktree
    for name, m in list(sys.modules.items()):
        f = getattr(m, "__file__", None) or ""
        if f.startswith(os.path.expanduser("~/isr_project/")) and not os.path.abspath(f).startswith(code + os.sep):
            raise RuntimeError("module %s was loaded from %s, not from %s" % (name, f, code))


# ---------------------------------------------------------------- running a commit

def lpi_periods(ch):
    per = set(int(p) for p in BENCH["lpi"][ch])
    a, b = BENCH["fit_lpi"][ch]
    return sorted(per | set(range(a, b + 1)))


def jobs_for(stage):
    if stage == "lpi":
        return [("lpi", ch, per[i:i + CHUNK]) for ch in ("zenith-l", "misa-l")
                for per in [lpi_periods(ch)] for i in range(0, len(per), CHUNK)]
    if stage == "rd300":
        return [("rd300", ch, [p]) for ch, ps in BENCH["range_doppler_300"].items() for p in ps]
    if stage == "fit_lpi":
        return [("fit_lpi", ch, []) for ch in BENCH["fit_lpi"]]
    if stage == "fit_lp":
        return [("fit_lp", ch, []) for ch, on in BENCH["fit_lp"].items() if on]


def environment_key():
    """Hash of everything besides the commit that the outputs depend on."""
    h = hashlib.sha256()
    for f in (os.path.abspath(__file__), os.path.join(HERE, "benchmark.json")):
        h.update(open(f, "rb").read())
    h.update(subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True).stdout.encode())
    h.update(sys.version.encode())
    for f in sorted(glob.glob(os.path.join(BENCH["table_dir"], "*"))):
        st = os.stat(f)
        h.update(("%s %d %d" % (os.path.basename(f), st.st_size, int(st.st_mtime))).encode())
    return h.hexdigest()[:12]


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
                   HDF5_USE_FILE_LOCKING="FALSE", MPLBACKEND="Agg", PYTHONUNBUFFERED="1")
        running.append((subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env), (kind, ch, per)))
    for p, job in running:
        if p.wait() != 0:
            failed.append(job)
    return failed


def run_commit(sha, key, n_jobs):
    """Outputs of commit `sha` on the benchmark, cached under (sha, key)."""
    out = os.path.join(ROOT, "runs", "%s-%s" % (sha, key))
    os.makedirs(os.path.join(ROOT, "runs"), exist_ok=True)
    lock = open(out + ".lock", "w")
    fcntl.flock(lock, fcntl.LOCK_EX)          # one run per commit at a time
    done = os.path.join(out, "DONE")
    if os.path.exists(done) and json.load(open(done)).get("key") == key:
        print("%s: cached outputs in %s" % (sha[:7], out), flush=True)
        return out, []
    if os.path.exists(out):
        shutil.rmtree(out)                    # nothing from an earlier, failed run may be compared
    code = os.path.join(ROOT, "worktrees", sha)
    if not os.path.isdir(code):
        os.makedirs(os.path.dirname(code), exist_ok=True)
        git("worktree", "add", "--detach", code, sha)
    logdir = os.path.join(out, "logs")
    os.makedirs(logdir)
    t0 = time.time()
    print("%s: running the benchmark with %d parallel jobs" % (sha[:7], n_jobs), flush=True)
    failed = run_jobs(jobs_for("lpi") + jobs_for("rd300"), code, out, logdir, n_jobs)
    if not failed:
        failed = run_jobs(jobs_for("fit_lpi") + jobs_for("fit_lp"), code, out, logdir, n_jobs)
    if not failed:
        json.dump(dict(sha=sha, key=key, seconds=time.time() - t0,
                       when=datetime.datetime.now(datetime.UTC).isoformat()), open(done, "w"))
        git("worktree", "remove", "--force", code)
    print("%s: done in %.0f s, %d failed jobs" % (sha[:7], time.time() - t0, len(failed)), flush=True)
    return out, failed


# ---------------------------------------------------------------- expected outputs

def expected(out):
    """The outputs the benchmark must produce; returns a list of what is missing."""
    from digital_rf import DigitalMetadataReader
    b0 = DigitalMetadataReader("%s/metadata/id_metadata" % DATA).get_bounds()[0]
    t = lambda ai: int((b0 + ai * 10 * 10**6) / 1e6)
    missing = []
    for ch, per in BENCH["lpi"].items():
        for ai in per:
            f = os.path.join(out, "lpi_%d" % LPI["rg"], ch, "lpi-%d.h5" % t(int(ai)))
            if not os.path.exists(f):
                missing.append(os.path.relpath(f, out))
    for ch, (a, b) in BENCH["fit_lpi"].items():
        k = sum(os.path.exists(os.path.join(out, "lpi_%d" % LPI["rg"], ch, "lpi-%d.h5" % t(ai))) for ai in range(a, b + 1))
        if k < MIN_STRETCH_FILES:
            missing.append("%s: only %d LPI files in the fit stretch (need %d)" % (ch, k, MIN_STRETCH_FILES))
        pp = [int(os.path.basename(f)[3:-3]) for f in glob.glob(os.path.join(out, "lpi_%d" % LPI["rg"], ch, "pp-*.h5"))]
        if not any(t(a) - 300 <= x <= t(b) for x in pp):
            missing.append("%s: no fit_lpi output for the fit stretch" % ch)
    for ch, per in BENCH["range_doppler_300"].items():
        for ai in per:
            f = os.path.join(out, "range_doppler_300_outlier", ch, "il_%d.h5" % t(ai))
            if not os.path.exists(f):
                missing.append(os.path.relpath(f, out))
    for ch, on in BENCH["fit_lp"].items():
        if on and not glob.glob(os.path.join(out, "range_doppler_300_outlier", ch, "pp-*.h5")):
            missing.append("%s: no fit_lp output" % ch)
    return missing


# ---------------------------------------------------------------- comparing

def datasets(h):
    """Every dataset of an HDF5 file, by its full path."""
    out = {}
    h.visititems(lambda name, obj: out.__setitem__(name, obj) if hasattr(obj, "shape") else None)
    return out


def identical(a, b):
    """Byte-for-byte equality, including dtype and shape."""
    if a.dtype != b.dtype or a.shape != b.shape:
        return False
    if a.dtype.kind == "O":
        return all(x == y for x, y in zip(a.ravel().tolist(), b.ravel().tolist()))
    return a.tobytes() == b.tobytes()


def same_attrs(x, y):
    kx, ky = set(x.attrs.keys()), set(y.attrs.keys())
    return kx == ky and all(identical(n.asarray(x.attrs[k]), n.asarray(y.attrs[k])) for k in kx)


def describe(name, a, b, base_sets):
    """How a dataset changed, in the most meaningful unit available."""
    try:
        if a.dtype != b.dtype:
            return "dtype %s -> %s" % (a.dtype, b.dtype)
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
        var = {"acfs_e": "acfs_var", "acfs_g": "acfs_var"}.get(name.split("/")[-1])
        s = n.sqrt(n.abs(n.asarray(base_sets[var][()]))) if var in base_sets else None
        if s is not None and s.shape == a.shape and both.any():
            z = (n.abs(b - a) / s)[both]
            z = z[n.isfinite(z)]
            if len(z):
                txt.append("|change| / sigma: median %.3g, max %.3g" % (n.median(z), n.max(z)))
        elif both.any():
            d = n.abs(b - a)[both]
            scale = n.abs(a)[both]
            nz = scale > 0
            if nz.any():
                r = d[nz] / scale[nz]
                txt.append("relative change: median %.3g, max %.3g" % (n.median(r), n.max(r)))
            if (~nz).any() and (d[~nz] > 0).any():
                txt.append("%d values changed from zero" % (d[~nz] > 0).sum())
        if not txt:
            txt.append("bits differ (NaN payload, signed zero or roundoff below display)")
        return "; ".join(txt)
    except Exception as e:
        return "changed (could not describe: %s)" % e


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
            db, dt = datasets(hb), datasets(ht)
            changed = False
            if not same_attrs(hb, ht):
                rows.append((f, "/", "file attributes differ"))
                changed = True
            for k in sorted(set(db) ^ set(dt)):
                rows.append((f, k, "dataset only in %s" % ("base" if k in db else "the commit under review")))
                changed = True
            for k in sorted(set(db) & set(dt)):
                a, b = n.asarray(db[k][()]), n.asarray(dt[k][()])
                if not identical(a, b):
                    rows.append((f, k, describe(k, a, b, db)))
                    changed = True
                elif not same_attrs(db[k], dt[k]):
                    rows.append((f, k, "attributes differ"))
                    changed = True
            n_ident += not changed
    return rows, n_ident, len(fb | ft)


# ---------------------------------------------------------------- main

def run(args):
    base = git("rev-parse", "--verify", args.base + "^{commit}")
    test = git("rev-parse", "--verify", args.target + "^{commit}")
    print("base %s (%s), under review %s (%s)" % (base[:7], args.base, test[:7], args.target), flush=True)
    if test == base:
        raise Invalid("the commit under review is the base itself: nothing to test")
    if subprocess.run(["git", "merge-base", "--is-ancestor", base, test], cwd=REPO).returncode != 0:
        raise Invalid("%s does not contain %s: merge or rebase %s into the branch first, so that only its own "
                      "changes are compared" % (args.target, args.base, args.base))
    for line in git("worktree", "list", "--porcelain").split("\n\n"):
        if "branch refs/heads/%s" % args.target in line:
            wt = line.split("\n")[0].split(" ", 1)[1]
            if git("status", "--porcelain", "--untracked-files=no", cwd=wt):
                print("WARNING: the worktree %s of %s has uncommitted changes; they are not tested" % (wt, args.target))
    key = environment_key()
    base_out, fail_b = run_commit(base, key, args.jobs)
    if fail_b:
        raise Invalid("the base failed on the benchmark: %s; see %s/logs" % (fail_b, base_out))
    test_out, fail_t = run_commit(test, key, args.jobs)
    miss = {"base": expected(base_out), "under review": expected(test_out)}
    rows, n_ident, n_files = compare(base_out, test_out)
    lines = ["# Regression benchmark: %s (%s) against %s (%s)" % (args.target, test[:7], args.base, base[:7]), "",
             "Run %s UTC, environment key %s. %d output files compared; %d bit-identical." % (
                 datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%d %H:%M"), key, n_files, n_ident), ""]
    for kind, ch, per in fail_t:
        lines.append("- FAILED job (under review): %s %s periods %s; see %s/logs" % (kind, ch, per, test_out))
    for side, m in miss.items():
        for x in m:
            lines.append("- MISSING (%s): %s" % (side, x))
    if n_files == 0:
        lines.append("- NO OUTPUTS were compared.")
    if rows:
        lines += ["", "| file | dataset | change |", "|---|---|---|"] + ["| %s | %s | %s |" % r for r in rows]
    elif n_files:
        lines.append("Every output is bit-identical.")
    txt = "\n".join(lines) + "\n"
    print(txt)
    os.makedirs(os.path.join(ROOT, "reports"), exist_ok=True)
    name = "%s-%s_vs_%s-%s.md" % (args.target.replace("/", "_"), test[:7], args.base.replace("/", "_"), base[:7])
    open(os.path.join(ROOT, "reports", name), "w").write(txt)
    print("report: %s" % os.path.join(ROOT, "reports", name))
    if fail_t or n_files == 0 or any(miss.values()):
        return 2
    return 1 if rows else 0


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
    try:
        status = run(args)
    except Invalid as e:
        print("INVALID RUN: %s" % e)
        status = 2
    except Exception:
        traceback.print_exc()
        print("INVALID RUN: the tool failed")
        status = 2
    sys.exit(status)


if __name__ == "__main__":
    main()
