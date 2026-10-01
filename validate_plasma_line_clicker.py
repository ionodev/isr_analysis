#!/usr/bin/env python3
"""
plasma_line_clicker.build_spectra (gates chosen by range, 185-560 km) against
the previous version of the script (gate indices 50:150), on a recording: the
per-minute spectra of both antennas must be bit identical.

The previous version is read from git (main before the change, 0e93f51) and run
with its interactive click_spec replaced by a function that keeps the arrays.
The two run as separate processes, at the same time.

    python3 validate_plasma_line_clicker.py <recording dir>
"""
import multiprocessing
import subprocess
import sys

import numpy as n

OLD_COMMIT = "0e93f51"


def old(dirname, q):
    src = subprocess.run(["git", "show", "%s:plasma_line_clicker.py" % OLD_COMMIT], capture_output=True, text=True, check=True).stdout
    src = src.replace("from mpl_point_clicker import clicker", "")
    got = {}

    def capture(t, s, freqs, dirname, calname="cal.h5"):
        got[calname] = (n.array(t), n.array(s))
    ns = {"__name__": "old_clicker"}
    sys.argv = ["plasma_line_clicker.py", dirname]
    code = compile(src, "old_plasma_line_clicker.py", "exec")
    # run the definitions, then swap click_spec before the module-level calls use it
    src_def, src_run = src.split("dirname=sys.argv[1]", 1)
    exec(compile(src_def, "old_def", "exec"), ns)
    ns["click_spec"] = capture
    exec(compile("dirname=sys.argv[1]" + src_run, "old_run", "exec"), ns)
    q.put(got)


def new(dirname, q):
    import plasma_line_clicker as plc
    zt, zs, mt, ms, f = plc.build_spectra(dirname)
    q.put({"zenithcal": (n.array(zt), n.vstack(zs)), "misacal": (n.array(mt), n.vstack(ms)) if len(ms) else None})


def main():
    dirname = sys.argv[1]
    ctx = multiprocessing.get_context("fork")
    qa, qb = ctx.Queue(), ctx.Queue()
    pa, pb = ctx.Process(target=old, args=(dirname, qa)), ctx.Process(target=new, args=(dirname, qb))
    pa.start(); pb.start()
    a, b = qa.get(), qb.get()
    pa.join(); pb.join()
    ok = True
    for k in ("zenithcal", "misacal"):
        ta, sa = a[k]
        tb, sb = b[k]
        same = n.array_equal(ta, tb) and n.array_equal(sa, sb, equal_nan=True)
        ok &= same
        print("%s: %d minutes, times and spectra bit identical: %s" % (k, len(ta), same))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
