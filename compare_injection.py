#!/usr/bin/env python3
"""
Bias of the LPI variants against a known truth: the injection test.

For each variant and scenario of inject_satellites.py, the injected run is
compared with the same variant's run on the same periods without the
injection. The raw noise is the same in both, so their difference is what
the synthetic satellite does to the plasma estimate; a perfect method would
leave it at zero. Measures, on the real part of lags 1-10 in standard
deviations of the injected run's real part (sqrt(acfs_var/2)):

    bias_echo        median |difference| at the echo's gate and the gates
                     next to it
    signed_echo      median signed difference there (a systematic shift)
    bias_echo_clean  as bias_echo, in the clean run's standard deviation:
                     the same yardstick for every variant (a variant whose
                     error bars grow with the echo looks less biased in
                     its own units)
    leak_profile     median |difference| at 150-700 km, 5 or more gates away
                     from the echo
    leak_frac3       fraction of those cells beyond 3 sigma
    sigma_ratio_echo median error bar with / without the echo, at the echo
    sigma_ratio_far  the same away from the echo

Usage:
    python3 compare_injection.py <injection_dir> <clean_variants_dir> <out_dir>
"""

import glob
import json
import os
import sys

import h5py
import numpy as n

T0_RECORDING = 1712484000
PERIOD_S = 10
LAGS = slice(1, 11)


def load(d, ai):
    f = os.path.join(d, "lpi_60", "zenith-l", "lpi-%d.h5" % (T0_RECORDING + PERIOD_S * ai))
    if not os.path.exists(f):
        return None
    with h5py.File(f, "r") as h:
        return h["acfs_e"][()], h["acfs_var"][()], h["rgs_km"][()]


def compare(inj_dir, clean_dir):
    meta = json.load(open(os.path.join(inj_dir, "injection.json")))
    periods = [p["period"] for p in meta["periods"]]
    echo, far, s_echo, s_far, signed, echo_c = [], [], [], [], [], []
    for ai in periods:
        a, b = load(inj_dir, ai), load(clean_dir, ai)
        if a is None or b is None:
            continue
        (xi, vi, r), (xc, vc, _) = a, b
        g = int(n.argmin(n.abs(r - meta["range_km"])))
        si = n.sqrt(vi[:, LAGS] / 2.0)
        sc = n.sqrt(vc[:, LAGS] / 2.0)
        z = (xi[:, LAGS] - xc[:, LAGS]).real / si
        near = n.zeros(len(r), bool)
        near[max(0, g - 1):g + 2] = True
        away = (n.abs(n.arange(len(r)) - g) >= 5) & (r >= 150) & (r <= 700)
        echo.append(z[near].ravel())
        echo_c.append(((xi[:, LAGS] - xc[:, LAGS]).real / sc)[near].ravel())
        far.append(z[away].ravel())
        s_echo.append((si / sc)[near].ravel())
        s_far.append((si / sc)[away].ravel())
    cat = lambda x: n.concatenate(x) if x else n.array([n.nan])
    echo, far, s_echo, s_far, echo_c = cat(echo), cat(far), cat(s_echo), cat(s_far), cat(echo_c)
    return dict(variant=meta["variant"], scenario=meta["scenario"], snr_db=meta["snr_db"],
                pulses_injected=meta["pulses_injected"], periods=len(periods),
                bias_echo=float(n.nanmedian(n.abs(echo))),
                signed_echo=float(n.nanmedian(echo)),
                bias_echo_clean=float(n.nanmedian(n.abs(echo_c))),
                leak_profile=float(n.nanmedian(n.abs(far))),
                leak_frac3=float(n.nanmean(n.abs(far) > 3)),
                sigma_ratio_echo=float(n.nanmedian(s_echo)),
                sigma_ratio_far=float(n.nanmedian(s_far)))


def main():
    inj_root, clean_root, outd = sys.argv[1:4]
    os.makedirs(outd, exist_ok=True)
    rows = []
    for d in sorted(glob.glob(os.path.join(inj_root, "*", "v*"))):
        v = os.path.basename(d)
        rows.append(compare(d, os.path.join(clean_root, v)))
        print(rows[-1], flush=True)
    json.dump(rows, open(os.path.join(outd, "injection_measures.json"), "w"), indent=1)
    keys = ["snr_db", "bias_echo", "signed_echo", "bias_echo_clean", "leak_profile", "leak_frac3", "sigma_ratio_echo", "sigma_ratio_far"]
    with open(os.path.join(outd, "injection_measures.md"), "w") as f:
        f.write("| scenario | variant | " + " | ".join(keys) + " |\n|" + "---|" * (len(keys) + 2) + "\n")
        for m in rows:
            f.write("| %s | %s | %s |\n" % (m["scenario"], m["variant"],
                                           " | ".join("%.3g" % m[k] for k in keys)))


if __name__ == "__main__":
    main()
