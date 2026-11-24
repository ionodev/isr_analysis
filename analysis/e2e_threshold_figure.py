#!/usr/bin/env python3
"""
Figure of the end-to-end test at three detection thresholds: the detection
probability per pulse against the pulse's SNR, and the ACF shift at the
echo's gate by the pass's peak SNR for tracked passes (v1, P_FA 1e-10) and
for confirmed passes at P_FA 1e-6 and 1e-4.

    python3 e2e_threshold_figure.py <v1_dir> <pfa1e-6_dir> <pfa1e-4_dir> <figure.png>
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as n


def curve(d):
    c = json.load(open(os.path.join(d, "detections.json")))["curve"]
    return [0.5 * (x["snr_lo"] + x["snr_hi"]) for x in c], [x["detected"] for x in c]


def shifts(d, name):
    for r in json.load(open(os.path.join(d, "shift_by_peak.json"))):
        if r["name"] == name:
            k = sorted(r["shift_by_peak"], key=float)
            return [float(x) for x in k], [r["shift_by_peak"][x] for x in k]


def main():
    v1, p6, p4, out = sys.argv[1:5]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for d, lab, c in ((v1, r"$P_{\rm FA}=10^{-10}$ (v1)", "C0"), (p6, r"$10^{-6}$", "C1"), (p4, r"$10^{-4}$", "C2")):
        x, y = curve(d)
        ax[0].plot(x, y, "o-", color=c, label=lab)
    ax[0].set_xlabel("matched-filter SNR of the pulse (dB)")
    ax[0].set_ylabel("fraction of pulses detected")
    ax[0].set_title("Detection probability")
    ax[0].legend(fontsize=8)
    for (d, name), lab, c in (((v1, "none"), "no satellite columns", "0.5"),
                              ((v1, "realistic_tracked_ext2_min1"), r"tracked, $10^{-10}$ (v1)", "C0"),
                              ((p6, "realistic_confirmed"), r"confirmed, $10^{-6}$", "C1"),
                              ((p4, "realistic_confirmed"), r"confirmed, $10^{-4}$", "C2"),
                              ((v1, "perfect"), "every pulse (best case)", "k")):
        x, y = shifts(d, name)
        ax[1].plot(x, y, "o-", color=c, label=lab, ls=":" if c in ("0.5", "k") else "-")
    ax[1].axhspan(-1, 1, color="0.9")
    ax[1].set_xlabel("peak SNR of the pass (dB)")
    ax[1].set_ylabel(r"ACF shift at the echo's gate ($\sigma$)")
    ax[1].set_title("Bias left by the pass")
    ax[1].set_ylim(-1.5, 13)
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)


if __name__ == "__main__":
    main()
