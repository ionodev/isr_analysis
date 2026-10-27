#!/usr/bin/env python3
"""
Altitude-time plots of n_e, T_e, T_i and v_i with their uncertainties and
quality flags (weeks 13-15 of the project plan), for one LPI variant's fits.

For one day it draws two figures:
    <name>_params.png       the four parameters; gates that are not good
                            (fit_quality: failed, at a fit limit, or noisy)
                            in grey; gates with a strong catalogued
                            satellite echo outlined in black dots
    <name>_uncertainty.png  dne/ne, dTe, dTi and dvi
and, given a second variant's fits, <name>_difference.png: the difference of
n_e and T_e between the two in the first variant's standard deviation (with
against without the satellite columns).

The zenith antenna transmits in cycles of about three minutes between the
MISA cycles, so the fits come in groups of three one-minute windows with
gaps; the time axis keeps the gaps. Altitude is taken as range (the zenith
beam is at 88 degrees elevation, Memo 37: altitude = 0.9995 range).

With --raw, the parameter figure shows every fitted value as it is: no
flagged gates greyed, no echo dots, each fit's column extended to the next
fit however long the gap, and colour maps without black (written as
<name>_params_raw.png).

Usage:
    python3 eclipse_plots.py [--raw] <fits_dir> <day YYYY-MM-DD> <out_prefix> [<echoes.npy>] [<other_fits_dir>]
"""

import datetime
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as n

import fit_quality as fq

# the partial eclipse at Westford, 8 April 2024 (contacts of nearby Boston)
ECLIPSE = (1712600160, 1712608680)
ALT = (150, 700)


MAX_COLUMN_S = 420


def grid(F, values, t_lo=None, t_hi=None):
    """One column per fit, from its start to the next fit's start (at most
    MAX_COLUMN_S, the usual way to draw ISR profiles), with an empty column
    over longer gaps. Returns (column edges in unix s, values per column)."""
    order = n.argsort(F["t0"])
    t0 = F["t0"][order]
    edges, cols = [t0[0]], []
    for i, k in enumerate(order):
        nxt = t0[i + 1] if i + 1 < len(t0) else t0[i] + 60.0
        right = min(nxt, t0[i] + MAX_COLUMN_S)
        cols.append(values[k])
        edges.append(right)
        if right < nxt:
            cols.append(n.full(values.shape[1], n.nan))
            edges.append(nxt)
    return n.array(edges), n.array(cols)


def hours(t, day0):
    return (n.asarray(t) - day0) / 3600.0


def panel(ax, tg, r, z, day0, cmap, vmin, vmax, label, log=False, overlay=None, hatch=None):
    zz = n.log10(z) if log else z
    m = ax.pcolormesh(hours(tg, day0), r, zz.T, cmap=cmap, vmin=vmin, vmax=vmax, shading="flat", rasterized=True)
    if overlay is not None:
        ax.pcolormesh(hours(tg, day0), r, n.where(overlay, 1.0, n.nan).T, cmap="Greys", vmin=0, vmax=2,
                      shading="flat", rasterized=True)
    if hatch is not None:
        tt, rr = n.nonzero(n.asarray(hatch) > 0.5)
        rc = 0.5 * (r[:-1] + r[1:])
        ax.plot(hours(0.5 * (tg[tt] + tg[tt + 1]), day0), rc[rr], "k.", ms=0.4)
    cb = plt.colorbar(m, ax=ax, pad=0.01)
    cb.set_label(label)
    ax.set_ylim(*ALT)
    ax.set_ylabel("altitude (km)")
    if ECLIPSE[0] < tg[-1] and ECLIPSE[1] > tg[0]:
        for x in ECLIPSE:
            ax.axvline(hours(x, day0), color="k", lw=1.0, ls="--")


def main():
    global MAX_COLUMN_S
    raw = "--raw" in sys.argv
    argv = [a for a in sys.argv if a != "--raw"]
    if raw:
        MAX_COLUMN_S = n.inf
    fdir, day, prefix = argv[1:4]
    echoes = n.load(argv[4]) if len(argv) > 4 and argv[4] != "-" else None
    other = argv[5] if len(argv) > 5 else None
    day0 = datetime.datetime.fromisoformat(day).replace(tzinfo=datetime.timezone.utc).timestamp()
    t_lo, t_hi = day0 + 11 * 3600, day0 + 23 * 3600
    F = fq.load_fits(fdir, (t_lo, t_hi))
    f = fq.flags(F, echoes)
    bad = ~fq.good(f)
    sat = (f & fq.SATELLITE) > 0
    r = F["r"]
    redges = n.concatenate([r - 0.5 * (r[1] - r[0]), [r[-1] + 0.5 * (r[1] - r[0])]])
    G = {k: grid(F, F[k]) for k in ("ne", "Te", "Ti", "vi", "dne", "dTe", "dTi", "dvi")}
    tg = G["ne"][0]
    gb = n.nan_to_num(grid(F, bad.astype(float))[1]) > 0.5
    gs = n.nan_to_num(grid(F, sat.astype(float))[1]) > 0.5
    title = "%s, zenith, %s" % (day, os.path.basename(os.path.normpath(fdir)))
    clean = lambda k: n.where(gb, n.nan, G[k][1])

    fig, ax = plt.subplots(4, 1, figsize=(11, 12), sharex=True)
    if raw:
        # every fitted value as it is: no flags, no dots, no black
        panel(ax[0], tg, redges, G["ne"][1], day0, "viridis", 10.8, 12.2, "log$_{10}$ $n_e$ (m$^{-3}$)", log=True)
        panel(ax[1], tg, redges, G["Te"][1], day0, "Spectral_r", 0, 3500, "$T_e$ (K)")
        panel(ax[2], tg, redges, G["Ti"][1], day0, "Spectral_r", 0, 2000, "$T_i$ (K)")
        panel(ax[3], tg, redges, G["vi"][1], day0, "RdBu_r", -150, 150, "$v_i$ (m/s, positive away)")
        ax[0].set_title(title + ": every fitted value, unflagged; dashed = eclipse")
    else:
        panel(ax[0], tg, redges, clean("ne"), day0, "viridis", 10.8, 12.2, "log$_{10}$ $n_e$ (m$^{-3}$)", log=True, overlay=gb, hatch=gs)
        panel(ax[1], tg, redges, clean("Te"), day0, "inferno", 500, 3500, "$T_e$ (K)", overlay=gb, hatch=gs)
        panel(ax[2], tg, redges, clean("Ti"), day0, "inferno", 500, 2000, "$T_i$ (K)", overlay=gb, hatch=gs)
        panel(ax[3], tg, redges, clean("vi"), day0, "RdBu_r", -150, 150, "$v_i$ (m/s, positive away)", overlay=gb, hatch=gs)
        ax[0].set_title(title + ": grey = flagged (failed, fit limit, noisy); dots = catalogued satellite echo; dashed = eclipse")
    ax[-1].set_xlabel("hours (UTC)")
    fig.tight_layout()
    fig.savefig(prefix + ("_params_raw.png" if raw else "_params.png"), dpi=110)
    if raw:
        return

    fig, ax = plt.subplots(4, 1, figsize=(11, 12), sharex=True)
    with n.errstate(invalid="ignore", divide="ignore"):
        panel(ax[0], tg, redges, n.abs(G["dne"][1] / G["ne"][1]), day0, "magma", 0, 0.5, "$\\delta n_e/n_e$")
    panel(ax[1], tg, redges, G["dTe"][1], day0, "magma", 0, 1000, "$\\delta T_e$ (K)")
    panel(ax[2], tg, redges, G["dTi"][1], day0, "magma", 0, 500, "$\\delta T_i$ (K)")
    panel(ax[3], tg, redges, G["dvi"][1], day0, "magma", 0, 100, "$\\delta v_i$ (m/s)")
    ax[0].set_title(title + ": uncertainties")
    ax[-1].set_xlabel("hours (UTC)")
    fig.tight_layout()
    fig.savefig(prefix + "_uncertainty.png", dpi=110)

    if other:
        O = fq.load_fits(other, (t_lo, t_hi))
        fo = fq.flags(O)
        if O is not None and len(O["t0"]) == len(F["t0"]):
            both = fq.good(f) & fq.good(fo)
            fig, ax = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True)
            for a, k, lab in ((ax[0], "ne", "$n_e$"), (ax[1], "Te", "$T_e$")):
                d = n.where(both, (O[k] - F[k]) / F["d" + k], n.nan)
                panel(a, tg, redges, grid(F, d)[1], day0, "RdBu_r", -3, 3,
                      "%s: other $-$ this ($\\sigma$)" % lab, hatch=gs)
            ax[0].set_title("%s: %s minus %s" % (day, os.path.basename(os.path.normpath(other)),
                                                 os.path.basename(os.path.normpath(fdir))))
            ax[-1].set_xlabel("hours (UTC)")
            fig.tight_layout()
            fig.savefig(prefix + "_difference.png", dpi=110)
    print("good gates at %d-%d km: %.2f" % (ALT[0], ALT[1], fq.good(f)[:, (r >= ALT[0]) & (r <= ALT[1])].mean()))


if __name__ == "__main__":
    main()
