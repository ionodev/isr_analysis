#!/usr/bin/env python3
"""
Plot raw (undecoded, unfiltered) complex voltage samples, one figure per
interpulse period.

The interpulse period boundaries are taken from the transmit pulse metadata
(metadata/id_metadata, "sweepid"), so every trace starts exactly at a transmit
pulse and ends where the next one begins. Two PNGs are written per IPP per
channel:

    <outdir>/<channel>/raw voltage/ipp-<index>-<unix_t>.png   |z|, real, imag
    <outdir>/<channel>/power/ipp-<index>-<unix_t>.png         |z|^2 in dB

plus, when more than one channel is plotted, all channels' power in one figure:

    <outdir>/combined/ipp-<index>-<unix_t>.png                |z|^2 in dB
    <outdir>/combined/close_up/ipp-<index>-<unix_t>.png       the same, transmit pulse only
    <outdir>/combined/tx_delay-<unix_t>.png                   matched filter delay vs IPP

--fix-tx-delay shifts each echo channel onto the tx-h time base by its measured
channel delay (see tx_delay.py) and writes the same set of plots to
IPP_offset_fix instead of IPP. The delay plot then shows the residual, which
should sit at zero.

The dataset is read from the pipeline config (config/millstone_eclipse2024.json
next to this script): "data_dir" says where the raw data is and "output_dir",
if present, is where the plots go.

Usage:
    python3 plot_raw_voltage.py [options]

Examples:
    # first 10 IPPs of every channel -> <output_dir>/IPP/<channel>/{raw voltage,power}
    python3 plot_raw_voltage.py

    # 20 IPPs of one channel only, starting 100 IPPs into the recording
    python3 plot_raw_voltage.py --channel zenith-l --n-ipp 20 --start-ipp 100
"""

import os
import sys
import json
import argparse

import numpy as n

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from digital_rf import DigitalRFReader, DigitalMetadataReader

import stuffr
from tx_delay import TX_CHANNEL, CODED_SWEEPIDS, matched_filter_delay, fractional_shift, tmm

# USRP DC offset caused by truncation instead of rounding in the FPGA.
# Only removed when --remove-dc is given: the default is genuinely raw data.
Z_DC = n.complex64(-0.212 - 0.221j)


# pipeline config, resolved next to this script so the cwd does not matter
CONFIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config", "millstone_eclipse2024.json")


def load_config():
    """Raw data directory and plot output base from the pipeline config."""
    if not os.path.isfile(CONFIG):
        print("ERROR: config %s not found" % (CONFIG))
        sys.exit(1)
    with open(CONFIG) as f:
        cfg = json.load(f)

    dirname = cfg.get("data_dir", None)
    if dirname is None:
        print("ERROR: no \"data_dir\" in %s" % (CONFIG))
        sys.exit(1)
    if not os.path.isdir(dirname):
        print("ERROR: data_dir %s does not exist" % (dirname))
        sys.exit(1)

    # plots go under output_dir when the config sets one, else next to the script
    outbase = cfg.get("output_dir", os.path.dirname(os.path.abspath(__file__)))
    return dirname, outbase


def get_ipps(dirname, n_ipp, start_ipp, t0_unix):
    """
    Return (keys, sweepids, ipp_len) for n_ipp interpulse periods.

    keys are sample indices (1 MHz) of the transmit pulse starts and ipp_len is
    the median spacing between consecutive pulses, in samples.
    """
    id_read = DigitalMetadataReader("%s/metadata/id_metadata" % (dirname))
    idsr = int(id_read.get_samples_per_second())
    idb = id_read.get_bounds()

    if t0_unix is not None:
        i0 = int(t0_unix * idsr)
        if i0 < idb[0] or i0 > idb[1]:
            print("ERROR: --t0 %1.3f is outside the data (%s to %s)"
                  % (t0_unix, stuffr.unix2datestr(float(idb[0]) / idsr), stuffr.unix2datestr(float(idb[1]) / idsr)))
            sys.exit(1)
    else:
        i0 = idb[0]

    # read a generous window and cut it down: we don't know the IPP length yet.
    # 100 ms per IPP is a safe upper bound for any Millstone mode.
    read_len = int((start_ipp + n_ipp + 2) * 0.1 * idsr)
    sid = id_read.read(i0, min(i0 + read_len, idb[1]), "sweepid")

    keys = n.array(list(sid.keys()), dtype=n.int64)
    keys.sort()
    if len(keys) < 2:
        print("ERROR: fewer than 2 transmit pulses found after %s" % (stuffr.unix2datestr(float(i0) / idsr)))
        sys.exit(1)

    ipp_len = int(n.median(n.diff(keys)))

    if start_ipp + n_ipp > len(keys):
        print("WARNING: only %d interpulse periods available, plotting those" % (max(0, len(keys) - start_ipp)))
    keys = keys[start_ipp:(start_ipp + n_ipp)]
    sweepids = n.array([sid[k] for k in keys], dtype=int)

    return keys, sweepids, ipp_len


def new_ipp_figure(sweepid, channel, ipp_idx, t_unix, what):
    """Figure with the transmit pulse and noise injection windows shaded."""
    fig, ax = plt.subplots(figsize=(12, 5))

    # transmit pulse (red) and noise injection (green) windows
    if sweepid in tmm:
        t = tmm[sweepid]
        ax.axvspan(t["tx0"], t["tx1"], color="red", alpha=0.12, lw=0, zorder=0)
        ax.axvspan(t["noise0"], t["noise1"], color="green", alpha=0.12, lw=0, zorder=0)

    ax.set_xlabel(r"Time since transmit pulse start ($\mu$s)")
    ax.set_title("%s %s, IPP %d, sweepid %d, %s\n(red = transmit pulse, green = noise injection)%s"
                 % (channel, what, ipp_idx, sweepid, stuffr.unix2datestr(t_unix), TITLE_NOTE), fontsize=11)
    return fig, ax


def plot_ipp_voltage(z, t_us, sweepid, channel, ipp_idx, t_unix, ofname):
    """One interpulse period: |z|, real and imaginary in a single linear plot."""
    fig, ax = new_ipp_figure(sweepid, channel, ipp_idx, t_unix, "raw voltage")

    ax.plot(t_us, z.real, color="C0", lw=0.5, label="real")
    ax.plot(t_us, z.imag, color="C1", lw=0.5, label="imag")
    ax.plot(t_us, n.abs(z), color="k", lw=0.5, label="|z|")

    ax.set_xlim(0, t_us[-1])
    ax.set_ylabel("Raw voltage (ADC units)")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(ofname, dpi=150)
    plt.close(fig)


def power_db(z):
    """Received power in dB. Samples that are exactly zero would be -inf, drop those."""
    pwr = n.abs(z) ** 2.0
    pwr[pwr <= 0] = n.nan
    return 10.0 * n.log10(pwr)


def plot_ipp_power(z, t_us, sweepid, channel, ipp_idx, t_unix, ofname):
    """One interpulse period: received power 10*log10(|z|^2), in dB."""
    fig, ax = new_ipp_figure(sweepid, channel, ipp_idx, t_unix, "power")

    ax.plot(t_us, power_db(z), color="k", lw=0.5, label=r"$10\log_{10}|z|^2$")

    ax.set_xlim(0, t_us[-1])
    ax.set_ylabel("Power (dB, uncalibrated)")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(ofname, dpi=150)
    plt.close(fig)


# how much to show on either side of the transmit pulse in the close up plots
CLOSEUP_MARGIN_US = 100.0

# subdirectory of output_dir, with and without the channel delay correction
OUT_SUBDIR = "IPP"
OUT_SUBDIR_FIXED = "IPP_offset_fix"

# appended to figure titles, set when the correction is applied
TITLE_NOTE = ""


def tx_window(sweepid, ipp_len, margin=0.0):
    """Sample indices of the transmit pulse within an interpulse period."""
    if sweepid not in tmm:
        return None
    i0 = int(max(0, tmm[sweepid]["tx0"] - margin))
    i1 = int(min(ipp_len, tmm[sweepid]["tx1"] + margin))
    return i0, i1


def measure_channel_delays(d_il, dirname, key0, channels, srs, ipp_len,
                           n_coded=20, oversample=100, max_lag_us=100.0):
    """
    Median matched filter delay of each echo channel relative to tx-h, measured
    from the coded pulses starting at key0.

    Only coded pulses are used: an uncoded long pulse has no sharp correlation
    peak and biases the estimate by ~0.2 us. See tx_delay.py.
    """
    id_read = DigitalMetadataReader("%s/metadata/id_metadata" % (dirname))
    idb = id_read.get_bounds()
    sid = id_read.read(int(key0), min(int(key0) + 20 * n_coded * ipp_len, idb[1]), "sweepid")
    keys = sorted(sid.keys())

    delays = {}
    for channel in channels:
        if channel == TX_CHANNEL:
            continue
        ds = []
        for key in keys:
            if len(ds) >= n_coded:
                break
            if sid[key] not in CODED_SWEEPIDS:
                continue
            t = tmm[sid[key]]
            n_s = t["tx1"] - t["tx0"]
            try:
                z_tx = d_il.read_vector_1d(int(key) + t["tx0"], n_s, TX_CHANNEL).astype("c8", casting="unsafe", copy=False)
                z_echo = d_il.read_vector_1d(int(key) + t["tx0"], n_s, channel).astype("c8", casting="unsafe", copy=False)
            except Exception:
                continue
            d_samples, amp = matched_filter_delay(z_tx, z_echo, oversample=oversample,
                                                  max_lag=max_lag_us * srs[channel] / 1e6)
            ds.append(d_samples / srs[channel] * 1e6)
        if len(ds) == 0:
            print("could not measure the delay of %s, leaving it uncorrected" % (channel))
            continue
        delays[channel] = float(n.median(n.array(ds)))
        print("%s: shifting by the measured delay %1.3f us (%d coded pulses)"
              % (channel, delays[channel], len(ds)))
    return delays


def plot_ipp_power_closeup(zs, t_us, sweepid, ipp_idx, t_unix, i0, i1, delays, ofname):
    """The transmit portion of the combined power plot, with the fitted delays."""
    fig, ax = new_ipp_figure(sweepid, " + ".join(zs.keys()), ipp_idx, t_unix, "transmit pulse power")

    for ci, channel in enumerate(zs.keys()):
        ax.plot(t_us[channel][i0:i1], power_db(zs[channel][i0:i1]), color="C%d" % (ci % 10),
                lw=0.7, alpha=0.8, label=channel)

    if len(delays) > 0:
        ax.set_title("%s\n%s" % (ax.get_title(),
                                 "  ".join(["%s delay %1.3f $\\mu$s" % (ch, delays[ch]) for ch in delays])),
                     fontsize=11)

    ax.set_xlim(t_us[list(zs.keys())[0]][i0], t_us[list(zs.keys())[0]][i1 - 1])
    ax.set_ylabel("Power (dB, uncalibrated)")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(ofname, dpi=150)
    plt.close(fig)


def plot_tx_delay(ipp_idxs, delays, ofname, residual=False):
    """Matched filter delay of each echo channel relative to the transmit channel."""
    fig, ax = plt.subplots(figsize=(10, 5))

    for ci, channel in enumerate(delays.keys()):
        d = n.array(delays[channel])
        ax.plot(ipp_idxs, d, ".-", color="C%d" % (ci % 10), lw=0.8, ms=8,
                label="%s (mean %1.3f, std %1.3f $\\mu$s)" % (channel, n.nanmean(d), n.nanstd(d)))

    ax.set_xlabel("IPP number")
    ax.set_ylabel(r"%s relative to %s ($\mu$s)" % ("Residual delay" if residual else "Delay", TX_CHANNEL))
    ax.set_title("%s of the transmit pulse, %d IPPs"
                 % ("Residual matched filter delay after correction" if residual
                    else "Matched filter delay", len(ipp_idxs)), fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(ofname, dpi=150)
    plt.close(fig)


def plot_ipp_power_combined(zs, t_us, sweepid, ipp_idx, t_unix, ofname):
    """One interpulse period: power in dB for every channel in one plot.

    zs maps channel name -> voltage samples, t_us maps channel name -> time axis.
    All channels are read from the same sample index, so the time axes line up.
    """
    fig, ax = new_ipp_figure(sweepid, " + ".join(zs.keys()), ipp_idx, t_unix, "power")

    t_max = 0.0
    for ci, channel in enumerate(zs.keys()):
        ax.plot(t_us[channel], power_db(zs[channel]), color="C%d" % (ci % 10),
                lw=0.5, alpha=0.8, label=channel)
        t_max = max(t_max, t_us[channel][-1])

    ax.set_xlim(0, t_max)
    ax.set_ylabel("Power (dB, uncalibrated)")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(ofname, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Plot raw voltage, one figure per interpulse period.")
    ap.add_argument("-c", "--channel", action="append", default=None,
                    help="rf_data channel, repeatable (default: all channels in the dataset)")
    ap.add_argument("-n", "--n-ipp", type=int, default=10, help="number of interpulse periods to plot (default 10)")
    ap.add_argument("-s", "--start-ipp", type=int, default=0, help="skip this many IPPs before plotting (default 0)")
    ap.add_argument("--t0", type=float, default=None, help="unix time to start at (default: start of the recording)")
    ap.add_argument("--remove-dc", action="store_true", help="subtract the known USRP DC offset")
    ap.add_argument("--fix-tx-delay", action="store_true",
                    help="shift each echo channel onto the tx-h time base by its measured "
                         "channel delay, and write to %s instead of %s" % (OUT_SUBDIR_FIXED, OUT_SUBDIR))
    ap.add_argument("--oversample", type=int, default=100,
                    help="cross correlation interpolation factor for the delay estimate (default 100)")
    ap.add_argument("--max-lag-us", type=float, default=100.0,
                    help="only search for the matched filter peak within this many us (default 100)")
    ap.add_argument("-o", "--outdir", default=None,
                    help="output directory (default: <output_dir>/IPP from the config)")
    args = ap.parse_args()

    global TITLE_NOTE

    dirname, outbase = load_config()
    outdir = os.path.join(outbase, OUT_SUBDIR_FIXED if args.fix_tx_delay else OUT_SUBDIR)
    if args.outdir is not None:
        outdir = args.outdir
    print("config %s\ndata   %s\nplots  %s" % (CONFIG, dirname, outdir))

    d_il = DigitalRFReader("%s/rf_data/" % (dirname))
    available = d_il.get_channels()
    channels = args.channel if args.channel is not None else available
    for ch in channels:
        if ch not in available:
            print("ERROR: channel %s not in %s" % (ch, str(available)))
            sys.exit(1)

    keys, sweepids, ipp_len = get_ipps(dirname, args.n_ipp, args.start_ipp, args.t0)

    z_dc = Z_DC if args.remove_dc else n.complex64(0.0)

    srs = {}
    t_us = {}
    vdirs = {}
    pdirs = {}
    for channel in channels:
        srs[channel] = float(d_il.get_properties(channel)["samples_per_second"])
        t_us[channel] = n.arange(ipp_len) / srs[channel] * 1e6

        vdirs[channel] = os.path.join(outdir, channel, "raw voltage")
        pdirs[channel] = os.path.join(outdir, channel, "power")
        os.makedirs(vdirs[channel], exist_ok=True)
        os.makedirs(pdirs[channel], exist_ok=True)

        print("%s: %d IPPs of %d samples (%1.3f ms) starting %s -> %s"
              % (channel, len(keys), ipp_len, ipp_len / srs[channel] * 1e3,
                 stuffr.unix2datestr(float(keys[0]) / srs[channel]), os.path.join(outdir, channel)))

    # power of every channel in one plot, only makes sense with more than one
    cdir = None
    udir = None
    if len(channels) > 1:
        cdir = os.path.join(outdir, "combined")
        udir = os.path.join(cdir, "close_up")
        os.makedirs(cdir, exist_ok=True)
        os.makedirs(udir, exist_ok=True)
        print("combined power -> %s" % (cdir))

    # put every channel on the tx-h time base before plotting
    delay_fix = {}
    if args.fix_tx_delay:
        if TX_CHANNEL not in channels:
            print("ERROR: --fix-tx-delay needs the %s channel" % (TX_CHANNEL))
            sys.exit(1)
        delay_fix = measure_channel_delays(d_il, dirname, keys[0], channels, srs, ipp_len,
                                           oversample=args.oversample, max_lag_us=args.max_lag_us)
        TITLE_NOTE = "\necho channels shifted onto the %s time base" % (TX_CHANNEL)

    # matched filter delay of every other channel against the transmit channel
    echo_channels = [ch for ch in channels if ch != TX_CHANNEL] if TX_CHANNEL in channels else []
    delays = {ch: [] for ch in echo_channels}
    ipp_idxs = []

    for i, key in enumerate(keys):
        ipp_idx = args.start_ipp + i
        zs = {}

        for channel in channels:
            z = d_il.read_vector_1d(int(key), ipp_len, channel).astype("c8", casting="unsafe", copy=False) - z_dc
            if channel in delay_fix:
                # negative: the echo channel lags tx-h, so move it earlier
                z = fractional_shift(z, -delay_fix[channel] * srs[channel] / 1e6)
            zs[channel] = z
            t_unix = float(key) / srs[channel]
            fname = "ipp-%04d-%d.png" % (ipp_idx, int(t_unix))

            plot_ipp_voltage(z, t_us[channel], sweepids[i], channel, ipp_idx, t_unix,
                             os.path.join(vdirs[channel], fname))
            plot_ipp_power(z, t_us[channel], sweepids[i], channel, ipp_idx, t_unix,
                           os.path.join(pdirs[channel], fname))
            print("wrote %s and %s" % (os.path.join(vdirs[channel], fname),
                                       os.path.join(pdirs[channel], fname)))

        # matched filter the isolated transmit pulse of each echo channel
        # against the transmit channel
        ipp_idxs.append(ipp_idx)
        tw = tx_window(sweepids[i], ipp_len)
        ipp_delays = {}
        for channel in echo_channels:
            if tw is None:
                ipp_delays[channel] = n.nan
            else:
                d_samples, amp = matched_filter_delay(zs[TX_CHANNEL][tw[0]:tw[1]],
                                                      zs[channel][tw[0]:tw[1]],
                                                      oversample=args.oversample,
                                                      max_lag=args.max_lag_us * srs[channel] / 1e6)
                ipp_delays[channel] = d_samples / srs[channel] * 1e6
            delays[channel].append(ipp_delays[channel])

        if cdir is not None:
            t_unix = float(key) / srs[channels[0]]
            fname = "ipp-%04d-%d.png" % (ipp_idx, int(t_unix))
            plot_ipp_power_combined(zs, t_us, sweepids[i], ipp_idx, t_unix, os.path.join(cdir, fname))
            print("wrote %s" % (os.path.join(cdir, fname)))

            # the same thing, zoomed onto the transmit pulse
            uw = tx_window(sweepids[i], ipp_len, margin=CLOSEUP_MARGIN_US * srs[channels[0]] / 1e6)
            if uw is not None:
                plot_ipp_power_closeup(zs, t_us, sweepids[i], ipp_idx, t_unix, uw[0], uw[1],
                                       ipp_delays, os.path.join(udir, fname))
                print("wrote %s%s" % (os.path.join(udir, fname),
                                      "".join(["  %s delay %1.3f us" % (ch, ipp_delays[ch]) for ch in ipp_delays])))

    # delay as a function of interpulse period
    if len(echo_channels) > 0 and cdir is not None:
        fname = "tx_delay-%d.png" % (int(float(keys[0]) / srs[channels[0]]))
        plot_tx_delay(ipp_idxs, delays, os.path.join(cdir, fname), residual=args.fix_tx_delay)
        print("wrote %s" % (os.path.join(cdir, fname)))
        for channel in echo_channels:
            d = n.array(delays[channel])
            print("%s %s relative to %s: mean %1.4f us, std %1.4f us"
                  % (channel, "residual delay" if args.fix_tx_delay else "delay",
                     TX_CHANNEL, n.nanmean(d), n.nanstd(d)))


if __name__ == "__main__":
    main()
