#!/usr/bin/env python3
"""
Collect the system noise temperature time series that the analysis already
produced, and attach the antenna pointing to it.

Every integration period written by outlier_lpi.py and avg_range_doppler_spec.py
carries a noise injection calibration:

    alpha = (median(P_inj+noise) - median(P_noise)) / T_INJECTION
    T_sys = median(P_noise) / alpha

The background gate is the last 500 microseconds before the injection pulse,
which is beyond the ionosphere (about 1080 to 1155 km for the 8910 microsecond
interpulse period), so T_sys there is receiver plus spillover plus sky, with no
ionospheric echo in it. That makes it directly comparable to a sky model.

What the analysis does not record is where the antenna was looking, because the
plasma parameter fits do not need it. This module reads the pointing back out of
the antenna control metadata and joins it onto each integration period, so the
T_sys series can be compared against a sky brightness model. See
sky_noise_model.py for the model and compare_tsys_sky.py for the comparison.

The output is one HDF5 file with a row per integration period. Every output set
found is harvested and tagged, including the low elevation horizon scan, so that
selecting between them is a decision made when fitting rather than when reading.

    python3 tsys_harvest.py config/millstone_eclipse2024.json [-o tsys.h5]
"""

import numpy as n
import h5py
import glob
import os
import sys
import json
import argparse
from multiprocessing import Pool

import millstone_radar_state as mrs

# The zenith antenna does not move, and nothing in the recording states where it
# points the way the antenna control metadata states MISA's. It is not at the
# geodetic zenith: the hard target work measured the beam at elevation 88.16,
# azimuth 172.9, a 1.84 degree displacement towards the south, from 251 tight
# tolerance satellite matches with Rayleigh Z = 367, using misa-l and its
# recorded encoder pointing as a control that showed no offset.
#
# The 440 MHz data confirm it independently. A tilt of 1.84 degrees towards
# azimuth 172.9 puts the boresight at declination
#
#     42.619 + 1.84 * cos(172.9) = 40.79 degrees
#
# and Cygnus A, the second brightest radio source in the sky, sits at 40.734.
# The two agree to 0.06 degrees, a twelfth of a beamwidth. Cyg A accordingly
# transits through the centre of this beam once per sidereal day and drives
# T_sys from 170 K to 1600 K, which is impossible for a beam pointed at the
# geodetic zenith: it would put the source 1.89 degrees off boresight, some 2.7
# beamwidths out, where even a generous sidelobe envelope leaves a few kelvin.
#
# Whether the displacement is mechanical, a survey error or a convention in the
# metadata is still unresolved; it is measured here, not explained.
ZENITH_EL_DEG = 88.16
ZENITH_AZ_DEG = 172.9

# Which receiving antenna each recorded channel belongs to. The noise gate of a
# channel measures that channel's own antenna regardless of which antenna was
# transmitting, so this mapping, and not the tx_antenna state, is what selects
# the beam to model.
CHANNEL_ANTENNA = {"zenith-l": "zenith", "misa-l": "misa"}

# Scalars to copy out of each file when they are present. T_sys and alpha are
# the measurement; the rest are kept because they are the obvious covariates to
# check a suspicious T_sys against.
# pass_band and filter_len are written only by outlier_lpi.py. They are needed
# to undo the bias its noise estimate carries, which depends on how many
# independent samples the filtered background window holds; see
# compare_tsys_sky.py. avg_range_doppler_spec.py does not record its passband,
# which is the module constant pass_band = 0.05e6 there.
SCALARS = ["T_sys", "alpha", "P_tx", "n_pulses", "mode", "pass_band", "filter_len"]


def period_files(path):
    """
    The per integration period files of one output directory.

    The fitting steps write their pp-*.h5 results into the same directory. Those
    carry a T_sys too, but it is an average over the periods that went into the
    fit and it has no i0, so it is neither a new measurement nor placeable in
    time. Only the per period files are harvested.
    """
    return sorted(f for f in glob.glob(os.path.join(path, "*.h5"))
                  if not os.path.basename(f).startswith("pp"))


def read_one(fname):
    """One integration period, or None if the file carries no T_sys."""
    try:
        with h5py.File(fname, "r") as h:
            if "T_sys" not in h.keys():
                return None
            row = {}
            for k in SCALARS:
                row[k] = float(h[k][()]) if k in h.keys() else n.nan

            # avg_range_doppler_spec.py writes i0 in microseconds, outlier_lpi.py
            # writes it in seconds. Both are unix epoch, so the magnitude tells
            # them apart without having to know which script wrote the file.
            i0 = float(h["i0"][()])
            t0 = i0 / 1e6 if i0 > 1e12 else i0
            row["t0"] = t0

            if "i1" in h.keys():
                i1 = float(h["i1"][()])
                row["t1"] = i1 / 1e6 if i1 > 1e12 else i1
            else:
                row["t1"] = n.nan
        return row
    except Exception as e:
        print("could not read %s: %s" % (fname, e))
        return None


def harvest_set(datadir, subdir, channel, pool):
    """Every output file of one analysis set for one channel."""
    path = os.path.join(datadir, subdir, channel)
    fl = period_files(path)
    if len(fl) == 0:
        return []
    rows = [r for r in pool.map(read_one, fl, chunksize=64) if r is not None]
    for r in rows:
        r["set"] = subdir
        r["channel"] = channel
    print("%-34s %-9s %6d periods" % (subdir, channel, len(rows)))
    return rows


def find_sets(datadir):
    """Analysis output directories, as (subdir, channel) pairs."""
    out = []
    for sub in sorted(os.listdir(datadir)):
        if not (sub.startswith("range_doppler") or sub.startswith("lpi_")):
            continue
        d = os.path.join(datadir, sub)
        if not os.path.isdir(d):
            continue
        for channel in sorted(os.listdir(d)):
            if channel not in CHANNEL_ANTENNA:
                continue
            if len(period_files(os.path.join(d, channel))) > 0:
                out.append((sub, channel))
    return out


def attach_pointing(rows, datadir):
    """
    Azimuth and elevation of the receiving antenna at the middle of each
    integration period, and the transmit and receive antenna selection.

    MISA pointing comes from the antenna control metadata. The zenith antenna
    does not move, so its pointing is a constant, but not the geodetic zenith;
    see ZENITH_EL_DEG. The pointing is interpolated at the
    centre of the period rather than its start, because the horizon scan moves
    several beamwidths in the ten seconds a period lasts.
    """
    meta = os.path.join(datadir, "metadata", "antenna_control_metadata")
    azf, elf, bounds = mrs.get_misa_az_el_model(meta)
    tx_sel, rx_sel = mrs.get_antenna_select(meta)

    t0 = n.array([r["t0"] for r in rows])
    t1 = n.array([r["t1"] for r in rows])
    # a period with no recorded end is one of the fixed length ones; the pointing
    # is then taken at its start, which for the zenith antenna changes nothing
    # and for MISA is half an integration period early.
    tmid = n.where(n.isfinite(t1), 0.5 * (t0 + t1), t0)

    # outside the metadata bounds the interpolation would raise; those periods
    # get no pointing rather than an extrapolated one.
    inside = (tmid >= bounds[0]) & (tmid <= bounds[1])
    az = n.full(len(rows), n.nan)
    el = n.full(len(rows), n.nan)
    az[inside] = azf(tmid[inside])
    el[inside] = elf(tmid[inside])

    # the antenna selection holds its first and last value outside the recorded
    # range, so it needs no such guard
    tx_ant = tx_sel(tmid * 1e6)
    rx_ant = rx_sel(tmid * 1e6)

    for i, r in enumerate(rows):
        if CHANNEL_ANTENNA[r["channel"]] == "zenith":
            r["az"] = ZENITH_AZ_DEG
            r["el"] = ZENITH_EL_DEG
        else:
            r["az"] = az[i]
            r["el"] = el[i]
        r["tmid"] = tmid[i]
        # 1 = MISA, -1 = zenith, 0 = unknown (around an antenna change), as
        # millstone_radar_state defines it
        r["tx_antenna"] = tx_ant[i]
        r["rx_antenna"] = rx_ant[i]
    return rows


def write(rows, outname):
    rows = sorted(rows, key=lambda r: (r["set"], r["channel"], r["t0"]))
    with h5py.File(outname, "w") as ho:
        for k in SCALARS + ["t0", "t1", "tmid", "az", "el", "tx_antenna", "rx_antenna"]:
            ho[k] = n.array([r[k] for r in rows])
        # variable length strings, so that a row identifies its own provenance
        st = h5py.string_dtype()
        ho.create_dataset("set", data=n.array([r["set"] for r in rows], dtype=object), dtype=st)
        ho.create_dataset("channel", data=n.array([r["channel"] for r in rows], dtype=object), dtype=st)
        ho.create_dataset("antenna",
                          data=n.array([CHANNEL_ANTENNA[r["channel"]] for r in rows], dtype=object),
                          dtype=st)
    print("wrote %d integration periods to %s" % (len(rows), outname))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config", help="analysis config json, read for data_dir")
    ap.add_argument("-o", "--output", default="tsys.h5")
    ap.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    args = ap.parse_args()

    with open(args.config, "r") as f:
        datadir = json.load(f)["data_dir"]
    print("reading %s" % (datadir))

    sets = find_sets(datadir)
    if len(sets) == 0:
        print("no analysis output found under %s" % (datadir))
        sys.exit(1)

    rows = []
    with Pool(args.jobs) as pool:
        for sub, channel in sets:
            rows += harvest_set(datadir, sub, channel, pool)

    rows = attach_pointing(rows, datadir)
    write(rows, args.output)


if __name__ == "__main__":
    main()
