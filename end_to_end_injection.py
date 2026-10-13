#!/usr/bin/env python3
"""
End-to-end injection test of the satellite detection and removal.

Synthetic satellite passes are added to the raw voltage of the quiet zenith
cycle (as in inject_satellites.py), one pass per 10-s period. Unlike
inject_satellites.py, the passes have a light curve and the satellite
columns get only what the real detector finds:

    detect   the single-pulse CFAR detector (single_pulse_satellite_cfar.
             detect_pulse, unchanged) runs on every pulse of the injected
             data; an injected echo counts as detected when a detection on
             its pulse lies within MATCH_DELAY samples and MATCH_DOPPLER_HZ
             of it. Writes the detection probability against per-pulse SNR.
    lpi CASE the lag profile inversion of the injected data, variant 5
             (variance weights, background column) with satellite columns
             for the catalogue's real echoes plus, for the injected ones:
               perfect    every injected pulse, at its exact delay/Doppler
               realistic  the injected pulses the detector found, at the
                          detector's delay and Doppler
               realistic_refined  the same, with each delay refined below
                          the detector's 8-sample grid
                          (satellite_columns.refine_delay; the real
                          catalogue echoes stay unrefined, as in the clean
                          run)
               realistic_tracked  the detected ones linked into passes, with
                          columns for every pulse of a pass
                          (satellite_columns.TrackedDetections)
               none       none (variant 3: no satellite columns at all)
             Extra keywords for lpi_files (e.g. outlier_rejection=True)
             can follow as key=value and are added to the case's name.
    measure  the shift of the ACF at the echo's gate against the clean
             runs of the same variant, by the pass's peak SNR

The truth is the same variant on the same data without the injection, as in
Memo 42. The pass: 450 km at its centre, receding at -500 Hz (about 170 m/s),
in the beam from 3 to 7 s of the period, the matched-filter SNR falling by
12 dB from the centre to the edges (a Gaussian beam crossing), peak SNRs
8-23 dB in steps of 3 dB over the periods.

Usage:
    python3 end_to_end_injection.py detect  <out_dir>
    mpirun -np 6 python3 end_to_end_injection.py lpi <case> <out_dir> [key=value ...]
    python3 end_to_end_injection.py measure <out_dir> <clean_variants_dir>
"""

import argparse
import glob
import json
import os
import sys
from multiprocessing import Pool

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import numpy as n
import scipy.constants as sc
from digital_rf import DigitalRFReader, DigitalMetadataReader

import inject_satellites as inj
from radar_timing import TMM
from tx_delay import fractional_shift

PEAK_SNR_DB = [8, 11, 14, 17, 20, 23]
WINDOW_S = (3.0, 7.0)
EDGE_DROP_DB = 12.0
RANGE_KM = 450.03
DOPPLER_HZ = -500.0
MATCH_DELAY = 8.0
MATCH_DOPPLER_HZ = 1500.0


def plan_passes(seed=2):
    """{pulse key: echo} for one pass per period, with its light curve."""
    from raw_reader import RawReader
    reader = RawReader("%s/rf_data/" % inj.DATA)
    rng = n.random.default_rng(seed)
    rdot = -DOPPLER_HZ * inj.WAVELENGTH / 2.0
    d0 = RANGE_KM * 1e3 * 2 / sc.c * 1e6
    plan, sid_all, meta = {}, {}, []
    half = 0.5 * (WINDOW_S[1] - WINDOW_S[0])
    for j, ai in enumerate(inj.PERIODS):
        peak = PEAK_SNR_DB[j % len(PEAK_SNR_DB)]
        i0 = inj.period_start(ai)
        sid = inj.pulse_codes(i0)
        sid_all.update(sid)
        keys = sorted(k for k in sid if sid[k] in TMM and TMM[sid[k]]["read_length"] <= 10000)
        before = [k for k in keys if (k - i0) / 1e6 < WINDOW_S[0] - 0.2]
        p_noise = inj.noise_power(reader, before, sid)
        tc = i0 + 1e6 * 0.5 * (WINDOW_S[0] + WINDOW_S[1])
        for k in keys:
            ts = (k - i0) / 1e6
            if WINDOW_S[0] <= ts < WINDOW_S[1]:
                x = (ts - 0.5 * (WINDOW_S[0] + WINDOW_S[1])) / half
                snr = peak - EDGE_DROP_DB * x ** 2
                plan[k] = dict(delay=d0 + 2 * rdot / sc.c * 1e6 * (k - tc) / 1e6, doppler=DOPPLER_HZ,
                               p_true=10 ** (snr / 10) * p_noise, snr_db=snr, peak_db=peak, period=ai,
                               phase=float(rng.uniform(0, 2 * n.pi)))
        meta.append(dict(period=ai, peak_snr_db=peak, noise_power=p_noise))
    return plan, sid_all, meta


def add_echo(v, key, plan, sid, n_samples):
    p = plan[key]
    tm = TMM[sid[key]]
    u = n.zeros(n_samples, dtype=n.complex128)
    u[tm["tx0"]:tm["tx1"]] = v[tm["tx0"]:tm["tx1"]] - inj.Z_DC
    amp = n.sqrt(p["p_true"] / n.sum(n.abs(u) ** 2))
    t = n.arange(n_samples)
    return v + amp * n.exp(1j * p["phase"]) * fractional_shift(u, p["delay"]) \
        * n.exp(2j * n.pi * p["doppler"] * t / 1e6)


class DetectorReader:
    """DigitalRFReader.read_vector for the detector, with the injection on
    the echo channel (tx-h is left as recorded)."""

    def __init__(self, plan, sid):
        self.reader = DigitalRFReader("%s/rf_data/" % inj.DATA)
        self.plan, self.sid = plan, sid

    def read_vector(self, sample, n_samples, channel, *a, **kw):
        v = self.reader.read_vector(sample, n_samples, channel, *a, **kw).astype(n.complex64)
        if channel != inj.CHANNEL or sample not in self.plan:
            return v
        if n_samples < TMM[self.sid[sample]]["tx1"]:
            return v
        return add_echo(v.astype(n.complex128), sample, self.plan, self.sid, n_samples).astype(n.complex64)


def _detect_period(args):
    ai, plan, sid = args
    import single_pulse_satellite_cfar as cfar
    dargs = cfar.parser().parse_args(["--data", inj.DATA, "--output", "/dev/null"])
    bank = cfar.FFTBank(dargs.fft_batch)
    rf = DetectorReader(plan, sid)
    out = []
    for k in sorted(x for x in plan if plan[x]["period"] == ai):
        if sid[k] not in cfar.MODES:
            continue
        try:
            dets = cfar.detect_pulse(rf, int(k), int(sid[k]), inj.CHANNEL, dargs, bank)
        except Exception as e:
            print("pulse %d: %s" % (k, e), file=sys.stderr)
            continue
        p = plan[k]
        best = None
        for d in dets:
            delay = d["raw_delay_sample"] - TMM[sid[k]]["tx0"] - inj.CHANNEL_DELAY_US
            if abs(delay - p["delay"]) <= MATCH_DELAY and abs(d["doppler_hz"] - p["doppler"]) <= MATCH_DOPPLER_HZ:
                if best is None or d["snr_db"] > best[2]:
                    best = (delay, d["doppler_hz"], d["snr_db"])
        out.append((int(k), p["snr_db"], p["peak_db"], best))
    return out


class PlanDetections:
    """The catalogue's real echoes plus the injected ones of a case."""

    def __init__(self, extra, tracker=None):
        import satellite_columns as satcol
        self.extra = extra
        self.tracker = tracker
        self.catalogue = satcol.CatalogueDetections(inj.CATALOGUE, inj.CHANNEL, inj.CHANNEL_DELAY_US)

    def for_period(self, i0, i1):
        out = self.catalogue.for_period(i0, i1)
        if self.tracker is not None:
            for k, v in self.tracker.for_period(i0, i1).items():
                out.setdefault(k, []).extend(v)
            print("tracked passes:", self.tracker.last_passes, flush=True)
            return out
        for k, v in self.extra.items():
            if i0 <= k < i1:
                out.setdefault(k, []).append(v)
        return out


class _Rows:
    def __init__(self, rows):
        self.r = rows

    def rows(self, i0, i1):
        return [x for x in self.r if i0 <= x[0] < i1]


def cmd_detect(out_dir):
    plan, sid, meta = plan_passes()
    with Pool(len(inj.PERIODS)) as pool:
        res = [r for rr in pool.map(_detect_period, [(ai, plan, sid) for ai in inj.PERIODS]) for r in rr]
    os.makedirs(out_dir, exist_ok=True)
    snr = n.array([r[1] for r in res])
    found = n.array([r[3] is not None for r in res])
    bins = n.arange(-6, 27, 3)
    curve = []
    for a, b in zip(bins[:-1], bins[1:]):
        m = (snr >= a) & (snr < b)
        if m.any():
            curve.append(dict(snr_lo=float(a), snr_hi=float(b), pulses=int(m.sum()), detected=float(found[m].mean())))
    detected = {str(r[0]): [float(r[3][0]), float(r[3][1]), float(r[3][2])] for r in res if r[3] is not None}
    json.dump(dict(passes=meta, curve=curve, pulses=len(res), detected_fraction=float(found.mean()),
                   detected=detected), open(os.path.join(out_dir, "detections.json"), "w"), indent=1)
    for c in curve:
        print("SNR %5.1f-%5.1f dB: %4d pulses, detected %.2f" % (c["snr_lo"], c["snr_hi"], c["pulses"], c["detected"]))


def cmd_lpi(case, out_dir, extra_kw):
    import outlier_lpi as olpi
    plan, sid, meta = plan_passes()
    kw = dict(inj.LPI)
    name = case
    for x in extra_kw:
        k, v = x.split("=")
        kw[k] = json.loads(v.lower()) if v.lower() in ("true", "false") else json.loads(v)
        name += "_%s-%s" % (k, v)
    if case == "perfect":
        dets = PlanDetections({k: (p["delay"], p["doppler"]) for k, p in plan.items()})
    elif case == "realistic":
        found = json.load(open(os.path.join(out_dir, "detections.json")))["detected"]
        dets = PlanDetections({int(k): (v[0], v[1]) for k, v in found.items()})
    elif case == "realistic_refined":
        import satellite_columns as satcol
        found = json.load(open(os.path.join(out_dir, "detections.json")))["detected"]
        rf = DetectorReader(plan, sid)
        extra = {}
        for k, v in found.items():
            k = int(k)
            tm = TMM[sid[k]]
            echo = rf.read_vector(k, 10000, inj.CHANNEL)
            tx = rf.read_vector(k, 10000, "tx-h")
            raw = v[0] + tm["tx0"] + inj.CHANNEL_DELAY_US
            raw_ref, _ = satcol.refine_delay(echo, tx, raw, v[1], tm["tx0"], tm["tx1"])
            extra[k] = (raw_ref - tm["tx0"] - inj.CHANNEL_DELAY_US, v[1])
        dets = PlanDetections(extra)
    elif case.startswith("realistic_tracked"):
        # realistic_tracked, optionally _ext<seconds> (the default extension
        # of passes too short to fit a light curve) and _min<n> (the fewest
        # detections that make a pass)
        import satellite_columns as satcol
        found = json.load(open(os.path.join(out_dir, "detections.json")))["detected"]
        rows = [(int(k), v[0], v[1], v[2]) for k, v in found.items()]
        tk = {}
        for part in case.split("_")[2:]:
            if part.startswith("ext"):
                tk["default_extend_s"] = float(part[3:])
            elif part.startswith("min"):
                tk["min_detections"] = int(part[3:])
        dets = PlanDetections({}, satcol.TrackedDetections(_Rows(rows), satcol.zenith_pulse_keys(inj.DATA), **tk))
    elif case == "none":
        dets = None
    else:
        raise ValueError(case)
    out = os.path.join(out_dir, name)
    os.makedirs(out, exist_ok=True)
    json.dump(dict(case=case, name=name, kw={k: v for k, v in kw.items() if k not in inj.LPI or kw[k] != inj.LPI[k]},
                   passes=meta), open(os.path.join(out, "case.json"), "w"), indent=1)
    olpi.RawReader = lambda path: inj.InjectingReader(path, plan, sid)
    olpi.lpi_files(dirname=inj.DATA, channel=inj.CHANNEL, output_base=out, periods=inj.PERIODS,
                   weighting="variance", background_column=True, satellite_detections=dets, **kw)


def cmd_measure(out_dir, clean_root):
    rows = []
    for d in sorted(glob.glob(os.path.join(out_dir, "*", "case.json"))):
        case = json.load(open(d))
        cd = os.path.dirname(d)
        clean = os.path.join(clean_root, "v3_variance_background" if case["case"] == "none"
                             else "v5_variance_background_satellites")
        by_peak = {}
        for p in case["passes"]:
            f = "lpi-%d.h5" % (1712484000 + 10 * p["period"])
            fi, fc = os.path.join(cd, "lpi_60", "zenith-l", f), os.path.join(clean, "lpi_60", "zenith-l", f)
            if not (os.path.exists(fi) and os.path.exists(fc)):
                continue
            with h5py.File(fi, "r") as hi, h5py.File(fc, "r") as hc:
                r = hc["rgs_km"][()]
                g = int(n.argmin(n.abs(r - RANGE_KM)))
                z = ((hi["acfs_e"][g, 1:6] - hc["acfs_e"][g, 1:6]).real / n.sqrt(hc["acfs_var"][g, 1:6] / 2)).mean()
            by_peak.setdefault(p["peak_snr_db"], []).append(float(z))
        row = dict(name=case["name"], shift_by_peak={str(k): float(n.mean(v)) for k, v in sorted(by_peak.items())})
        rows.append(row)
        print("%-40s %s" % (row["name"], "  ".join("%sdB %+6.1f" % (k, v) for k, v in row["shift_by_peak"].items())))
    json.dump(rows, open(os.path.join(out_dir, "shift_by_peak.json"), "w"), indent=1)


def main():
    cmd = sys.argv[1]
    if cmd == "detect":
        cmd_detect(sys.argv[2])
    elif cmd == "lpi":
        cmd_lpi(sys.argv[2], sys.argv[3], sys.argv[4:])
    elif cmd == "measure":
        cmd_measure(sys.argv[2], sys.argv[3])


if __name__ == "__main__":
    main()
