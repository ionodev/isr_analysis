#!/usr/bin/env python3
"""
Injection test of the impulse-row removal (impulse_rows.py, Memo 47):
synthetic power-line spark impulses added to the raw voltage of the quiet
zenith cycle, where the truth is the same inversion without them.

The impulses have the properties Memo 21 measured on zenith-l:
    timing   locked to the 60 Hz mains on both half-cycles: each of
             N_GAPS spark gaps fires at its own fixed phase, every 1/120 s
             (the phases uniform at random)
    shape    one sample (80 %) or two consecutive samples (20 %), each with
             a random phase (a nanosecond spark is the receiver's impulse
             response at 1 MHz sampling)
    power    log-normal: median 53 times the noise power per sample, 1 %
             above 430 times (sigma of ln 0.90)
N_GAPS = 8 gives about 960 impulses per second, the rate of the driest spark
cycle of Memo 47 (8 April 15:03 UTC: 4.8e-3 of the echo-window samples
detected as impulses, about five detected samples per impulse with the
guard, and about a quarter of the impulses below the detector's threshold).
The impulses go into every sample the reader returns on zenith-l (echo,
noise and calibration windows alike, as real sparks do); tx-h is left as
recorded.

Four runs of the final product's inversion (variant 5, tracked columns at
P_FA 1e-6 with confirmation, outlier rejection off) on the quiet cycle:
clean_off, clean_on, inj_off, inj_on (impulses injected or not; impulse
rows removed or not). The effect of the injected impulses with the method
on or off is inj_X - clean_X: the existing impulses of the quiet cycle are
in both and cancel.

Usage:
    mpirun -np 9 python3 inject_impulses.py lpi <case> <out_dir>
    python3 inject_impulses.py measure <out_dir>
(the fits: run_analysis.py on a config with the lpi step disabled and
output_dir <out_dir>/<case>)
"""

import glob
import json
import os
import sys

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import numpy as n

import inject_satellites as inj

N_GAPS = 8
MAINS_HZ = 60.0
P_SINGLE = 0.8
POWER_MEDIAN = 53.0
POWER_SIGMA_LN = 0.90
SEED = 21
CASES = {"clean_off": (False, False), "clean_on": (False, True), "inj_off": (True, False), "inj_on": (True, True)}
QUIET_CATALOGUE = "/mnt/data/henrik/sat_v2_2026-10/quiet_catalogue_pfa1e-6"


def lpi_config(out, impulse_rows):
    return {"channel": inj.CHANNEL, "range_gate_us": 60, "avg_dur": 10, "min_tx_frac": 0.5,
            "pass_band_hz": 100000, "filter_len": 20, "max_range_delay_us": 7000, "save_acf_images": False,
            "lag_avg": 1, "reanalyze": False, "weighting": "variance", "background_column": True,
            "outlier_rejection": False, "periods": inj.PERIODS, "impulse_rows": impulse_rows,
            "satellites": {"catalogue": QUIET_CATALOGUE, "channel_delay_us": inj.CHANNEL_DELAY_US, "track": True,
                           "strong_pfa": 1e-10, "track_options": {"min_weak_pulses": 3}}}


class ImpulseReader:
    """A raw reader that adds the synthetic spark impulses to zenith-l."""

    def __init__(self, path, noise_power, seed=SEED):
        self.reader = inj.RawReader(path)
        self.noise_power = noise_power
        rng = n.random.default_rng(seed)
        # each gap's phase within the half-cycle, in microseconds
        self.phase_us = rng.uniform(0, 1e6 / (2 * MAINS_HZ), N_GAPS)
        self.seed = seed

    def __getattr__(self, name):
        return getattr(self.reader, name)

    def impulses(self, s0, n_samples):
        """(sample offsets, complex amplitudes) of the impulses in [s0, s0+n)."""
        half = 1e6 / (2 * MAINS_HZ)
        offs, amps = [], []
        for g, ph in enumerate(self.phase_us):
            k0 = int(n.ceil((s0 - ph) / half))
            k1 = int(n.floor((s0 + n_samples - 1 - ph) / half))
            for k in range(k0, k1 + 1):
                s = int(round(ph + k * half))
                # the impulse's own random numbers, the same whichever read
                # returns it
                r = n.random.default_rng([self.seed, g, k & 0x7FFFFFFF])
                width = 1 if r.uniform() < P_SINGLE else 2
                p = POWER_MEDIAN * n.exp(POWER_SIGMA_LN * r.standard_normal()) * self.noise_power
                for w in range(width):
                    if 0 <= s + w - s0 < n_samples:
                        offs.append(s + w - s0)
                        amps.append(n.sqrt(p) * n.exp(2j * n.pi * r.uniform()))
        return n.array(offs, dtype=int), n.array(amps, dtype=n.complex128)

    def read_vector_1d(self, key, n_samples, channel, *args, **kw):
        v = self.reader.read_vector_1d(key, n_samples, channel, *args, **kw).astype(n.complex64, casting="unsafe")
        if channel != inj.CHANNEL:
            return v
        o, a = self.impulses(int(key), int(n_samples))
        if len(o):
            v = v.astype(n.complex128)
            v[o] += a
            v = v.astype(n.complex64)
        return v


def cmd_lpi(case, out_dir):
    import outlier_lpi as olpi
    import run_analysis as ra
    inject, rows = CASES[case]
    out = os.path.join(out_dir, case)
    os.makedirs(out, exist_ok=True)
    if inject:
        # the noise power per sample: the median over the interval's
        # periods, measured on the clean data as in inject_satellites
        from raw_reader import RawReader
        reader = RawReader("%s/rf_data/" % inj.DATA)
        p = []
        for ai in inj.PERIODS[::6]:
            sid = inj.pulse_codes(inj.period_start(ai))
            keys = sorted(k for k in sid if sid[k] in inj.TMM and inj.TMM[sid[k]]["read_length"] <= 10000)
            p.append(inj.noise_power(reader, keys, sid))
        noise = float(n.median(p))
        olpi.RawReader = lambda path: ImpulseReader(path, noise)
        json.dump(dict(case=case, noise_power=noise, n_gaps=N_GAPS, mains_hz=MAINS_HZ, p_single=P_SINGLE,
                       power_median=POWER_MEDIAN, power_sigma_ln=POWER_SIGMA_LN, seed=SEED),
                  open(os.path.join(out, "injection.json"), "w"), indent=1)
    ra.run_lpi(inj.DATA, out, lpi_config(out, rows))


def fit_config(out_dir, case):
    return {"experiment": "inject_impulses_" + case, "data_dir": inj.DATA, "output_dir": os.path.join(out_dir, case),
            "steps": {"lpi": dict(lpi_config(out_dir, CASES[case][1]), enabled=False),
                      "fit_lpi": {"enabled": True, "channel": inj.CHANNEL, "max_dt": 60, "first_lag": 0,
                                  "range_avg": [1, 3, 5], "reanalyze": False}}}


def _lpi(d):
    out = {}
    for f in sorted(glob.glob(os.path.join(d, "lpi_60", "zenith-l", "lpi-*.h5"))):
        with h5py.File(f, "r") as h:
            out[os.path.basename(f)] = dict(acf=h["acfs_e"][()], var=h["acfs_var"][()], r=h["rgs_km"][()],
                                            removed=h["impulse_rows_removed_fraction"][()] if "impulse_rows_removed_fraction" in h else None,
                                            samples=float(h["impulse_sample_fraction"][()]) if "impulse_sample_fraction" in h else n.nan)
    return out


def cmd_measure(out_dir):
    import fit_quality as fq
    L = {c: _lpi(os.path.join(out_dir, c)) for c in CASES}
    files = sorted(set.intersection(*[set(v) for v in L.values()]))
    r = L["clean_off"][files[0]]["r"]
    band = (r >= 150) & (r <= 600)
    res = dict(periods=len(files))
    lags = None
    for m in ("off", "on"):
        z, sr = [], []
        for f in files:
            a, b = L["inj_" + m][f], L["clean_" + m][f]
            s = n.sqrt(b["var"] / 2)
            z.append(((a["acf"] - b["acf"]).real / s)[band])
            sr.append(n.sqrt(a["var"] / b["var"])[band])
        z, sr = n.array(z), n.array(sr)
        nl = z.shape[2]
        lags = n.arange(nl) * 10
        # mean over the periods and the gates, per lag, with its standard error
        zm = n.nanmean(z, axis=(0, 1))
        ze = n.nanstd(n.nanmean(z, axis=1), axis=0) / n.sqrt(len(files))
        res[m] = dict(shift_by_lag=zm.tolist(), shift_err=ze.tolist(),
                      shift_lags_1_6=float(n.nanmean(zm[1:7])), shift_lags_1_6_err=float(n.sqrt(n.nanmean(ze[1:7] ** 2) / 6)),
                      shift_lag0=float(zm[0]), sigma_ratio_by_lag=n.nanmedian(sr, axis=(0, 1)).tolist(),
                      sigma_ratio_lag1=float(n.nanmedian(sr[:, :, 1])))
    for c in ("clean_on", "inj_on"):
        res[c + "_rows_removed"] = float(n.nanmean([L[c][f]["removed"] for f in files]))
        res[c + "_impulse_sample_fraction"] = float(n.nanmean([L[c][f]["samples"] for f in files]))
    # the 60 s fits: inj minus clean, in the clean fit's standard deviation
    fits = {}
    for c in CASES:
        F = fq.load_fits(os.path.join(out_dir, c))
        if F is not None:
            fits[c] = F
    if len(fits) == 4:
        rf = fits["clean_off"]["r"]
        for m in ("off", "on"):
            A, B = fits["clean_" + m], fits["inj_" + m]
            out = {}
            for lo, hi in ((150, 250), (250, 450), (450, 600)):
                g = (rf >= lo) & (rf < hi)
                with n.errstate(invalid="ignore", divide="ignore"):
                    out["%d-%d" % (lo, hi)] = dict(
                        ne_ratio=float(n.nanmedian((B["ne"] / A["ne"])[:, g])),
                        ne_z=float(n.nanmedian(((B["ne"] - A["ne"]) / A["dne"])[:, g])),
                        Te_z=float(n.nanmedian(((B["Te"] - A["Te"]) / A["dTe"])[:, g])),
                        Ti_z=float(n.nanmedian(((B["Ti"] - A["Ti"]) / A["dTi"])[:, g])),
                        dTe=float(n.nanmedian((B["Te"] - A["Te"])[:, g])),
                        dTi=float(n.nanmedian((B["Ti"] - A["Ti"])[:, g])))
            res[m]["fits"] = out
    json.dump(res, open(os.path.join(out_dir, "measures.json"), "w"), indent=1)
    for m in ("off", "on"):
        x = res[m]
        print("removal %-3s: ACF shift lag 0 %+.3f sigma, lags 10-60 us %+.3f +- %.3f sigma; sigma ratio lag 1 %.3f" % (
            m, x["shift_lag0"], x["shift_lags_1_6"], x["shift_lags_1_6_err"], x["sigma_ratio_lag1"]))
        for k, v in x.get("fits", {}).items():
            print("   %s km: ne ratio %.3f (%+.2f sigma), Te %+.0f K (%+.2f sigma), Ti %+.0f K (%+.2f sigma)" % (
                k, v["ne_ratio"], v["ne_z"], v["dTe"], v["Te_z"], v["dTi"], v["Ti_z"]))
    print("rows removed: clean %.3f, injected %.3f; impulse sample fraction clean %.2e, injected %.2e" % (
        res["clean_on_rows_removed"], res["inj_on_rows_removed"],
        res["clean_on_impulse_sample_fraction"], res["inj_on_impulse_sample_fraction"]))


def main():
    if sys.argv[1] == "lpi":
        cmd_lpi(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "fitconfig":
        json.dump(fit_config(sys.argv[3], sys.argv[2]), open(sys.argv[4], "w"), indent=1)
    elif sys.argv[1] == "measure":
        cmd_measure(sys.argv[2])


if __name__ == "__main__":
    main()
