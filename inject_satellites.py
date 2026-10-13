#!/usr/bin/env python3
"""
Injection test of the LPI variants: a synthetic satellite pass added to the
raw voltage of a quiet interval, where the truth is known.

Each LPI variant (weeks 7-9: equal weights, variance weights, with the
background column, the satellite columns, or both) is run on the injected
data. The same variant on the same periods without the injection (made
separately, with the same settings) is the truth, so the difference between
the two is what the echo does to that variant's plasma estimate: its bias,
and the change of its error bars.

The echo is built as in Memo 15: each pulse's own transmit pulse, as it leaks
into the echo channel, delayed by the echo's delay (with a fractional part),
Doppler shifted, and scaled so that every pulse has the same echo energy, a
chosen matched-filter SNR times the noise power per sample. It is added by
wrapping the raw reader, so the inversion itself carries no test code. The
satellite variants get the catalogue's real echoes plus the injected ones,
at their exact delay and Doppler: a perfect detection, the best case.

Usage (one variant and scenario per call; run with mpirun for speed):

    mpirun -np 6 python3 inject_satellites.py <variant> <scenario> <out_dir>

variant: v1_equal, v2_variance, v3_variance_background,
         v4_variance_satellites, v5_variance_background_satellites
scenario: bright (30 dB) or faint (15 dB); both at 450 km, -500 Hz, in the
          beam from 2 to 8 s of every period
"""

import json
import os
import sys

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import numpy as n
import scipy.constants as sc
from digital_rf import DigitalMetadataReader

import outlier_lpi as olpi
import satellite_columns as satcol
from radar_timing import TMM
from raw_reader import RawReader
from tx_delay import fractional_shift

DATA = "/mnt/data/juha/millstone_hill/isr/eclipse2024/usrp-rx0-r_20240407T100000_20240409T110000"
CATALOGUE = "/mnt/data/henrik/satellite_detections_2026-10-02"
CHANNEL = "zenith-l"
CHANNEL_DELAY_US = 10.62
# the quiet reference interval: one zenith cycle, 7 April 20:51:00-20:53:50 UTC
PERIODS = list(range(3906, 3924))
# the inversion settings of the variant runs
LPI = dict(rg=60, avg_dur=10, min_tx_frac=0.5, pass_band=1e5, filter_len=20,
           maximum_range_delay=7000, save_acf_images=False, lag_avg=1, reanalyze=True,
           outlier_rejection=False)
VARIANTS = {
    "v1_equal": dict(weighting="equal", background_column=False, sat=False),
    "v2_variance": dict(weighting="variance", background_column=False, sat=False),
    "v3_variance_background": dict(weighting="variance", background_column=True, sat=False),
    "v4_variance_satellites": dict(weighting="variance", background_column=False, sat=True),
    "v5_variance_background_satellites": dict(weighting="variance", background_column=True, sat=True),
}
SCENARIOS = {"bright": 30.0, "faint": 15.0}
# the synthetic satellite: away from the Starlink shells at 535-590 km
SAT_RANGE_KM = 450.03
SAT_DOPPLER_HZ = -500.0
SAT_WINDOW_S = (2.0, 8.0)
# the constant DC offset outlier_lpi subtracts from zenith-l
Z_DC = n.complex64(-0.212 - 0.221j)
WAVELENGTH = sc.c / 440.2e6


def period_start(ai):
    b0 = DigitalMetadataReader("%s/metadata/id_metadata" % DATA).get_bounds()[0]
    return b0 + ai * int(LPI["avg_dur"] * 1e6)


def pulse_codes(i0):
    r = DigitalMetadataReader("%s/metadata/id_metadata" % DATA)
    sid = r.read(i0, i0 + int(LPI["avg_dur"] * 1e6) + 40000, "sweepid")
    return {int(k): int(v) for k, v in sid.items()}


def noise_power(reader, keys, sid):
    """Noise power per raw sample, from the window before the noise injection,
    on the pulses before the satellite enters the beam."""
    pw = []
    for k in keys[:60]:
        le = TMM[sid[k]]["last_echo"]
        v = reader.read_vector_1d(k, 10000, CHANNEL).astype(n.complex64, casting="unsafe")
        pw.append(n.mean(n.abs(v[(le - 500):le] - Z_DC) ** 2))
    return float(n.median(pw))


def make_plan(snr_db, seed=1):
    """{pulse key: echo} for every pulse of the interval in the beam window."""
    reader = RawReader("%s/rf_data/" % DATA)
    plan, sid_all, meta = {}, {}, []
    rng = n.random.default_rng(seed)
    rdot = -SAT_DOPPLER_HZ * WAVELENGTH / 2.0
    d0 = SAT_RANGE_KM * 1e3 * 2 / sc.c * 1e6
    for ai in PERIODS:
        i0 = period_start(ai)
        sid = pulse_codes(i0)
        sid_all.update(sid)
        keys = sorted(k for k in sid if sid[k] in TMM and TMM[sid[k]]["read_length"] <= 10000)
        before = [k for k in keys if (k - i0) / 1e6 < SAT_WINDOW_S[0] - 0.2]
        p_noise = noise_power(reader, before, sid)
        p_true = 10 ** (snr_db / 10) * p_noise
        tc = i0 + 1e6 * 0.5 * (SAT_WINDOW_S[0] + SAT_WINDOW_S[1])
        for k in keys:
            ts = (k - i0) / 1e6
            if SAT_WINDOW_S[0] <= ts < SAT_WINDOW_S[1]:
                plan[k] = dict(delay=d0 + 2 * rdot / sc.c * 1e6 * (k - tc) / 1e6,
                               doppler=SAT_DOPPLER_HZ, p_true=p_true,
                               phase=float(rng.uniform(0, 2 * n.pi)))
        meta.append(dict(period=ai, noise_power=p_noise, echo_energy=p_true))
    return plan, sid_all, meta


class InjectingReader:
    """A raw reader that adds the planned echo to the planned pulses."""

    def __init__(self, path, plan, sid):
        self.reader = RawReader(path)
        self.plan = plan
        self.sid = sid

    def __getattr__(self, name):
        return getattr(self.reader, name)

    def read_vector_1d(self, key, n_samples, channel, *args, **kw):
        v = self.reader.read_vector_1d(key, n_samples, channel, *args, **kw).astype(n.complex64, casting="unsafe")
        if channel != CHANNEL or key not in self.plan:
            return v
        p = self.plan[key]
        tm = TMM[self.sid[key]]
        u = n.zeros(n_samples, dtype=n.complex128)
        u[tm["tx0"]:tm["tx1"]] = v[tm["tx0"]:tm["tx1"]] - Z_DC
        amp = n.sqrt(p["p_true"] / n.sum(n.abs(u) ** 2))
        t = n.arange(n_samples)
        echo = amp * n.exp(1j * p["phase"]) * fractional_shift(u, p["delay"]) \
            * n.exp(2j * n.pi * p["doppler"] * t / 1e6)
        return (v + echo).astype(n.complex64)


class Detections:
    """The catalogue's echoes plus the injected ones, as satellite_detections."""

    def __init__(self, plan):
        self.plan = plan
        self.catalogue = satcol.CatalogueDetections(CATALOGUE, CHANNEL, CHANNEL_DELAY_US)

    def for_period(self, i0, i1):
        out = self.catalogue.for_period(i0, i1)
        for k, p in self.plan.items():
            if i0 <= k < i1:
                out.setdefault(k, []).append((p["delay"], p["doppler"]))
        return out


def main():
    variant, scenario, out_dir = sys.argv[1:4]
    v = VARIANTS[variant]
    plan, sid, meta = make_plan(SCENARIOS[scenario])
    out = os.path.join(out_dir, scenario, variant)
    os.makedirs(out, exist_ok=True)
    json.dump(dict(variant=variant, scenario=scenario, snr_db=SCENARIOS[scenario],
                   range_km=SAT_RANGE_KM, doppler_hz=SAT_DOPPLER_HZ, window_s=SAT_WINDOW_S,
                   pulses_injected=len(plan), periods=meta),
              open(os.path.join(out, "injection.json"), "w"), indent=1)
    olpi.RawReader = lambda path: InjectingReader(path, plan, sid)
    olpi.lpi_files(dirname=DATA, channel=CHANNEL, output_base=out, periods=PERIODS,
                   weighting=v["weighting"], background_column=v["background_column"],
                   satellite_detections=Detections(plan) if v["sat"] else None, **LPI)


if __name__ == "__main__":
    main()
