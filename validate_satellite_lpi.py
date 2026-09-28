#!/usr/bin/env python3
"""
Test 1 of the satellite column method: a synthetic satellite echo injected into
the raw voltage of a clean integration period.

The echo is the transmit pulse leaking into the echo channel, delayed by a
fractional number of samples, Doppler shifted and scaled to a chosen matched
filter SNR. That is the same template the columns are built from, so this tests
the machinery (the handoff, the columns, the solve) with the truth known
exactly, not how well the leakage describes a real satellite's echo. The
injection is done by wrapping the DigitalRF reader, so the inversion itself
carries no test code.

The period is inverted several ways:

    clean          no injection, no outlier rejection  reference for the method
    clean_current  no injection, outlier rejection     the pipeline as it was
    nohandling     injected, no columns, no rejection  shows the bias
    fit            injected, columns, no rejection     the method
    current        injected, no columns, rejection     what the pipeline did
    fit_dop0       as fit, columns built with 0 Hz     range rate from the phase

    python3 validate_satellite_lpi.py [--snr-db 30] [--out plots/satfit_tests]
"""

import os
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
# one thread per process, as outlier_lpi.py intends. Set here because numpy is
# imported below before outlier_lpi is, and its thread pool would otherwise be
# sized to the whole machine in every one of the parallel runs.
for _v in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ[_v] = "1"
import argparse
import json

import numpy as n
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.constants as sc
from digital_rf import DigitalRFReader, DigitalMetadataReader

import outlier_lpi as olpi
import millstone_radar_state as mrs
from radar_timing import TMM
from tx_delay import fractional_shift

DATA = "/mnt/data/juha/millstone_hill/isr/eclipse2024/usrp-rx0-r_20240407T100000_20240409T110000"
CHANNEL = "zenith-l"
# inversion settings of config/millstone_eclipse2024.json
LPI = dict(avg_dur=10, rg=60, min_tx_frac=0.5, pass_band=1e5, filter_len=20,
           maximum_range_delay=7000, save_acf_images=True, lag_avg=1, reanalyze=True,
           # as Test 1 and the first Test 2 runs were made, before these became
           # the defaults; runs that differ say so
           precise_weights=False, satellite_template="pulse")
# the constant DC offset outlier_lpi subtracts from zenith-l
Z_DC = n.complex64(-0.212 - 0.221j)
WAVELENGTH = sc.c / 440.2e6
# outlier_lpi.lpi_files() skips pulses transmitted below this, its default
MIN_TX_PWR = 400e3

# the synthetic satellite: 550 km, receding at 170 m/s (-500 Hz), in the beam
# for one second in the middle of the period. The delay has a fractional part
# so the fractional shift is exercised.
SAT_RANGE_KM = 550.03
SAT_DOPPLER_HZ = -500.0
SAT_WINDOW_S = (4.5, 5.5)


_ANT = {}


def _zenith_chunk(args):
    """Is each period of a chunk on the zenith antenna throughout? Sampled five
    times per period."""
    b0, step, periods = args
    tx_ant, rx_ant, zpm = _ANT["tx"], _ANT["rx"], _ANT["zpm"]
    # the transmitter must be on as well: with it below the inversion's
    # min_tx_pwr every pulse is skipped, and the catalogue is empty for the
    # same reason, which would make such a period look clean
    return [all((tx_ant(k) <= -0.99) and (rx_ant(k) <= -0.99) and (zpm(k / 1e6) >= MIN_TX_PWR)
                for k in [b0 + ai * step + int(f * step) for f in (0.05, 0.3, 0.5, 0.7, 0.95)])
            for ai in periods]


def _catalogue_hour(args):
    """Number of catalogued echoes on this channel in each period of one hour."""
    b0, step, h0, h1 = args
    sat = DigitalMetadataReader("%s/metadata/satellite_detections" % DATA)
    count = {}
    for k, ch in sat.read(h0, h1, "channel").items():
        if ch in (CHANNEL, CHANNEL.encode()):
            ai = int((k - b0) // step)
            count[ai] = count.get(ai, 0) + 1
    return count


def survey_periods(cache, n_proc=None):
    """For every period of the recording: on the zenith antenna throughout, and
    how many catalogued echoes it has on this channel. Spread over all CPUs:
    the periods in chunks for the antenna state, the catalogue by the hour."""
    if os.path.exists(cache):
        z = n.load(cache)
        return z["zen"], z["cnt"], int(z["b0"])
    import multiprocessing
    n_proc = n_proc or multiprocessing.cpu_count()
    idr = DigitalMetadataReader("%s/metadata/id_metadata" % DATA)
    idb = idr.get_bounds()
    step = int(LPI["avg_dur"] * 1e6)
    n_per = int((idb[1] - idb[0]) // step)
    # built once here, inherited by the forked workers
    _ANT["tx"], _ANT["rx"] = mrs.get_antenna_select("%s/metadata/antenna_control_metadata" % DATA)
    _ANT["zpm"], _ = mrs.get_tx_power_model("%s/metadata/powermeter" % DATA)
    ctx = multiprocessing.get_context("fork")
    with ctx.Pool(n_proc) as pool:
        chunks = n.array_split(n.arange(n_per), n_proc * 4)
        zen = n.concatenate([n.array(z, dtype=bool) for z in
                             pool.map(_zenith_chunk, [(idb[0], step, list(c)) for c in chunks])])
        hours = [(idb[0], step, h, min(h + 3600 * 10**6, idb[1])) for h in range(idb[0], idb[1], 3600 * 10**6)]
        cnt = n.zeros(n_per, dtype=int)
        for c in pool.map(_catalogue_hour, hours):
            for ai, v in c.items():
                if 0 <= ai < n_per:
                    cnt[ai] += v
    n.savez(cache, zen=zen, cnt=cnt, b0=idb[0])
    return zen, cnt, idb[0]


def find_clean_period(t0_unix, t1_unix, cache, margin=0):
    """The first period in the span on the zenith antenna throughout, with no
    catalogued echo on this channel in it or within margin periods of it."""
    zen, cnt, b0 = survey_periods(cache)
    step = int(LPI["avg_dur"] * 1e6)
    a0 = max(margin, int((t0_unix * 1e6 - b0) // step))
    a1 = min(len(zen) - margin, int((t1_unix * 1e6 - b0) // step))
    for ai in range(a0, a1):
        near = slice(ai - margin, ai + margin + 1)
        if zen[near].all() and cnt[near].sum() == 0:
            return ai, ai * step + b0
    raise RuntimeError("no clean period found between %d and %d: %d zenith periods, %d of them clean"
                       % (t0_unix, t1_unix, zen[a0:a1].sum(), (zen[a0:a1] & (cnt[a0:a1] == 0)).sum()))


class InjectingReader:
    """DigitalRFReader that adds a synthetic satellite echo to chosen pulses."""

    def __init__(self, path, plan, sid):
        self.reader = DigitalRFReader(path)
        self.plan = plan
        self.sid = sid

    def __getattr__(self, name):
        return getattr(self.reader, name)

    def read_vector_1d(self, key, n_samples, channel):
        v = self.reader.read_vector_1d(key, n_samples, channel).astype(n.complex64, casting="unsafe")
        if channel != CHANNEL or key not in self.plan:
            return v
        p = self.plan[key]
        tm = TMM[self.sid[key]]
        u = n.zeros(n_samples, dtype=n.complex128)
        u[tm["tx0"]:tm["tx1"]] = v[tm["tx0"]:tm["tx1"]] - Z_DC
        e_u = n.sum(n.abs(u)**2)
        # scale so every pulse has the same echo energy, P_true
        amp = n.sqrt(p["p_true"] / e_u)
        t = n.arange(n_samples)
        echo = amp * n.exp(1j * p["phase"]) * fractional_shift(u, p["delay"]) \
            * n.exp(2j * n.pi * p["doppler"] * t / 1e6)
        return (v + echo).astype(n.complex64)


class PlannedDetections:
    """The injected echoes, as the inversion's satellite_detections."""

    def __init__(self, plan, doppler_override=None):
        self.plan = plan
        self.doppler_override = doppler_override

    def for_period(self, i0, i1):
        out = {}
        for k, p in self.plan.items():
            if i0 <= k < i1:
                fd = p["doppler"] if self.doppler_override is None else self.doppler_override
                out[k] = [(p["delay"], fd)]
        return out


def raw_noise_power(i0, keys, sid):
    """Noise power per raw sample, from the quiet window before the injection."""
    r = DigitalRFReader("%s/rf_data/" % DATA)
    pw = []
    for k in keys[:60]:
        le = TMM[sid[k]]["last_echo"]
        v = r.read_vector_1d(k, 10000, CHANNEL).astype(n.complex64, casting="unsafe")
        pw.append(n.mean(n.abs(v[(le - 500):le] - Z_DC)**2))
    return float(n.median(pw))


def make_plan(i0, snr_db, seed=1):
    idr = DigitalMetadataReader("%s/metadata/id_metadata" % DATA)
    sid = idr.read(i0, i0 + int(LPI["avg_dur"] * 1e6) + 40000, "sweepid")
    sid = {int(k): int(v) for k, v in sid.items()}
    keys = sorted(k for k in sid if sid[k] in TMM)
    n_raw = raw_noise_power(i0, [k for k in keys if (k - i0) / 1e6 < SAT_WINDOW_S[0] - 0.5], sid)
    # matched filter SNR = echo energy / noise power per sample
    p_true = 10**(snr_db / 10) * n_raw
    rdot = -SAT_DOPPLER_HZ * WAVELENGTH / 2.0
    d0 = SAT_RANGE_KM * 1e3 * 2 / sc.c * 1e6
    tc = i0 + 1e6 * 0.5 * (SAT_WINDOW_S[0] + SAT_WINDOW_S[1])
    rng = n.random.default_rng(seed)
    plan = {}
    for k in keys:
        ts = (k - i0) / 1e6
        if SAT_WINDOW_S[0] <= ts < SAT_WINDOW_S[1]:
            plan[k] = dict(delay=d0 + 2 * rdot / sc.c * 1e6 * (k - tc) / 1e6,
                           doppler=SAT_DOPPLER_HZ, p_true=p_true,
                           phase=rng.uniform(0, 2 * n.pi))
    return plan, sid, n_raw, p_true, d0


def run(name, out, ai, reader_factory=None, sat=None, rejection=False, **kw):
    h5 = [f for f in os.listdir("%s/%s/lpi_%d/%s" % (out, name, LPI["rg"], CHANNEL))] \
        if os.path.isdir("%s/%s/lpi_%d/%s" % (out, name, LPI["rg"], CHANNEL)) else []
    if any(f.endswith(".h5") for f in h5):
        print("%s: done already" % name)
        return
    olpi.DigitalRFReader = reader_factory if reader_factory is not None else DigitalRFReader
    args = dict(LPI)
    args.update(kw)
    olpi.lpi_files(dirname=DATA, channel=CHANNEL, output_base="%s/%s" % (out, name),
                   periods=[ai], satellite_detections=sat, outlier_rejection=rejection, **args)
    olpi.DigitalRFReader = DigitalRFReader


def load(out, name):
    d = "%s/%s/lpi_%d/%s" % (out, name, LPI["rg"], CHANNEL)
    f = [x for x in os.listdir(d) if x.endswith(".h5")][0]
    with h5py.File("%s/%s" % (d, f), "r") as h:
        return {k: h[k][()] for k in h.keys()}


def report(out, meta):
    R = {k: load(out, k) for k in ["clean", "clean_current", "nohandling", "fit", "current", "fit_dop0"]}
    rgs = R["clean"]["rgs_km"]
    lags = R["clean"]["lags"]
    d0 = meta["d0"]
    gate_delay = n.arange(len(rgs)) * LPI["rg"]
    # gates whose lagged products the echo reaches: within a pulse length of it
    foot = n.abs(gate_delay - d0) < 600
    away = (n.abs(gate_delay - d0) > 1200) & n.isfinite(R["clean"]["acfs_e"]).all(axis=1)
    lines = []

    def compare(x, ref, what):
        a, b, v = R[x]["acfs_e"], R[ref]["acfs_e"], R[ref]["acfs_var"]
        z = (a - b) / n.sqrt(v)
        zf = z[foot][n.isfinite(z[foot])]
        za = z[away][n.isfinite(z[away])]
        lines.append("%-11s vs %-13s  footprint: mean Re %+8.2f sigma, rms |.| %8.2f sigma, max %9.1f | "
                     "elsewhere rms %.3f sigma   (%s)"
                     % (x, ref, n.mean(zf.real), n.sqrt(n.mean(n.abs(zf)**2)), n.max(n.abs(zf)),
                        n.sqrt(n.mean(n.abs(za)**2)), what))

    lines.append("period %d, t0 %d, injected %d pulses, MF SNR %.1f dB, %.2f km, %.0f Hz"
                 % (meta["ai"], meta["i0"] / 1e6, meta["n_inj"], meta["snr_db"], SAT_RANGE_KM, SAT_DOPPLER_HZ))
    lines.append("ACF differences at the %d gates the echo reaches, in units of the reference's standard deviation:"
                 % foot.sum())
    compare("nohandling", "clean", "bias from the satellite, nothing done")
    compare("current", "clean_current", "the pipeline with its outlier rejection")
    compare("fit", "clean", "the method")
    compare("fit_dop0", "clean", "the method, columns at 0 Hz")
    infl = n.sqrt(R["fit"]["acfs_var"][foot] / R["clean"]["acfs_var"][foot])
    lines.append("error bar of the fit at those gates / clean: median %.3f, max %.3f"
                 % (n.nanmedian(infl), n.nanmax(infl)))

    for name in ["fit", "fit_dop0"]:
        x = R[name]["sat_amp_e"]
        keys = R[name]["sat_keys"]
        inj = n.array([k in meta["plan_keys"] for k in keys])
        ratio = n.abs(x[inj]) / meta["p_true"]
        lines.append("%s: %d satellite columns (%d injected pulses). |amplitude| / true echo energy: "
                     "median %.4f, 16-84%% %.4f-%.4f, median over pulses of the lag 1 / lag %d ratio %.4f"
                     % (name, len(keys), inj.sum(), n.nanmedian(ratio), *n.nanpercentile(ratio, [16, 84]),
                        len(lags), n.nanmedian(n.abs(x[inj, 0]) / n.abs(x[inj, -1]))))
        ph = n.unwrap(n.angle(n.nanmean(x[inj] / n.abs(x[inj]), axis=0)))
        slope = n.polyfit(lags, ph, 1)[0]
        f_col = SAT_DOPPLER_HZ if name == "fit" else 0.0
        f_est = f_col - slope / (2 * n.pi)
        lines.append("%s: phase slope %.1f rad/s -> Doppler %.1f Hz (true %.1f), range rate %.2f m/s (true %.2f)"
                     % (name, slope, f_est, SAT_DOPPLER_HZ, -f_est * WAVELENGTH / 2, -SAT_DOPPLER_HZ * WAVELENGTH / 2))
    txt = "\n".join(lines)
    print(txt)
    with open("%s/report.txt" % out, "w") as f:
        f.write(txt + "\n")

    # figure: ACF real part near the satellite at two lags, and the amplitudes
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    near = n.abs(rgs - SAT_RANGE_KM) < 350
    for j, li in enumerate([1, 8]):
        a = ax[0, j]
        sd = n.sqrt(R["clean"]["acfs_var"][near, li])
        a.fill_between(rgs[near], R["clean"]["acfs_e"][near, li].real - sd,
                       R["clean"]["acfs_e"][near, li].real + sd, color="0.8", label="clean ±1σ")
        for name, st in [("clean", "k-"), ("nohandling", "r-"), ("current", "m--"), ("fit", "b.-")]:
            a.plot(rgs[near], R[name]["acfs_e"][near, li].real, st, label=name, lw=1)
        a.axvline(SAT_RANGE_KM, color="g", ls=":", label="satellite")
        lo, hi = n.nanpercentile(R["clean"]["acfs_e"][near, li].real, [2, 98])
        pad = 2 * (hi - lo) + 1e-30
        a.set_ylim(lo - pad, hi + pad)
        a.set_title("Re ACF, lag %.0f µs" % (lags[li] * 1e6))
        a.set_xlabel("range (km)")
        a.legend(fontsize=7)
    x = R["fit"]["sat_amp_e"]
    inj = n.array([k in meta["plan_keys"] for k in R["fit"]["sat_keys"]])
    ax[1, 0].plot(lags * 1e6, (n.abs(x[inj]) / meta["p_true"]).T, color="b", alpha=0.08)
    ax[1, 0].plot(lags * 1e6, n.nanmedian(n.abs(x[inj]) / meta["p_true"], axis=0), "k", lw=2, label="median")
    ax[1, 0].axhline(1, color="g", ls=":")
    ax[1, 0].set_ylim(0, 2)
    ax[1, 0].set_xlabel("lag (µs)")
    ax[1, 0].set_ylabel("|amplitude| / true echo energy")
    ax[1, 0].set_title("fitted satellite amplitude, every pulse")
    ax[1, 0].legend()
    x0 = R["fit_dop0"]["sat_amp_e"]
    inj0 = n.array([k in meta["plan_keys"] for k in R["fit_dop0"]["sat_keys"]])
    ax[1, 1].plot(lags * 1e6, n.unwrap(n.angle(x0[inj0]), axis=1).T, color="b", alpha=0.08)
    ax[1, 1].plot(lags * 1e6, -2 * n.pi * SAT_DOPPLER_HZ * lags, "g:", lw=2, label="expected, %.0f Hz" % SAT_DOPPLER_HZ)
    ax[1, 1].set_xlabel("lag (µs)")
    ax[1, 1].set_ylabel("phase (rad)")
    ax[1, 1].set_title("amplitude phase, columns built at 0 Hz")
    ax[1, 1].legend()
    fig.tight_layout()
    fig.savefig("%s/test1.png" % out, dpi=120)
    print("wrote %s/test1.png" % out)


def metrics(out, meta):
    """The headline numbers of one period and SNR, for the summary table."""
    R = {k: load(out, k) for k in ["clean", "nohandling", "fit", "fit_dop0"]}
    gate_delay = n.arange(len(R["clean"]["rgs_km"])) * LPI["rg"]
    foot = n.abs(gate_delay - meta["d0"]) < 600

    def dev(x):
        z = (R[x]["acfs_e"] - R["clean"]["acfs_e"]) / n.sqrt(R["clean"]["acfs_var"])
        z = z[foot][n.isfinite(z[foot])]
        return n.sqrt(n.mean(n.abs(z)**2)), n.max(n.abs(z))

    x = R["fit_dop0"]["sat_amp_e"]
    inj = n.array([k in meta["plan_keys"] for k in R["fit_dop0"]["sat_keys"]])
    lags = R["clean"]["lags"]
    ph = n.unwrap(n.angle(n.nanmean(x[inj] / n.abs(x[inj]), axis=0)))
    f_est = -n.polyfit(lags, ph, 1)[0] / (2 * n.pi)
    # coherent mean over pulses of the fitted amplitude, per lag, relative to
    # the truth, and its standard deviation. Averaging the complex values
    # avoids the upward bias the magnitude of a noisy number has.
    inj_f = n.array([k in meta["plan_keys"] for k in R["fit"]["sat_keys"]])
    xc = n.abs(n.nanmean(R["fit"]["sat_amp_e"][inj_f], axis=0)) / meta["p_true"]
    sd = n.sqrt(n.nanmean(R["fit"]["sat_amp_var"][inj_f], axis=0) / inj_f.sum()) / meta["p_true"]
    return dict(nh=dev("nohandling"), fit=dev("fit"),
                infl=n.nanmedian(n.sqrt(R["fit"]["acfs_var"][foot] / R["clean"]["acfs_var"][foot])),
                amp=n.nanmedian(xc), flat=n.nanmax(n.abs(xc - 1) / sd), rdot=-f_est * WAVELENGTH / 2)


def summary(root):
    """One line per period and SNR, for every test directory under root."""
    import glob
    rows = []
    for mf in sorted(glob.glob("%s/**/meta.json" % root, recursive=True)):
        out = os.path.dirname(mf)
        meta = json.load(open(mf))
        if "plan_keys" not in meta:
            continue   # not an injection test
        meta["plan_keys"] = set(meta["plan_keys"])
        try:
            m = metrics(out, meta)
        except (FileNotFoundError, IndexError, KeyError):
            continue
        rows.append((meta["ai"], meta["i0"], meta["snr_db"], m))
    rdot_true = -SAT_DOPPLER_HZ * WAVELENGTH / 2
    lines = ["ACF deviation from the clean run at the gates the echo reaches, in σ; error bar of the fit / clean;",
             "fitted echo energy / truth (median over lags), and its largest deviation from 1 over the 45 lags in σ",
             "%-6s %-17s %4s | %-19s %-19s | %6s | %6s %8s | %s"
             % ("period", "start UTC", "SNR", "no handling rms/max", "fit rms/max", "σ fit",
                "energy", "flat", "range rate (true %.2f m/s)" % rdot_true)]
    import datetime
    for ai, i0, snr, m in sorted(rows, key=lambda r: (r[2], r[0])):
        t = datetime.datetime.fromtimestamp(i0 / 1e6, datetime.timezone.utc).strftime("%m-%d %H:%M:%S")
        lines.append("%-6d %-17s %4.0f | %6.2f / %6.2f σ     %6.2f / %6.2f σ     | %6.3f | %6.3f %6.2f σ | %.2f"
                     % (ai, t, snr, *m["nh"], *m["fit"], m["infl"], m["amp"], m["flat"], m["rdot"]))
    txt = "\n".join(lines)
    print(txt)
    with open("%s/summary.txt" % root, "w") as f:
        f.write(txt + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--snr-db", type=float, default=30.0, help="matched filter SNR of the injected echo")
    p.add_argument("--out", default="plots/satfit_tests")
    p.add_argument("--search", nargs=2, type=float, default=[1712491200, 1712509200],
                   help="unix time span to look for a clean period in")
    p.add_argument("--runs", default=None,
                   help="comma separated subset of the inversions, to run them in parallel; "
                        "the report is written when all are run in one call")
    p.add_argument("--summary", action="store_true",
                   help="only tabulate every finished test under --out")
    args = p.parse_args()
    if args.summary:
        summary(args.out)
        return
    out = "%s/snr%02.0f" % (args.out, args.snr_db)
    os.makedirs(out, exist_ok=True)

    mf = "%s/meta.json" % out
    if os.path.exists(mf):
        meta = json.load(open(mf))
        meta["plan_keys"] = set(meta["plan_keys"])
        ai, i0 = meta["ai"], meta["i0"]
        plan, sid, _, _, _ = make_plan(i0, args.snr_db)
    else:
        ai, i0 = find_clean_period(*args.search, cache="%s/period_survey.npz" % args.out)
        plan, sid, n_raw, p_true, d0 = make_plan(i0, args.snr_db)
        meta = dict(ai=int(ai), i0=int(i0), snr_db=args.snr_db, n_raw=float(n_raw), p_true=float(p_true), d0=float(d0),
                    n_inj=len(plan), plan_keys=sorted(int(k) for k in plan.keys()))
        json.dump(meta, open(mf, "w"), indent=1)
        meta["plan_keys"] = set(meta["plan_keys"])
    print("period %d starting %d, %d pulses injected" % (ai, i0 / 1e6, len(plan)))

    inject = lambda path: InjectingReader(path, plan, sid)
    runs = {
        "clean": lambda: run("clean", out, ai, rejection=False),
        "clean_current": lambda: run("clean_current", out, ai, rejection=True),
        "nohandling": lambda: run("nohandling", out, ai, inject, rejection=False),
        "fit": lambda: run("fit", out, ai, inject, sat=PlannedDetections(plan), rejection=False),
        "current": lambda: run("current", out, ai, inject, rejection=True),
        "fit_dop0": lambda: run("fit_dop0", out, ai, inject,
                                sat=PlannedDetections(plan, doppler_override=0.0), rejection=False),
    }
    for name in (runs if args.runs is None else [r for r in args.runs.split(",") if r]):
        runs[name]()
    if args.runs is None:
        report(out, meta)


if __name__ == "__main__":
    main()
