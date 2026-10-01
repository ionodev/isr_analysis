#!/usr/bin/env python3
"""
Test 2 of the satellite column method: real, catalogued satellite transits.

Two kinds of transit are chosen from the single pulse catalogue, each on the
zenith antenna with the transmitter on and lying wholly inside one integration
period:

    short   one object, 0.5-3 s, SNR >= 20 dB. The recording has no bright short
            transit: a bright satellite is detected far into the sidelobes, so
            it makes a long pass.
    long    one object, 3-8.5 s, SNR >= 30 dB: the common case for a bright
            satellite, and the hard one for a free amplitude per pulse, since
            the satellite is then in much of the period.

A transit matched to the space object catalogue by the hard target survey is
preferred. Delays are converted into the inversion's frame with the channel
delay tx_delay.py measures at the start of that period.

The 550 km shell is occupied almost continuously, so the neighbouring periods
always hold echoes near the transit's range. The reference is therefore those
four periods (10 and 20 s either side) masked: inverted with every pulse that
carries a catalogued echo dropped. The period itself is inverted five ways:

    nohandling  no columns, no outlier rejection
    current     no columns, outlier rejection            the pipeline as it was
    masked      pulses with a catalogued echo dropped    masking
    fit         a column for every catalogued echo       the method
    fit_dop0    as fit, columns built at 0 Hz             range rate from the fit alone

    python3 validate_satellite_real.py --kind short --select
    python3 validate_satellite_real.py --kind short --runs fit     (one process per run)
    python3 validate_satellite_real.py --kind short --report
"""

import os
for _v in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ[_v] = "1"
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
import argparse
import datetime
import json
import multiprocessing

import numpy as n
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from digital_rf import DigitalMetadataReader

import validate_satellite_lpi as v1
from validate_satellite_lpi import DATA, CHANNEL, LPI, WAVELENGTH, run, load
import outlier_lpi as olpi
import millstone_radar_state as mrs
import satellite_columns as satcol
from radar_timing import TMM
import tx_delay

TLE_MATCH = os.path.expanduser(
    "~/isr_project/isr_analysis_hardtarget_interference/plots/hard_targets/tle_match_zenith-l.h5")
STEP = int(LPI["avg_dur"] * 1e6)
SAT_MD = "%s/metadata/satellite_detections" % DATA

KINDS = {
    "short": dict(dur=(0.5, 3.0), margin=1.0, snr=20.0),
    "long": dict(dur=(3.0, 8.5), margin=0.5, snr=30.0),
}


def _echoes(rec):
    """(range, Doppler, SNR) of each echo of one catalogue record on this channel."""
    ch = rec["channel"].decode() if isinstance(rec["channel"], bytes) else rec["channel"]
    if ch != CHANNEL:
        return []
    return list(zip(rec["range_km"], rec["doppler_hz"], rec["snr_db"]))


def period_echo_keys(ai, b0):
    """Every pulse of a period with a catalogued echo on this channel."""
    i0 = b0 + ai * STEP
    sat = DigitalMetadataReader(SAT_MD)
    return sorted(int(k) for k, rec in sat.read(i0, i0 + STEP).items() if _echoes(rec))


def _transit_stats(args):
    """The catalogued echoes of one period, summarised as a single run."""
    ai, b0 = args
    i0 = b0 + ai * STEP
    rows = []
    for k, rec in DigitalMetadataReader(SAT_MD).read(i0, i0 + STEP).items():
        rows += [((k - i0) / 1e6, r, f, s) for r, f, s in _echoes(rec)]
    if len(rows) < 40:
        return None
    t, r, fd, snr = map(n.array, zip(*sorted(rows)))
    return dict(ai=ai, i0=i0, n=len(t), t0=t[0], t1=t[-1], max_gap=n.max(n.diff(t)),
                r_spread=n.percentile(r, 95) - n.percentile(r, 5), r_med=float(n.median(r)),
                fd_med=float(n.median(fd)), snr_med=float(n.median(snr)), snr_max=float(n.max(snr)))


def select_transit(out, kind):
    c = KINDS[kind]
    zen, cnt, b0 = v1.survey_periods("%s/period_survey.npz" % os.path.dirname(out))
    cand = [ai for ai in range(2, len(zen) - 2) if zen[ai - 2:ai + 3].all() and 40 <= cnt[ai] <= 1200]
    print("%d candidate periods; reading their echoes on %d CPUs" % (len(cand), multiprocessing.cpu_count()))
    with multiprocessing.get_context("fork").Pool(multiprocessing.cpu_count()) as pool:
        stats = [s for s in pool.map(_transit_stats, [(ai, b0) for ai in cand]) if s is not None]
    good = [s for s in stats if s["r_spread"] < 20 and s["max_gap"] < 0.5
            and s["t0"] > c["margin"] and s["t1"] < LPI["avg_dur"] - c["margin"]
            and c["dur"][0] <= s["t1"] - s["t0"] <= c["dur"][1] and s["snr_med"] >= c["snr"]]
    with h5py.File(TLE_MATCH, "r") as h:
        tt, tr, tm = h["t_unix"][()], h["r_km"][()], h["matched"][()]
        tn = [x.decode().strip() for x in h["match_name"][()]]
    for s in good:
        tc = s["i0"] / 1e6 + 0.5 * (s["t0"] + s["t1"])
        j = n.where((n.abs(tt - tc) < 10) & (n.abs(tr - s["r_med"]) < 30) & (tm == 1))[0]
        s["tle"] = tn[j[0]] if len(j) else None
    # a catalogue match first, then the brightest
    good.sort(key=lambda s: (s["tle"] is None, -s["snr_med"]))
    print("%d periods pass the %s criteria; the best:" % (len(good), kind))
    for s in good[:6]:
        print("  period %d  %.2f-%.2f s  %d echoes  %.1f km  %.0f Hz  SNR median %.1f max %.1f dB  %s"
              % (s["ai"], s["t0"], s["t1"], s["n"], s["r_med"], s["fd_med"], s["snr_med"], s["snr_max"], s["tle"]))
    if not good:
        raise RuntimeError("no transit passes")
    return good[0], b0


def plan(out, kind):
    s, b0 = select_transit(out, kind)
    ai, i0 = s["ai"], s["i0"]
    # the channel delay of this period, from the coded pulses at its start
    delay_us, spread_us, n_used = tx_delay.estimate_channel_delay(DATA, CHANNEL, t0_unix=i0 / 1e6, verbose=False)
    if delay_us is None:
        raise RuntimeError("channel delay could not be measured")
    idr = DigitalMetadataReader("%s/metadata/id_metadata" % DATA)
    sid = {int(k): int(x) for k, x in idr.read(i0, i0 + STEP + 40000, "sweepid").items()}
    keys = sorted(k for k in sid if sid[k] in TMM)
    quiet = [k for k in keys if (k - i0) / 1e6 < s["t0"] - 0.3] or keys
    n_raw = v1.raw_noise_power(i0, quiet, sid)
    tw = (i0 / 1e6 + s["t0"] - 0.05, i0 / 1e6 + s["t1"] + 0.05)
    rw = (s["r_med"] - 15, s["r_med"] + 15)
    transit = satcol.CatalogueDetections(SAT_MD, CHANNEL, delay_us, t_window=tw,
                                         range_window=rw).for_period(i0, i0 + STEP)
    meta = dict(kind=kind, ai=int(ai), i0=int(i0), b0=int(b0),
                transit={k: (x if isinstance(x, str) or x is None else float(x)) for k, x in s.items()},
                delay_us=float(delay_us), delay_spread_us=float(spread_us), delay_pulses=int(n_used),
                n_raw=float(n_raw), t_window=tw, range_window=rw,
                transit_keys=sorted(int(k) for k in transit),
                transit_delay=float(n.median([d for x in transit.values() for d, _ in x])),
                echo_keys=period_echo_keys(ai, b0),
                neighbour_echo_keys={str(d): period_echo_keys(ai + d, b0) for d in (-2, -1, 1, 2)})
    json.dump(meta, open("%s/meta.json" % out, "w"), indent=1)
    print("period %d: channel delay %.3f us (+- %.3f, %d pulses); transit on %d pulses; %d pulses with any "
          "echo; neighbours %s" % (ai, delay_us, spread_us, n_used, len(transit), len(meta["echo_keys"]),
                                   {d: len(x) for d, x in meta["neighbour_echo_keys"].items()}))
    return meta


def all_detections(meta, doppler_override=None, refine=False):
    """Every catalogued echo of the period, as a production run would take them;
    refine: with the delays refined below the catalogue's grid"""
    from raw_reader import RawReader
    det = satcol.CatalogueDetections(SAT_MD, CHANNEL, meta["delay_us"],
                                     refine_reader=RawReader("%s/rf_data/" % DATA) if refine else None)
    if doppler_override is None:
        return det

    class Override:
        def for_period(self, i0, i1):
            return {k: [(d, doppler_override) for d, _ in x] for k, x in det.for_period(i0, i1).items()}
    return Override()


def run_masked(name, out, ai, mask):
    """Masking: pulses carrying a catalogued echo are dropped, by reporting them
    as not transmitted on this antenna, which the inversion then skips."""
    mask = set(int(k) for k in mask)
    original = mrs.get_antenna_select

    def masked_select(path):
        tx_ant, rx_ant = original(path)
        return (lambda k: 1.0 if int(k) in mask else tx_ant(k)), rx_ant
    olpi.mrs.get_antenna_select = masked_select
    try:
        run(name, out, ai, rejection=False)
    finally:
        olpi.mrs.get_antenna_select = original


def transit_columns(R, meta):
    """Which satellite columns model the chosen transit."""
    tk = set(meta["transit_keys"])
    return n.array([int(k) in tk and abs(d - meta["transit_delay"]) < 60
                    for k, d in zip(R["sat_keys"], R["sat_delay_samples"])])


def report(out, meta, fit="fit", dop0="fit_dop0", nohandling="nohandling", current="current"):
    """fit, dop0: which of the fitted runs to report on, as the solve and the
    satellite template were tested in more than one version"""
    ai = meta["ai"]
    R = {k: load(out, v) for k, v in [("nohandling", nohandling), ("current", current),
                                      ("masked", "masked"), ("fit", fit), ("fit_dop0", dop0)]}
    N = [load(out, "neighbour%+d" % d) for d in (-2, -1, 1, 2)]
    rgs = R["fit"]["rgs_km"]
    lags = R["fit"]["lags"]
    gate_delay = n.arange(len(rgs)) * LPI["rg"]
    foot = n.abs(gate_delay - meta["transit_delay"]) < 600
    # control: gates of similar range the echo does not reach
    ctrl = (n.abs(gate_delay - meta["transit_delay"]) > 1200) & (n.abs(gate_delay - meta["transit_delay"]) < 3000)
    nb = n.array([x["acfs_e"] for x in N])
    nb_mean = n.nanmean(nb, axis=0)
    # uncertainty of the neighbour mean from their own scatter, which holds both
    # the noise and the real change of the ionosphere over 40 s
    nb_var = n.nanvar(nb, axis=0, ddof=1) / len(N)
    tr = meta["transit"]
    t_utc = datetime.datetime.fromtimestamp(meta["i0"] / 1e6, datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    lines = ["%s transit: period %d (%s), %.2f-%.2f s into it, %d echoes at %.1f km, %.0f Hz, catalogue SNR "
             "median %.1f dB (max %.1f), TLE match: %s"
             % (meta["kind"], ai, t_utc, tr["t0"], tr["t1"], tr["n"], tr["r_med"], tr["fd_med"],
                tr["snr_med"], tr["snr_max"], tr["tle"]),
             "channel delay %.3f us from %d coded pulses; %d pulses of the period carry a catalogued echo (%.0f%%)"
             % (meta["delay_us"], meta["delay_pulses"], len(meta["echo_keys"]),
                100 * len(meta["echo_keys"]) / 1121.0),
             "",
             "Re and Im of the ACF against the mean of the four masked neighbouring periods, in units of the",
             "combined uncertainty. Control gates are at similar range out of the echo's reach, so they show",
             "the spread to expect when nothing is wrong:"]
    for name in ["nohandling", "current", "masked", "fit"]:
        z = (R[name]["acfs_e"] - nb_mean) / n.sqrt(R[name]["acfs_var"] + nb_var)
        zf, zc = z[foot], z[ctrl]
        zf, zc = zf[n.isfinite(zf)], zc[n.isfinite(zc)]
        lines.append("  %-10s footprint rms %6.2f  max %7.2f  mean Re %+6.2f | control rms %5.2f  max %5.2f"
                     % (name, n.sqrt(n.mean(n.abs(zf)**2)), n.max(n.abs(zf)), n.mean(zf.real),
                        n.sqrt(n.mean(n.abs(zc)**2)), n.max(n.abs(zc))))
    z = (R["fit"]["acfs_e"] - R["masked"]["acfs_e"]) / n.sqrt(R["masked"]["acfs_var"])
    zf = z[foot][n.isfinite(z[foot])]
    lines.append("  fit - masked at the footprint: rms %.2f σ, max %.2f σ; error bar fit / masked: median %.3f"
                 % (n.sqrt(n.mean(n.abs(zf)**2)), n.max(n.abs(zf)),
                    n.nanmedian(n.sqrt(R["fit"]["acfs_var"][foot] / R["masked"]["acfs_var"][foot]))))

    # the satellite
    own = transit_columns(R["fit"], meta)
    x = R["fit"]["sat_amp_e"]
    keys = R["fit"]["sat_keys"]
    rec = DigitalMetadataReader(SAT_MD).read(meta["i0"], meta["i0"] + STEP)

    def cat(k, field):
        # the catalogue echo of this pulse that belongs to the transit
        e = [(r, f, s) for r, f, s in _echoes(rec.get(int(k), {"channel": b"", "range_km": [], "doppler_hz": [], "snr_db": []}))
             if meta["range_window"][0] <= r <= meta["range_window"][1]]
        return e[0][{"doppler": 1, "snr": 2}[field]] if e else n.nan
    snr = n.array([cat(k, "snr") if o else n.nan for k, o in zip(keys, own)])
    fdc = n.array([cat(k, "doppler") if o else n.nan for k, o in zip(keys, own)])
    # the catalogue's matched filter SNR is echo energy over the noise power per
    # raw sample, so it predicts the fitted energy independently of the fit
    e_cat = 10**(snr / 10) * meta["n_raw"]
    e_fit = n.abs(n.nanmean(x[:, 0:10], axis=1))
    ratio = e_fit[own] / e_cat[own]
    lines += ["",
              "satellite: %d columns in the period, %d of them the transit" % (len(keys), own.sum()),
              "  fitted echo energy / energy the catalogue SNR implies: median %.3f, 16-84%% %.3f-%.3f"
              % (n.nanmedian(ratio), *n.nanpercentile(ratio, [16, 84]))]
    # flatness and Doppler per pulse. Averaging over the pulses of a pass would
    # mix the different Dopplers of a sweeping pass and fake a decay with lag.
    ok = n.arange(len(lags)) < 36   # beyond ~360 us the pulse barely overlaps itself
    flat = n.nanmedian(n.abs(x[own][:, ok]) / n.abs(x[own][:, [0]]), axis=0)
    lines.append("  per pulse |amplitude(lag)| / |amplitude(first lag)|, median over pulses, lags <= %.0f us: "
                 "min %.3f max %.3f" % (lags[ok][-1] * 1e6, n.min(flat), n.max(flat)))
    own0 = transit_columns(R["fit_dop0"], meta)
    x0 = R["fit_dop0"]["sat_amp_e"]
    # the Doppler of each pulse from the fit alone, against the catalogue's
    fp = n.array([-n.polyfit(lags[ok], n.unwrap(n.angle(x0[i, ok] * n.conj(x0[i, 0]))), 1)[0] / (2 * n.pi)
                  if n.all(n.isfinite(x0[i, ok])) else n.nan for i in n.where(own0)[0]])
    fdc0 = n.array([cat(k, "doppler") for k in R["fit_dop0"]["sat_keys"][own0]])
    d = fp - fdc0
    lines.append("  per pulse, Doppler from the fit alone (0 Hz columns) minus catalogue Doppler: median %.1f Hz, "
                 "16-84%% %.1f to %.1f Hz; catalogue Doppler over the pass %.0f to %.0f Hz"
                 % (n.nanmedian(d), *n.nanpercentile(d, [16, 84]), n.nanmin(fdc0), n.nanmax(fdc0)))
    # error bar of the fit against masking, gate by gate across the footprint
    for gi in n.where(foot)[0]:
        lines.append("  gate %6.1f km (%+4.0f km from the satellite): error bar fit / masked %6.2f"
                     % (rgs[gi], (gate_delay[gi] - meta["transit_delay"]) * 0.1499,
                        n.nanmedian(n.sqrt(R["fit"]["acfs_var"][gi] / R["masked"]["acfs_var"][gi]))))
    txt = "\n".join(lines)
    print(txt)
    open("%s/report_%s.txt" % (out, fit), "w").write(txt + "\n")

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    near = n.abs(gate_delay - meta["transit_delay"]) < 2400
    for j, li in enumerate([1, 8]):
        a = ax[0, j]
        sd = n.sqrt(nb_var[near, li])
        a.fill_between(rgs[near], nb_mean[near, li].real - 2 * sd, nb_mean[near, li].real + 2 * sd,
                       color="0.8", label="masked neighbours ±2σ")
        for name, st in [("nohandling", "r-"), ("current", "m--"), ("masked", "g-"), ("fit", "b.-")]:
            a.plot(rgs[near], R[name]["acfs_e"][near, li].real, st, label=name, lw=1)
        a.axvline(tr["r_med"], color="k", ls=":", label="satellite")
        lo, hi = n.nanpercentile(nb_mean[near, li].real, [2, 98])
        a.set_ylim(lo - 2 * (hi - lo), hi + 2 * (hi - lo))
        a.set_title("Re ACF, lag %.0f µs" % (lags[li] * 1e6))
        a.set_xlabel("range (km)")
        a.legend(fontsize=7)
    tp = (keys - meta["i0"]) / 1e6
    ax[1, 0].semilogy(tp[own], e_fit[own], "b.", label="fitted echo energy")
    ax[1, 0].semilogy(tp[own], e_cat[own], "r.", ms=3, label="from catalogue SNR")
    ax[1, 0].set_xlabel("time into period (s)")
    ax[1, 0].set_title("light curve of the transit")
    ax[1, 0].legend()
    tq = (R["fit_dop0"]["sat_keys"][own0] - meta["i0"]) / 1e6
    ax[1, 1].plot(tq, fdc0, "r.", ms=3, label="catalogue (matched filter)")
    ax[1, 1].plot(tq, fp, "b.", ms=3, label="lag products, 0 Hz columns")
    ax[1, 1].set_xlabel("time into period (s)")
    ax[1, 1].set_ylabel("Doppler (Hz)")
    ax[1, 1].set_title("Doppler of each pulse")
    ax[1, 1].legend()
    fig.suptitle("%s transit, %s, %s" % (meta["kind"], t_utc, fit))
    fig.tight_layout()
    fig.savefig("%s/test2_%s.png" % (out, fit), dpi=120)
    print("wrote %s/test2_%s.png" % (out, fit))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--kind", choices=sorted(KINDS), required=True)
    p.add_argument("--root", default="plots/satfit_tests")
    p.add_argument("--select", action="store_true")
    p.add_argument("--runs", default=None)
    p.add_argument("--report", action="store_true")
    p.add_argument("--fit", default="fit", help="which fitted run to report on: fit, fit_eq, fit_avg")
    args = p.parse_args()
    out = "%s/real_%s" % (args.root, args.kind)
    os.makedirs(out, exist_ok=True)
    if args.select:
        plan(out, args.kind)
        return
    meta = json.load(open("%s/meta.json" % out))
    ai = meta["ai"]
    runs = {
        "nohandling": lambda: run("nohandling", out, ai, rejection=False),
        "current": lambda: run("current", out, ai, rejection=True),
        "masked": lambda: run_masked("masked", out, ai, meta["echo_keys"]),
        "fit": lambda: run("fit", out, ai, sat=all_detections(meta), rejection=False),
        "fit_dop0": lambda: run("fit_dop0", out, ai, sat=all_detections(meta, 0.0), rejection=False),
        # after the solve was made numerically stable, with each pulse's own
        # template and with the per code average template
        "fit_eq": lambda: run("fit_eq", out, ai, sat=all_detections(meta), rejection=False),
        "fit_dop0_eq": lambda: run("fit_dop0_eq", out, ai, sat=all_detections(meta, 0.0), rejection=False),
        "fit_avg": lambda: run("fit_avg", out, ai, sat=all_detections(meta), rejection=False,
                               satellite_template="average"),
        "fit_dop0_avg": lambda: run("fit_dop0_avg", out, ai, sat=all_detections(meta, 0.0), rejection=False,
                                    satellite_template="average"),
        # with the noise weights in double precision
        "nohandling_p": lambda: run("nohandling_p", out, ai, rejection=False, precise_weights=True),
        "current_p": lambda: run("current_p", out, ai, rejection=True, precise_weights=True),
        "fit_avg_p": lambda: run("fit_avg_p", out, ai, sat=all_detections(meta), rejection=False,
                                 satellite_template="average", precise_weights=True),
        "fit_dop0_avg_p": lambda: run("fit_dop0_avg_p", out, ai, sat=all_detections(meta, 0.0), rejection=False,
                                      satellite_template="average", precise_weights=True),
        # with the catalogue's delays refined below its grid
        "fit_avg_p_refined": lambda: run("fit_avg_p_refined", out, ai, sat=all_detections(meta, refine=True),
                                         rejection=False, satellite_template="average", precise_weights=True),
    }
    for d in (-2, -1, 1, 2):
        runs["neighbour%+d" % d] = (lambda d=d: run_masked("neighbour%+d" % d, out, ai + d,
                                                            meta["neighbour_echo_keys"][str(d)]))
    if args.runs:
        for name in args.runs.split(","):
            runs[name]()
    if args.report:
        dop0 = {"fit": "fit_dop0", "fit_eq": "fit_dop0_eq", "fit_avg": "fit_dop0_avg",
                "fit_avg_p": "fit_dop0_avg_p", "fit_avg_p_refined": "fit_dop0_avg_p"}[args.fit]
        p_ = args.fit.endswith("_p") or args.fit.endswith("_p_refined")
        report(out, meta, fit=args.fit, dop0=dop0, nohandling="nohandling_p" if p_ else "nohandling",
               current="current_p" if p_ else "current")


if __name__ == "__main__":
    main()
