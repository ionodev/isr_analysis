#!/usr/bin/env python3
"""
Compare the measured Millstone Hill system temperature with a 440 MHz sky model.

Takes the output of tsys_harvest.py and the sky model of sky_noise_model.py and
fits, for each antenna,

    T_sys = T_0 + eta * T_sky

T_0 collects everything that does not move with the sky: receiver noise, ground
and far sidelobe spillover, atmospheric emission, and any steady interference.
eta is the slope of the measurement against the model.

What eta means
--------------
eta is not the main beam efficiency alone. T_sys is derived from the injected
noise temperature as T_sys = noise / alpha with alpha proportional to
1 / T_INJECTION, so the whole measured series scales linearly with the assumed
T_INJECTION = 1172 K. If that constant is wrong by a factor k, the fitted slope
is k times the true beam efficiency, and the two cannot be separated by this fit
alone. What the fit tests is whether the product is consistent with unity, and
if not, by how much.

Which T_sys is used
-------------------
range_doppler_300_outlier, not lpi_30. The lag profile inversion's noise
estimate is biased low by a signal dependent factor; see lpi_bias_factor below.
The LPI series is still read, corrected with that formula, and compared, because
the two then agreeing is a strong check on both.

Which antenna the result rests on
---------------------------------
The zenith antenna. MISA is fitted and reported, but its slope is not a
measurement of anything; see UNUSABLE_ANTENNAS.

Bright sources
--------------
The model is a diffuse sky map at 1 degree resolution and contains no discrete
sources at their true brightness. Periods when a strong source is within a few
beamwidths are therefore masked rather than fitted. This matters most for the
zenith antenna, whose boresight passes 0.06 degrees from Cygnus A: at transit
T_sys reaches 1600 K against a 170 K baseline, ten times the diffuse signal the
fit is trying to measure.

Rejecting interference
----------------------
Interference is rejected by its shape in time, not by the size of its residual.
A robust fit against residual size would discard the galactic plane transits,
which are large, real, and carry nearly all the information about the slope. The
test used here instead asks whether a period disagrees with its immediate
neighbours: a satellite echo or an RFI burst lasts one or two ten second periods
and leaves its neighbours untouched, whereas a source transit takes minutes and
carries its neighbours with it. After that cut the fit is ordinary least
squares, so the transits keep full weight.

    python3 compare_tsys_sky.py tsys.h5 [-o plots/] [--model gsm2016]
"""

import numpy as n
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py
import os
import argparse

import sky_noise_model as snm
from radar_timing import T_INJECTION

# The estimator to fit. range_doppler duplicates range_doppler_300_outlier
# exactly and range_doppler_300 is a partial rerun of the same periods, so
# including either would weight those periods more than once.
PRIMARY_SET = "range_doppler_300_outlier"
# Carried through the bias correction and compared, not stacked with the above.
SECONDARY_SET = "lpi_30"

# An integration period outside this range did not measure a temperature, it
# failed. Wide enough to keep the Cygnus A transit, which is real.
T_SYS_VALID_K = (50.0, 5000.0)

# Below this elevation the ground fills enough of the beam that a model of the
# main beam alone is not a useful description. This also excludes the mode 800
# horizon scan, which sits at 6 degrees.
MIN_ELEVATION_DEG = 20.0

# Strong 440 MHz sources, and how close to the beam they may come before the
# period is masked, in beamwidths. The model map does not carry them at their
# true brightness, so any period where one is in the beam is unfittable.
BRIGHT_SOURCES = {
    "Cyg A": (299.8681, 40.7339),
    "Cas A": (350.8575, 58.8115),
    "Tau A": (83.6331, 22.0145),
    "Vir A": (187.7059, 12.3911),
}
SOURCE_MASK_BEAMWIDTHS = 4.0

# Neighbour test for transient interference. A period is rejected if it differs
# from the median of its neighbours by more than NEIGHBOUR_SIGMA times the local
# scale. The window is short so that a real transit, which neighbours follow, is
# never rejected; it is a few minutes long and the window is under a minute.
NEIGHBOUR_HALFWIDTH = 3          # periods either side
NEIGHBOUR_SIGMA = 6.0
SCALE_HALFWIDTH = 30             # periods either side, for the local scale
# Samples further apart than this are not neighbours; there are data gaps.
NEIGHBOUR_MAX_GAP_S = 60.0

# Antennas whose fit is reported but must not be read as a result, and why.
UNUSABLE_ANTENNAS = {
    "misa": """MISA alternates its azimuth between integration periods, and the
    antenna control metadata does not locate the switches well enough to say
    which azimuth a given ten second period belongs to. The model therefore
    swings by tens of kelvin from one period to the next while the measurement
    stays flat, which is visible directly in the time series plot. Error in the
    model sky temperature is error in the regressor, and that drags a least
    squares slope towards zero, which is the direction MISA's slope moves. Its
    two estimators disagree (0.48 against 0.63) where the zenith antenna's agree
    to one per cent, and range_doppler_300_outlier additionally spans elevation
    37 to 45 degrees, mid slew, where lpi_30 is parked at 45. Fixing this needs
    the pointing resolved per pulse rather than per integration period, which is
    not attempted here.""",
}

# Block length for the bootstrap, seconds, and number of resamples.
BOOTSTRAP_BLOCK_S = 1800.0
N_BOOTSTRAP = 500

# Length of the background window used by both noise estimates, seconds.
BG_WINDOW_S = 500e-6


def lpi_bias_factor(t_est, pass_band_hz, t_inj=T_INJECTION, bg_window_s=BG_WINDOW_S):
    """
    Undo the bias in outlier_lpi.py's system temperature.

    That estimator takes the mean of the 500 microsecond background window as a
    DC offset and subtracts it from both the background and the injection
    window. The mean is not a DC offset, it is an estimate of one from N
    independent samples, where N is the filtered bandwidth times the window
    length. Subtracting it removes 1/N of the background window's own variance,
    and adds an uncorrelated |mean|^2 = P_bg/N to the injection window, which
    the background mean has nothing to do with. So

        T_est = T_inj P_bg (1 - 1/N) / [ (P_inj - P_bg) + 2 P_bg / N ]

    and writing y = T_inj / T_true for the true ratio (P_inj - P_bg) / P_bg,

        T_true = T_inj / [ T_inj (1 - 1/N) / T_est  -  2/N ]

    The filter cutoff is 1.2 * pass_band, and the passband of a complex baseband
    signal is twice the cutoff, so N = 2 * 1.2 * pass_band * 500 us. At the
    18 kHz passband used for this experiment N is only 21.6, not the 500 samples
    the window holds, which is why the bias is large.

    Measured against range_doppler_300_outlier on the same periods, this turns
    159 K into 169 K at baseline where the unbiased value is 170 K, and 1365 K
    into 1613 K at the Cygnus A peak where the unbiased value is 1607 K.
    """
    n_eff = 2.0 * 1.2 * pass_band_hz * bg_window_s
    denom = t_inj * (1.0 - 1.0 / n_eff) / t_est - 2.0 / n_eff
    return n.where(denom > 0, t_inj / n.where(denom > 0, denom, 1.0), n.nan)


def load(fname):
    with h5py.File(fname, "r") as h:
        d = {k: h[k][()] for k in h.keys()}
    for k in ["set", "channel", "antenna"]:
        d[k] = n.array(d[k]).astype(str)
    return d


def source_separations(t_unix, az_deg, el_deg):
    """Angular distance from the beam to each bright source, degrees."""
    from astropy.coordinates import EarthLocation, AltAz, SkyCoord
    from astropy.time import Time
    import astropy.units as u

    loc = EarthLocation(lat=snm.RADAR_LAT * u.deg, lon=snm.RADAR_LON * u.deg,
                        height=snm.RADAR_HGT * u.m)
    beam = SkyCoord(AltAz(az=az_deg * u.deg, alt=el_deg * u.deg,
                          obstime=Time(t_unix, format="unix"),
                          location=loc)).icrs
    out = {}
    for name, (ra, dec) in BRIGHT_SOURCES.items():
        out[name] = beam.separation(SkyCoord(ra=ra * u.deg, dec=dec * u.deg)).deg
    return out


def neighbour_outliers(t, y):
    """
    Periods that disagree with their immediate neighbours.

    Returns a boolean mask of the periods to reject. Both inputs must be sorted
    in time.
    """
    m = len(y)
    bad = n.zeros(m, dtype=bool)
    if m < 2 * SCALE_HALFWIDTH + 1:
        return bad
    for i in range(m):
        lo = max(0, i - NEIGHBOUR_HALFWIDTH)
        hi = min(m, i + NEIGHBOUR_HALFWIDTH + 1)
        idx = n.array([j for j in range(lo, hi)
                       if j != i and abs(t[j] - t[i]) <= NEIGHBOUR_MAX_GAP_S])
        if len(idx) < 2:
            continue
        slo = max(0, i - SCALE_HALFWIDTH)
        shi = min(m, i + SCALE_HALFWIDTH + 1)
        # scale from successive differences, which a smooth transit does not
        # inflate the way a standard deviation about a local mean would
        dy = n.abs(n.diff(y[slo:shi]))
        scale = 1.4826 * n.median(dy) / n.sqrt(2.0)
        if not n.isfinite(scale) or scale <= 0:
            continue
        if n.abs(y[i] - n.median(y[idx])) > NEIGHBOUR_SIGMA * scale:
            bad[i] = True
    return bad


def select(d, which_set, antenna):
    """Rows of one antenna and one output set worth fitting, and what was cut."""
    keep = (d["set"] == which_set) & (d["antenna"] == antenna)
    cuts = {}

    def drop(mask, why):
        nonlocal keep
        cuts[why] = int((keep & ~mask).sum())
        keep = keep & mask

    drop(n.isfinite(d["T_sys"]), "T_sys not finite")
    drop((d["T_sys"] > T_SYS_VALID_K[0]) & (d["T_sys"] < T_SYS_VALID_K[1]),
         "T_sys outside %.0f-%.0f K" % T_SYS_VALID_K)
    drop(n.isfinite(d["az"]) & n.isfinite(d["el"]), "no pointing")
    drop(d["el"] > MIN_ELEVATION_DEG, "elevation below %.0f deg" % MIN_ELEVATION_DEG)
    # 1 = MISA, -1 = zenith, as millstone_radar_state defines it. A channel whose
    # own antenna was not the selected receiver is not measuring that antenna's
    # sky, whatever the channel is named.
    drop(d["rx_antenna"] == (1.0 if antenna == "misa" else -1.0),
         "receiver switched to the other antenna")
    return keep, cuts


def fit_antenna(d, keep, antenna, model, freq_mhz, correct_lpi_pass_band=None):
    """Sky model, masking, interference rejection and least squares."""
    sm = snm.SkyModel(freq_mhz=freq_mhz, model=model, antenna=antenna)

    idx = n.where(keep)[0]
    idx = idx[n.argsort(d["tmid"][idx])]
    t, az, el = d["tmid"][idx], d["az"][idx], d["el"][idx]
    t_sys = d["T_sys"][idx].astype(float)

    if correct_lpi_pass_band is not None:
        t_sys = lpi_bias_factor(t_sys, correct_lpi_pass_band)

    t_sky = sm.at_horizontal(t, az, el)

    cuts = {}
    seps = source_separations(t, az, el)
    masked = n.zeros(len(t), dtype=bool)
    for name, sep in seps.items():
        hit = sep < SOURCE_MASK_BEAMWIDTHS * sm.beam_fwhm_deg
        cuts["within %.1f beamwidths of %s" % (SOURCE_MASK_BEAMWIDTHS, name)] = int(hit.sum())
        masked |= hit

    transient = neighbour_outliers(t, t_sys)
    cuts["disagrees with its neighbours"] = int((transient & ~masked).sum())

    good = n.isfinite(t_sys) & ~masked & ~transient
    tg, xg, yg = t[good], t_sky[good], t_sys[good]

    A = n.vstack([n.ones(len(xg)), xg]).T
    beta = n.linalg.lstsq(A, yg, rcond=None)[0]
    resid = yg - A.dot(beta)

    # uncertainty from a moving block bootstrap: consecutive periods see almost
    # the same sky, so the residuals are strongly autocorrelated and the least
    # squares covariance would be optimistic by more than an order of magnitude
    rng = n.random.default_rng(0)
    block = n.floor((tg - tg.min()) / BOOTSTRAP_BLOCK_S).astype(int)
    ids = n.unique(block)
    idx_of = {b: n.where(block == b)[0] for b in ids}
    boot = []
    for _ in range(N_BOOTSTRAP):
        sel = n.concatenate([idx_of[b] for b in rng.choice(ids, len(ids), replace=True)])
        Ab = n.vstack([n.ones(len(sel)), xg[sel]]).T
        boot.append(n.linalg.lstsq(Ab, yg[sel], rcond=None)[0])
    boot = n.array(boot)

    return dict(antenna=antenna, sky=sm, cuts=cuts,
                t=tg, t_sky=xg, t_sys=yg, resid=resid,
                t_all=t, t_sys_all=t_sys, t_sky_all=t_sky, good=good,
                T_0=beta[0], dT_0=boot[:, 0].std(),
                eta=beta[1], deta=boot[:, 1].std(),
                rms=float(n.std(resid)),
                var_explained=1.0 - n.var(resid) / n.var(yg),
                n=len(yg), n_blocks=len(ids))


def report(f, label):
    print("\n%s / %s: %d periods fitted, %d bootstrap blocks"
          % (f["antenna"], label, f["n"], f["n_blocks"]))
    if f["antenna"] in UNUSABLE_ANTENNAS:
        print("  *** not a measurement:")
        for line in UNUSABLE_ANTENNAS[f["antenna"]].split("\n"):
            print("      %s" % (line.strip()))
    for why, cnt in f["cuts"].items():
        if cnt > 0:
            print("    masked %5d: %s" % (cnt, why))
    print("  T_0 = %7.1f +- %4.1f K    receiver + spillover + atmosphere" % (f["T_0"], f["dT_0"]))
    print("  eta = %7.3f +- %5.3f      calibration scale x beam efficiency" % (f["eta"], f["deta"]))
    print("  model T_sky %.1f to %.1f K, residual rms %.1f K, variance explained %.2f"
          % (f["t_sky"].min(), f["t_sky"].max(), f["rms"], f["var_explained"]))
    if f["eta"] > 0 and f["antenna"] not in UNUSABLE_ANTENNAS:
        print("  eta is the product of beam efficiency and any error in"
              " T_INJECTION = %.0f K;" % (T_INJECTION))
        print("  this fit cannot separate them. A perfect beam would require"
              " T_INJECTION = %.0f K." % (T_INJECTION / f["eta"]))


def plot(groups, outdir, model):
    os.makedirs(outdir, exist_ok=True)
    ncol = len(groups)

    fig, axes = plt.subplots(2, ncol, figsize=(7.0 * ncol, 9.0), squeeze=False)
    for j, (label, f) in enumerate(groups):
        ax = axes[0, j]
        th = (f["t_all"] - f["t_all"].min()) / 3600.0
        ax.plot(th[~f["good"]], f["t_sys_all"][~f["good"]], ".", ms=2.0,
                color="tab:orange", alpha=0.7, label="masked")
        ax.plot(th[f["good"]], f["t_sys_all"][f["good"]], ".", ms=1.2,
                color="0.3", alpha=0.4, label="measured")
        ax.plot(th, f["T_0"] + f["eta"] * f["t_sky_all"], "-", lw=1.0,
                color="crimson", label="%.1f K + %.2f T$_{sky}$" % (f["T_0"], f["eta"]))
        ax.set_yscale("log")
        ax.set_xlabel("hours since start of recording")
        ax.set_ylabel("T (K)")
        ax.set_title("%s  (%s)" % (label, model))
        ax.legend(loc="upper right", markerscale=6, fontsize=8)
        ax.grid(alpha=0.3)

        ax = axes[1, j]
        ax.plot(f["t_sky"], f["t_sys"], ".", ms=1.5, alpha=0.3, color="0.3")
        xs = n.linspace(f["t_sky"].min(), f["t_sky"].max(), 10)
        ax.plot(xs, f["T_0"] + f["eta"] * xs, "-", color="crimson", lw=2,
                label="T$_0$=%.1f$\\pm$%.1f K\n$\\eta$=%.3f$\\pm$%.3f"
                      % (f["T_0"], f["dT_0"], f["eta"], f["deta"]))
        ax.plot(xs, f["T_0"] + xs, "--", color="tab:blue", lw=1.2, label="$\\eta$ = 1")
        ax.set_xlabel("model T$_{sky}$ (K)")
        ax.set_ylabel("measured T$_{sys}$ (K)")
        ax.legend(loc="upper left", fontsize=9)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    p = os.path.join(outdir, "tsys_sky_%s.png" % model)
    fig.savefig(p, dpi=140)
    plt.close(fig)

    # folded on sidereal time: the sky repeats every sidereal day, the receiver
    # and any solar or diurnal interference do not, so the scatter between the
    # two days at one sidereal hour measures the non sky terms
    fig, axes = plt.subplots(1, ncol, figsize=(6.5 * ncol, 5.0), squeeze=False)
    for j, (label, f) in enumerate(groups):
        ax = axes[0, j]
        lst = snm.local_sidereal_time_h(f["t"])
        ax.plot(lst, f["t_sys"], ".", ms=1.5, alpha=0.3, color="0.3", label="measured")
        ax.plot(lst, f["T_0"] + f["eta"] * f["t_sky"], ".", ms=1.5, alpha=0.6,
                color="crimson", label="model")
        ax.set_xlabel("local apparent sidereal time (h)")
        ax.set_ylabel("T (K)")
        ax.set_title(label)
        ax.set_xlim([0, 24])
        ax.legend(loc="upper right", markerscale=6)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    q = os.path.join(outdir, "tsys_lst_%s.png" % model)
    fig.savefig(q, dpi=140)
    plt.close(fig)
    print("\nwrote %s and %s" % (p, q))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsys", help="output of tsys_harvest.py")
    ap.add_argument("-o", "--outdir", default="tsys_sky")
    ap.add_argument("--model", default="gsm2008", choices=["gsm2008", "gsm2016"])
    ap.add_argument("--freq-mhz", type=float, default=snm.DEFAULT_FREQ_MHZ)
    ap.add_argument("--no-lpi", action="store_true",
                    help="skip the bias corrected lag profile inversion comparison")
    args = ap.parse_args()

    d = load(args.tsys)
    print("read %d integration periods" % (len(d["T_sys"])))

    groups = []
    for antenna in ["zenith", "misa"]:
        keep, cuts = select(d, PRIMARY_SET, antenna)
        if keep.sum() < 100:
            print("\n%s: only %d periods in %s, skipping" % (antenna, keep.sum(), PRIMARY_SET))
            continue
        print("\n%s / %s: %d periods before masking" % (antenna, PRIMARY_SET, keep.sum()))
        for why, cnt in cuts.items():
            if cnt > 0:
                print("    dropped %5d: %s" % (cnt, why))
        f = fit_antenna(d, keep, antenna, args.model, args.freq_mhz)
        report(f, PRIMARY_SET)
        groups.append(("%s, %s" % (antenna, PRIMARY_SET), f))

    if not args.no_lpi:
        print("\n" + "=" * 72)
        print("Bias corrected lag profile inversion, as a check on both estimators")
        print("=" * 72)
        for antenna in ["zenith", "misa"]:
            keep, _ = select(d, SECONDARY_SET, antenna)
            if keep.sum() < 100:
                continue
            pb = n.nanmedian(d["pass_band"][keep])
            if not n.isfinite(pb):
                print("\n%s: no pass_band recorded, cannot correct" % antenna)
                continue
            f = fit_antenna(d, keep, antenna, args.model, args.freq_mhz,
                            correct_lpi_pass_band=pb)
            report(f, "%s corrected, pass_band %.0f Hz" % (SECONDARY_SET, pb))

    if groups:
        plot(groups, args.outdir, args.model)


if __name__ == "__main__":
    main()
