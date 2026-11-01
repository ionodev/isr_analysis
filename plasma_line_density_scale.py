#!/usr/bin/env python3
"""
The density scale from the plasma lines, with outlier rejection, and its
comparison with the digisonde (Henrik, 28 October: "we had actually decided
to use the autopicker with outlier rejection"; it replaces the digisonde
scale of 12 October, Memo 43).

The fits give n_e in relative units; the constant C turns them into m^-3
(fit_quality.DENSITY_CONSTANT). This script measures C from the picks of
plasma_line_picker.py and, on the same footing, from the digisonde, so that
the two constants differ only in their reference:

  1. each pick's frequency f becomes a density. The Langmuir waves' dispersion
     relation (Bohm-Gross, with the magnetic term for a beam at angle alpha to
     the field)
         f^2 = f_p^2 + 3 k^2 (k_B T_e / m_e) / (2 pi)^2 + f_ce^2 sin^2 alpha
     puts the line above the plasma frequency f_p, by 2-5 % in density at
     5-7 MHz and 2500 K; n_e = (f_p / 8.98)^2. T_e is the fit's, at the peak;
     k = 4 pi f_radar / c is the Bragg wave number. The constant is also given
     without the correction, n_e = (f / 8.98)^2, as Memo 43 had it (4-7 %
     above the digisonde).
  2. the fit nearest in time (within MATCH_PICK_S) gives the radar's peak
     density in relative units: a quadratic in log(n_e) against log(range)
     over the good gates (fit_quality.good) within -100..+100 km of the
     highest one, below MAX_PEAK_KM (estimate_magic_constant.py's estimator,
     with Memo 43's arguments -100 100 500). The quadratic must have a
     maximum inside the fitted gates.
  3. C = n_e(plasma line) / peak. Outlying C are rejected iteratively: in
     log C, more than K_SIGMA robust standard deviations (1.4826 times the
     median absolute deviation, the standard deviation for normal scatter)
     from the median. They are bad picks (interference, a ridge lost at
     dusk), bad fits, or a density that changed between the pick and the fit.
     C is the median of the rest; its uncertainty a bootstrap over blocks of
     BLOCK_S (neighbouring minutes are not independent).
  4. the digisonde (GIRO MHJ45, autoscaled): N_mF_2 = 1.24e10 foF2^2 m^-3
     against the same peak estimator, with Memo 43's selection (CS >= 50,
     its hmF2 within 60 km of the radar's, which rejects the autoscaler's
     wrong upper branch) and the same rejection; and, independent of the
     radar's scale, the plasma line's peak density over the digisonde's at
     the same times.

Usage:
    python3 plasma_line_density_scale.py <fits dir> <zenithcal.h5> <digisonde txt> <out dir>

<fits dir> holds lpi_60/zenith-l/pp-*.h5 (fit_quality.load_fits). Writes
summary.txt, picks.txt, digisonde.txt and plasma_line_density_scale.png.
"""

import os
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as n

import fit_quality as fq

RADAR_FREQ_HZ = 440.2e6
C_LIGHT = 299792458.0
K_BRAGG = 4 * n.pi * RADAR_FREQ_HZ / C_LIGHT      # rad/m
KB = 1.380649e-23
ME = 9.1093837015e-31
QE = 1.602176634e-19
# the field in the F region above Millstone (about 45 000 nT at 300 km,
# inclination about 67 degrees, so the zenith beam is about 23 degrees from
# it); approximate, but the whole magnetic term is 0.5 % in density at 7 MHz
B_F_REGION_T = 4.5e-5
ANGLE_TO_B_DEG = 23.0
F_CE_HZ = QE * B_F_REGION_T / (2 * n.pi * ME)

RMIN_KM, RMAX_KM, MAX_PEAK_KM = -100.0, 100.0, 500.0
MIN_GATES = 5
TE_HALF_KM = 30.0
MATCH_PICK_S = 60.0     # picks are minutes, fits 60 s windows
MATCH_DIG_S = 225.0     # half the digisonde's 7.5 min cadence
DIG_MIN_CS = 50
DIG_MAX_DH_KM = 60.0
K_SIGMA = 3.0
BLOCK_S = 1800.0
DAY0 = 1712448000       # 2024-04-07 00:00 UTC
# 8 April 18:16-20:38 UT, the partial eclipse at Westford (digisonde_density_scale.py)
ECLIPSE = (1712600160, 1712608680)


def line_offset_hz2(Te):
    """f^2 - f_p^2 of the Langmuir line (Hz^2): thermal plus magnetic term."""
    thermal = 3 * K_BRAGG ** 2 * KB * n.asarray(Te, dtype=float) / ME / (2 * n.pi) ** 2
    return thermal + (F_CE_HZ * n.sin(n.radians(ANGLE_TO_B_DEG))) ** 2


def plasma_density(f_hz, Te=None):
    """n_e (m^-3) from the line's frequency; without Te uncorrected, (f/8.98)^2."""
    f2 = n.asarray(f_hz, dtype=float) ** 2
    if Te is not None:
        f2 = f2 - line_offset_hz2(Te)
    with n.errstate(invalid="ignore"):
        return n.where(f2 > 0, f2 / 8.98 ** 2, n.nan)


def peak(r, ne, ok, rmin=RMIN_KM, rmax=RMAX_KM, max_km=MAX_PEAK_KM):
    """Peak (density, range km) of a profile: a quadratic in log n_e against
    log range around the highest good gate below max_km. NaN if there is no
    maximum inside the fitted gates."""
    sel = ok & n.isfinite(ne) & (ne > 0) & (r < max_km)
    if sel.sum() < MIN_GATES:
        return n.nan, n.nan
    r0 = r[sel][n.argmax(ne[sel])]
    w = sel & (r - r0 >= rmin) & (r - r0 < rmax)
    if w.sum() < MIN_GATES:
        return n.nan, n.nan
    x, y = n.log10(r[w]), n.log10(ne[w])
    A = n.column_stack([n.ones_like(x), x, x ** 2])
    a = n.linalg.lstsq(A, y, rcond=None)[0]
    if a[2] >= 0:
        return n.nan, n.nan
    xp = -a[1] / (2 * a[2])
    if not (x.min() <= xp <= x.max()):
        return n.nan, n.nan
    return 10 ** (a[0] + a[1] * xp + a[2] * xp ** 2), 10 ** xp


def te_at(r, Te, ok, r_peak, half_km=TE_HALF_KM):
    """Median fitted T_e of the good gates within half_km of r_peak."""
    m = ok & n.isfinite(Te) & (n.abs(r - r_peak) <= half_km)
    return n.median(Te[m]) if m.any() else n.nan


def reject_outliers(x, k=K_SIGMA, max_iter=20):
    """Keep mask: iteratively, |x - median| <= k * 1.4826 * MAD of the kept."""
    x = n.asarray(x, dtype=float)
    keep = n.isfinite(x)
    for _ in range(max_iter):
        med = n.median(x[keep])
        s = 1.4826 * n.median(n.abs(x[keep] - med))
        new = n.isfinite(x) & (n.abs(x - med) <= k * s)
        if (new == keep).all():
            break
        keep = new
    return keep


def block_bootstrap_median(t, x, block_s=BLOCK_S, n_boot=2000, seed=0):
    """Standard deviation of the median over resampled blocks of block_s."""
    b = n.floor((t - t.min()) / block_s).astype(int)
    groups = [x[b == u] for u in n.unique(b)]
    rng = n.random.default_rng(seed)
    med = [n.median(n.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))]))
           for _ in range(n_boot)]
    return n.std(med)


def read_digisonde(fn):
    """GIRO fastchar rows: t, CS, foF2 (MHz), hmF2 (km) (digisonde_density_scale.py)."""
    rows = []
    for line in open(fn):
        if not line.startswith("20"):
            continue
        p = line.split()
        t = n.datetime64(p[0].rstrip("Z")).astype("datetime64[s]").astype(float)
        f = lambda s: n.nan if s.startswith("-") else float(s)
        rows.append((t, float(p[1]), f(p[2]), f(p[4])))
    return n.array(rows)


def fit_peaks(F):
    """Per fit window: centre time, peak (relative), peak range km, T_e at the peak."""
    ok = fq.good(fq.flags(F))
    tc = 0.5 * (F["t0"] + F["t1"])
    pk, hm, te = (n.full(len(tc), n.nan) for _ in range(3))
    for w in range(len(tc)):
        pk[w], hm[w] = peak(F["r"], F["ne"][w], ok[w])
        if n.isfinite(hm[w]):
            te[w] = te_at(F["r"], F["Te"][w], ok[w], hm[w])
    return tc, pk, hm, te


def nearest(t_ref, t, max_dt):
    """Index into t_ref (sorted) of the nearest time to each t, -1 beyond max_dt."""
    j = n.clip(n.searchsorted(t_ref, t), 1, len(t_ref) - 1)
    j = n.where(n.abs(t_ref[j - 1] - t) <= n.abs(t_ref[j] - t), j - 1, j)
    return n.where(n.abs(t_ref[j] - t) <= max_dt, j, -1)


def constant(name, t, C, lines):
    """Rejection, median, bootstrap error and spread of C; one summary line."""
    keep = reject_outliers(n.log10(C))
    med = n.median(C[keep])
    err = block_bootstrap_median(t[keep], C[keep])
    q = n.percentile(C[keep], [16, 84])
    lines.append("%-44s n=%4d (%3d rejected)  C = %.3g +- %.2g (%.1f %%), 16-84 %%: %.3g-%.3g"
                 % (name, keep.sum(), (~keep).sum(), med, err, 100 * err / med, q[0], q[1]))
    return keep, med, err


def main():
    fits_dir, cal, dig_file, out = sys.argv[1:5]
    os.makedirs(out, exist_ok=True)
    F = fq.load_fits(fits_dir, constant=1.0)
    tc, pk, hm, te = fit_peaks(F)
    o = n.argsort(tc)
    tc, pk, hm, te = tc[o], pk[o], hm[o], te[o]

    with h5py.File(cal, "r") as h:
        pf = h["plasma_frequency"][()]
    tp, fp = pf[:, 0], n.abs(pf[:, 1])
    j = nearest(tc, tp, MATCH_PICK_S)
    m = (j >= 0) & n.isfinite(pk[n.maximum(j, 0)]) & n.isfinite(te[n.maximum(j, 0)])
    tp, fp, j = tp[m], fp[m], j[m]
    ne_raw = plasma_density(fp)
    ne_cor = plasma_density(fp, te[j])
    C_raw, C_cor = ne_raw / pk[j], ne_cor / pk[j]
    ok = n.isfinite(C_cor)
    tp, fp, j, ne_raw, ne_cor, C_raw, C_cor = (a[ok] for a in (tp, fp, j, ne_raw, ne_cor, C_raw, C_cor))

    dig = read_digisonde(dig_file)
    dig = dig[(dig[:, 1] >= DIG_MIN_CS) & n.isfinite(dig[:, 2]) & n.isfinite(dig[:, 3])]
    nm_dig = 1.24e10 * dig[:, 2] ** 2
    jd = nearest(dig[:, 0], tc, MATCH_DIG_S)
    md = (jd >= 0) & n.isfinite(pk) & (n.abs(dig[n.maximum(jd, 0), 3] - hm) < DIG_MAX_DH_KM)
    td, C_dig = tc[md], nm_dig[jd[md]] / pk[md]

    lines = ["fits: %s (%d windows, %d with a peak); picks: %s (%d matched)"
             % (fits_dir, len(tc), n.isfinite(pk).sum(), cal, len(tp))]
    k_cor, C_pl, e_pl = constant("plasma lines, corrected (thermal, magnetic)", tp, C_cor, lines)
    k_raw, C_pl_raw, _ = constant("plasma lines, uncorrected (f/8.98)^2", tp, C_raw, lines)
    day = n.floor((tp - DAY0) / 86400).astype(int)
    for d in sorted(set(day)):
        constant("  ... corrected, %d April" % (7 + d), tp[day == d], C_cor[day == d], lines)
    ecl = (tp > ECLIPSE[0]) & (tp < ECLIPSE[1])
    if ecl.sum() > 2:
        constant("  ... corrected, during the eclipse", tp[ecl], C_cor[ecl], lines)
    k_dig, C_d, e_d = constant("digisonde, heights agree", td, C_dig, lines)
    dday = n.floor((td - DAY0) / 86400).astype(int)
    for d in sorted(set(dday)):
        constant("  ... digisonde, %d April" % (7 + d), td[dday == d], C_dig[dday == d], lines)
    lines.append("plasma lines (corrected) / digisonde: %.3f +- %.3f; uncorrected / digisonde: %.3f"
                 % (C_pl / C_d, C_pl / C_d * n.hypot(e_pl / C_pl, e_d / C_d), C_pl_raw / C_d))

    # the two absolute peak densities at the same times, without the radar's scale
    jp = nearest(dig[:, 0], tp, MATCH_DIG_S)
    mp = k_cor & (jp >= 0) & (n.abs(dig[n.maximum(jp, 0), 3] - hm[j]) < DIG_MAX_DH_KM)
    ratio = ne_cor[mp] / nm_dig[jp[mp]]
    if mp.sum() > 2:
        q = n.percentile(ratio, [16, 50, 84])
        lines.append("N_mF_2 plasma line (corrected) / digisonde, same time: median %.3f "
                     "(16-84 %%: %.3f-%.3f, n=%d)" % (q[1], q[0], q[2], mp.sum()))
    open(os.path.join(out, "summary.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))

    n.savetxt(os.path.join(out, "picks.txt"),
              n.column_stack([tp, fp / 1e6, te[j], ne_raw, ne_cor, pk[j], hm[j], C_raw, C_cor, k_cor]),
              fmt="%.10g", header="t_unix f_MHz Te_K ne_uncorrected ne_corrected isr_peak_rel "
              "isr_hmF2_km C_uncorrected C_corrected kept")
    n.savetxt(os.path.join(out, "digisonde.txt"),
              n.column_stack([td, pk[md], hm[md], dig[jd[md], 2], dig[jd[md], 3], C_dig, k_dig]),
              fmt="%.10g", header="t_unix isr_peak_rel isr_hmF2_km dig_foF2_MHz dig_hmF2_km C kept")

    hrs = lambda t: (t - DAY0) / 3600.0
    fig, ax = plt.subplots(3, 1, figsize=(9, 9.5), sharex=True)
    ax[0].semilogy(hrs(tp[k_cor]), C_cor[k_cor], ".", ms=3, label="plasma line, corrected")
    ax[0].semilogy(hrs(tp[~k_cor]), C_cor[~k_cor], "x", ms=4, color="0.5", label="rejected")
    ax[0].semilogy(hrs(td[k_dig]), C_dig[k_dig], "+", ms=5, color="C2", label="digisonde")
    ax[0].axhline(C_pl, color="C0", label="plasma lines: %.3g" % C_pl)
    ax[0].axhline(C_pl_raw, color="C0", ls=":", label="uncorrected: %.3g" % C_pl_raw)
    ax[0].axhline(C_d, color="C2", ls="--", label="digisonde: %.3g" % C_d)
    ax[0].set_ylabel("C = $n_e$ (m$^{-3}$) / relative peak")
    ax[0].legend(fontsize=7, ncol=2)
    ax[1].semilogy(hrs(tp[k_cor]), ne_cor[k_cor], ".", ms=3, label="plasma line, corrected")
    ax[1].semilogy(hrs(td[k_dig]), nm_dig[jd[md]][k_dig], "+", ms=5, color="C2", label="digisonde $N_mF_2$")
    # no line across the night between the two days
    gap = n.concatenate([[False], n.diff(tc) > 3600])
    ax[1].semilogy(hrs(tc), n.where(gap, n.nan, pk * C_pl), "-", lw=0.6, color="0.3",
                   label="radar peak $\\times$ C")
    ax[1].set_ylabel("$N_mF_2$ (m$^{-3}$)")
    ax[1].legend(fontsize=7)
    if mp.sum() > 2:
        ax[2].plot(hrs(tp[mp]), ratio, ".", ms=3)
    ax[2].axhline(1.0, color="0.3", lw=0.8)
    ax[2].set_ylabel("plasma line / digisonde $N_mF_2$")
    ax[2].set_xlabel("hours since 2024-04-07 00 UT")
    for a in ax:
        a.axvspan(hrs(ECLIPSE[0]), hrs(ECLIPSE[1]), color="orange", alpha=0.15)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "plasma_line_density_scale.png"), dpi=130)


if __name__ == "__main__":
    main()
