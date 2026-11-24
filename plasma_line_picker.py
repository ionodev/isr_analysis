#!/usr/bin/env python3
"""
Automatic picker for the plasma-line ridge: the daily check of the density
scale (Henrik, 12 October: the digisonde sets the scale, the plasma lines
check it).

plasma_line_clicker.py shows, per minute and antenna, the plasma-line
spectrum's maximum over the gates at 185-560 km, and an operator clicks the
ridge. This module builds the same spectra (in parallel, cached), finds the
ridge itself, and writes the same file the clicker writes, so that
estimate_magic_constant.py reads it unchanged:

    <calname>.h5   plasma_frequency  (N, 2): unix time, frequency in Hz

The picker follows the ridge through time, on one side of zero at a time
(in this recording the downshifted line is clear on zenith and the
upshifted side carries discrete interference, so the two sides are not
alike):
  1. each minute's spectrum in units of its own noise (median removed,
     divided by 1.4826 times the median absolute deviation), clipped at
     clip_sigma so that a burst of interference cannot dominate, averaged in
     frequency bins of bin_hz, inside the search band;
     Each bin's median over the whole recording is then subtracted: a
     narrowband interference line stays at one frequency for hours (one at
     -6.25 MHz all night), the plasma line moves;
  2. the path through time with the highest summed score whose frequency
     moves by at most max_step_hz per minute (dynamic programming, as a
     Viterbi search): the ridge changes slowly, interference does not follow
     a smooth path;
  3. a minute's pick is kept where the path's score there is at least
     min_snr (the line is absent at night, where the path wanders in noise);
  4. only minutes with the sun less than min_sun_elev_deg below the horizon
     at Millstone are kept: the daytime plasma line is enhanced by
     photoelectrons and needs a sunlit F region (Memo 11), and at night
     interference outscores it;
  5. both sides are searched; the side with more kept picks is used.

The frequency is that of the line, which the clicker also records; the
conversion to density, n_e = (f/8.98)^2, and its caveats (no thermal or
Debye-length correction, Memo 11) are estimate_magic_constant.py's.

Usage:
    python3 plasma_line_picker.py <data_dir> <out_dir>
"""

import os
import sys
from multiprocessing import Pool

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import h5py
import numpy as n
from digital_rf import DigitalMetadataReader

from plasma_line_clicker import gate_window, PLASMA_LINE_RANGES_KM, MIN_GATES

FREQS = n.fft.fftshift(n.fft.fftfreq(12000, d=1 / 25e6))
# search band for |f|: above the region near zero offset, below the edges of
# the 25 MHz band, which carry interference above about 10.5 MHz
BAND_HZ = (2.0e6, 10.5e6)
BIN_HZ = 25e3
CLIP_SIGMA = 10.0
MAX_STEP_HZ = 0.15e6
MIN_SNR = 4.0
MIN_SUN_ELEV_DEG = -10.0
RADAR_LAT, RADAR_LON = 42.6193, -71.4912


def solar_elevation_deg(t_unix, lat=RADAR_LAT, lon=RADAR_LON):
    """Solar elevation from the standard low-precision formulas (about
    0.1 degree), enough to tell day from night."""
    d = n.asarray(t_unix, dtype=float) / 86400.0 + 2440587.5 - 2451545.0
    g = n.radians((357.529 + 0.98560028 * d) % 360)
    q = (280.459 + 0.98564736 * d) % 360
    lam = n.radians(q + 1.915 * n.sin(g) + 0.020 * n.sin(2 * g))
    eps = n.radians(23.439 - 0.00000036 * d)
    dec = n.arcsin(n.sin(eps) * n.sin(lam))
    ra = n.arctan2(n.cos(eps) * n.sin(lam), n.cos(lam))
    gmst = (18.697374558 + 24.06570982441908 * d) % 24
    ha = n.radians(gmst * 15 + lon) - ra
    la = n.radians(lat)
    return n.degrees(n.arcsin(n.sin(la) * n.sin(dec) + n.cos(la) * n.cos(dec) * n.cos(ha)))


def minute(args):
    """Both antennas' spectra for one minute, as plasma_line_clicker.build_spectra."""
    dirname, t0 = args
    plmd = DigitalMetadataReader("%s/metadata/integrated_plasma_line_metadata_hires/" % dirname)
    try:
        d = plmd.read_flatdict(t0, t0 + 60)
    except Exception:
        return None
    if "spec" not in d:
        return None
    out = {}
    for ant in ("zenith", "misa"):
        s, t, k = n.zeros(len(FREQS)), 0.0, 0
        for si in range(len(d["spec"])):
            is_misa = d["antenna"][si] == "MISA"
            if is_misa != (ant == "misa"):
                continue
            g = gate_window(d["ranges"][si], PLASMA_LINE_RANGES_KM)
            if (g.stop - g.start) < MIN_GATES:
                continue
            x = d["spec"][si][:, g]
            s += n.nanmax(x - n.nanmedian(x, axis=0), axis=1)
            t += d["t1"][si]
            k += 1
        if k > 0:
            s = s / k
            out[ant] = (t / k, s - n.nanmedian(s))
    return out


def build_spectra(dirname, cache, n_proc=48):
    """Per minute and antenna, the clicker's spectra; cached in an npz file."""
    if os.path.exists(cache):
        c = n.load(cache)
        return {a: (c[a + "_t"], c[a + "_s"]) for a in ("zenith", "misa")}
    plmd = DigitalMetadataReader("%s/metadata/integrated_plasma_line_metadata_hires/" % dirname)
    b = plmd.get_bounds()
    n_min = int((b[1] - b[0]) / plmd.get_file_cadence_secs())
    with Pool(n_proc) as p:
        res = p.map(minute, [(dirname, b[0] + i * 60) for i in range(n_min)], chunksize=4)
    S = {}
    for a in ("zenith", "misa"):
        rows = [r[a] for r in res if r is not None and a in r]
        t = n.array([x[0] for x in rows])
        s = n.vstack([x[1] for x in rows]) if rows else n.zeros((0, len(FREQS)))
        S[a] = (t, s)
    n.savez(cache, **{a + "_t": S[a][0] for a in S}, **{a + "_s": S[a][1] for a in S})
    return S


def side_scores(s, side, band=BAND_HZ, bin_hz=BIN_HZ, clip=CLIP_SIGMA):
    """Per minute, the clipped noise-normalised score in frequency bins of one
    side (side=+1 or -1), ordered by |f|. Returns (bin centres in Hz, score)."""
    sel = n.where((side * FREQS >= band[0]) & (side * FREQS <= band[1]))[0]
    sel = sel[n.argsort(n.abs(FREQS[sel]))]
    k = max(1, int(round(bin_hz / (FREQS[1] - FREQS[0]))))
    nb = len(sel) // k
    sel = sel[:nb * k]
    f = FREQS[sel].reshape(nb, k).mean(axis=1)
    sc = n.zeros((len(s), nb))
    for i, row in enumerate(s):
        med = n.nanmedian(row)
        sig = 1.4826 * n.nanmedian(n.abs(row - med))
        if not n.isfinite(sig) or sig <= 0:
            continue
        z = n.clip((row[sel] - med) / sig, -clip, clip)
        # the mean of k bins of unit noise has standard deviation 1/sqrt(k)
        sc[i] = z.reshape(nb, k).mean(axis=1) * n.sqrt(k)
    sc = n.nan_to_num(sc)
    # stationary lines: each bin's median over the recording
    sc = sc - n.median(sc, axis=0)
    return f, sc


def ridge(t, sc, f, max_step_hz=MAX_STEP_HZ):
    """The path of highest summed score whose frequency moves by at most
    max_step_hz per minute of elapsed time (Viterbi). Returns bin indices."""
    T, B = sc.shape
    df = abs(f[1] - f[0])
    acc = sc[0].copy()
    back = n.zeros((T, B), dtype=n.int32)
    for i in range(1, T):
        w = int(max(1, round(max_step_hz * max(1.0, (t[i] - t[i - 1]) / 60.0) / df)))
        # best predecessor within +-w bins: a running maximum by shifts
        best = acc.copy()
        arg = n.arange(B)
        for d in range(1, w + 1):
            for sh in (d, -d):
                cand = n.full(B, -n.inf)
                idx = n.arange(B) - sh
                ok = (idx >= 0) & (idx < B)
                cand[ok] = acc[idx[ok]]
                better = cand > best
                best[better] = cand[better]
                arg[better] = idx[better]
        back[i] = arg
        acc = best + sc[i]
    path = n.zeros(T, dtype=int)
    path[-1] = int(n.argmax(acc))
    for i in range(T - 1, 0, -1):
        path[i - 1] = back[i, path[i]]
    return path


def pick(t, s, min_snr=MIN_SNR):
    """Ridge picks on the better side: arrays (time, frequency in Hz, score)."""
    best = None
    for side in (-1, +1):
        f, sc = side_scores(s, side)
        path = ridge(t, sc, f)
        score = sc[n.arange(len(t)), path]
        keep = (score >= min_snr) & (solar_elevation_deg(t) >= MIN_SUN_ELEV_DEG)
        res = (t[keep], f[path][keep], score[keep])
        if best is None or keep.sum() > len(best[0]):
            best = res
    return best


def write_cal(out_dir, calname, t, f, score):
    """The file plasma_line_clicker.py writes, plus the scores and the method."""
    with h5py.File(os.path.join(out_dir, calname + ".h5"), "w") as h:
        h["plasma_frequency"] = n.column_stack([t, f])
        h["score_sigma"] = score
        h["method"] = "plasma_line_picker (automatic)"


def main():
    dirname, out_dir = sys.argv[1:3]
    os.makedirs(out_dir, exist_ok=True)
    S = build_spectra(dirname, os.path.join(out_dir, "plasma_line_spectra.npz"))
    for ant, cal in (("zenith", "zenithcal"), ("misa", "misacal")):
        t, s = S[ant]
        if len(t) == 0:
            continue
        pt, pf, ps = pick(t, s)
        write_cal(out_dir, cal, pt, pf, ps)
        print("%s: %d minutes, %d picks" % (ant, len(t), len(pt)))


if __name__ == "__main__":
    main()
