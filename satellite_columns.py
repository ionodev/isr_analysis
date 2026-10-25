#!/usr/bin/env python3
"""
Satellite columns for the lag profile inversion.

A satellite is a coherent point target. Within one pulse its echo is a
delayed, Doppler shifted copy of the transmitted pulse, so its contribution to
the lagged products of that pulse has a known shape, and one unknown complex
amplitude per lag. outlier_lpi.lpi_files() can add one such column to the
design matrix for every detected echo, beside the plasma columns and the range
independent background column, so the contaminated lagged products are
explained by the model rather than rejected.

The column is built the way the echo itself is processed: the inversion's own
transmit template (the pulse leaking into the echo channel, normalised to unit
energy) is delayed by the echo delay, Doppler shifted, zeroed outside the echo
window, low pass filtered with the echo's filter, and its lagged products are
decimated into range gates. Because the template is the leakage in the same
channel, the tx-h to echo channel delay does not appear here. It appears only
when a delay measured against tx-h, as the satellite catalogue's is, is
converted into this frame, which CatalogueDetections does.

The amplitude fitted to a column is the echo energy of that pulse, in the units
of the lagged products, times the phase exp(-j 2 pi df tau) left over from any
error df in the Doppler the column was built with. Its magnitude should not
change with lag, which is the check that the echo was a coherent target.

    satellite_template()   the synthetic, filtered echo of one detection
    template_lag_column()  its decimated lagged product at one lag
    refine_delay()         a catalogue delay refined below the catalogue's grid
    CatalogueDetections    detections from the single pulse catalogue
    TrackedDetections      the catalogue's detections linked into passes,
                           with columns for the pulses of a pass the
                           detector missed
    RandomControlColumns   columns at random places without an echo, the
                           pure-plasma control
"""

import numpy as n
from digital_rf import DigitalMetadataReader

from radar_timing import TMM as tmm
from tx_delay import fractional_shift


def satellite_template(z_tx, delay_samples, doppler_hz, lpf, gc, last_echo, sr=1e6):
    """
    The echo a point target would leave in this pulse after the processing the
    inversion applies to the echo.

    z_tx           the inversion's transmit template for this pulse: the
                   leakage in the echo channel, zero outside the transmit gate,
                   unit energy
    delay_samples  echo delay relative to that template, need not be an integer
    doppler_hz     Doppler shift; the echo is taken to rotate as exp(+j 2 pi f t)
    lpf            the filter applied to the echo (outlier_lpi.fft_lpf)
    gc, last_echo  the echo window. The echo is zeroed outside it before
                   filtering, so the template is too.
    """
    t = n.arange(len(z_tx))
    s = fractional_shift(n.array(z_tx, dtype=n.complex128), delay_samples)
    s = s * n.exp(2j * n.pi * doppler_hz * t / sr)
    s[0:gc] = 0.0
    s[last_echo:len(s)] = 0.0
    return lpf.lpf(s)


def average_templates(reader, keys, sid, channel, z_dc, tmm, keep, n_len=10000, n_iter=2):
    """
    The transmit template of each pulse code, averaged over the pulses of a
    period.

    The leakage of one pulse into the echo channel is only ~20 dB above the
    receiver noise, so a template taken from a single pulse cannot describe an
    echo much brighter than that: the template's own noise leaves a residual of
    the echo ~20 dB down. Averaging the pulses of one code lowers that noise by
    the number of pulses. Each pulse's leakage is first normalised to unit
    energy and rotated onto a common phase, as the transmitter phase need not
    be the same from pulse to pulse; neither matters to a lagged product.

    keys    the pulses of the period
    keep    keep(key) is True for a pulse the inversion uses
    returns {sweepid: (template, number of pulses, scatter)}. scatter is the
            mean energy of a pulse's difference from the average, relative to
            the average: the template noise if the transmitted pulses are all
            alike, more if they are not.
    """
    groups = {}
    for key in keys:
        s = sid[key]
        if s not in tmm or not keep(key):
            continue
        tx0, tx1 = tmm[s]["tx0"], tmm[s]["tx1"]
        v = reader.read_vector_1d(key, tx1, channel).astype("c8", casting="unsafe", copy=False) - z_dc
        u = n.zeros(n_len, dtype=n.complex128)
        u[tx0:tx1] = v[tx0:tx1]
        e = n.sqrt(n.sum(n.abs(u)**2))
        if e > 0 and n.isfinite(e):
            groups.setdefault(s, []).append(u / e)
    out = {}
    for s, us in groups.items():
        U = n.array(us)
        ref = U[0]
        for _ in range(n_iter):
            ph = n.angle(U.dot(n.conj(ref)))
            aligned = U * n.exp(-1j * ph)[:, None]
            ref = aligned.mean(axis=0)
        scatter = n.mean(n.sum(n.abs(aligned - ref)**2, axis=1)) / n.sum(n.abs(ref)**2)
        out[s] = (ref / n.sqrt(n.sum(n.abs(ref)**2)), len(us), float(scatter))
    return out


def template_lag_column(s, lag, decim, m0, m1):
    """
    Decimated lagged product of a satellite template at one lag, over the
    measurement rows m0:m1. Formed exactly as the echo's lagged products are.
    """
    return decim.decimate(s[0:(len(s) - lag)] * n.conj(s[lag:len(s)]))[m0:m1]


def refine_delay(echo, tx, raw_delay, doppler_hz, tx0, tx1, search=12, sr=1e6):
    """
    Refine a catalogue delay below its grid of 8 samples.

    The catalogue matched filters a copy of the echo summed over blocks of 8
    samples, so its raw_delay_sample is on that grid, +- 4 samples of the true
    delay. Here the full rate echo is matched filtered with the pulse's own
    tx-h samples, with the catalogue's Doppler removed, at every integer delay
    within +- search samples of the catalogue's, and a parabola through the
    peak power and its two neighbours gives the sub-sample delay.

    echo, tx     the pulse's samples on the echo channel and on tx-h, from the
                 pulse start
    raw_delay    the catalogue's raw_delay_sample: where in echo the copy of
                 tx[tx0:tx1] starts
    returns      (refined raw delay, float, in the catalogue's frame; peak
                 matched filter power over the power at the catalogue's delay)
    """
    t = n.asarray(tx[tx0:tx1], dtype=n.complex128)
    t = t - n.mean(tx[0:tx0])
    L = len(t)
    ramp = n.exp(-2j * n.pi * doppler_hz * n.arange(L) / sr)
    w = n.conj(t) * ramp
    s0 = max(0, int(raw_delay) - search)
    s1 = min(len(echo) - L, int(raw_delay) + search)
    shifts = n.arange(s0, s1 + 1)
    if len(shifts) < 3:
        return float(raw_delay), 1.0
    seg = n.lib.stride_tricks.sliding_window_view(n.asarray(echo, dtype=n.complex128), L)[s0:s1 + 1]
    p = n.abs(seg @ w)**2
    k = int(n.argmax(p))
    frac = 0.0
    if 0 < k < len(p) - 1:
        den = p[k - 1] - 2 * p[k] + p[k + 1]
        if den < 0:
            frac = float(n.clip(0.5 * (p[k - 1] - p[k + 1]) / den, -0.5, 0.5))
    i_cat = int(raw_delay) - s0
    gain = float(p[k] / p[i_cat]) if 0 <= i_cat < len(p) and p[i_cat] > 0 else 1.0
    return float(shifts[k] + frac), gain


class CatalogueDetections:
    """
    Echoes from the single pulse satellite catalogue
    (single_pulse_satellite_cfar.py) for one channel, with their delays
    converted into the frame of the lag profile inversion.

    The catalogue measures delay against tx-h. Its corrected_delay_samples is
    raw_delay_sample - tx0 - 11, the 11 being a fixed receiver delay. The
    inversion measures delay against the pulse leaking into the echo channel,
    so the conversion here subtracts the measured channel delay (tx_delay.py)
    instead of the fixed one.

    The catalogue searches a grid decimated by eight, so its delays are known
    to +- 4 samples, +- 0.6 km in range. Given a reader, each delay is refined
    below that grid from the raw samples (refine_delay); without one, they are
    used as the catalogue gives them.

    for_period() returns {pulse key: [(delay_samples, doppler_hz), ...]}, the
    form outlier_lpi.lpi_files() takes. The optional windows restrict it to one
    pass: t_window in unix seconds, range_window in km as the catalogue states
    range, min_snr_db on the catalogue's matched filter SNR.
    """

    def __init__(self, metadata_dir, channel, channel_delay_samples,
                 t_window=None, range_window=None, min_snr_db=None, refine_reader=None):
        self.reader = DigitalMetadataReader(metadata_dir)
        # an rf reader (read_vector_1d) to refine the delays with, or None
        self.refine_reader = refine_reader
        self.channel = channel
        self.channel_delay_samples = float(channel_delay_samples)
        self.t_window = t_window
        self.range_window = range_window
        self.min_snr_db = min_snr_db

    def for_period(self, i0, i1):
        return {k: [(d[0], d[1]) for d in v] for k, v in self._detections(i0, i1).items()}

    def _detections(self, i0, i1):
        if self.t_window is not None:
            i0 = max(i0, int(self.t_window[0] * 1e6))
            i1 = min(i1, int(self.t_window[1] * 1e6))
            if i1 <= i0:
                return {}
        out = {}
        for key, rec in self.reader.read(i0, i1).items():
            ch = rec["channel"]
            if isinstance(ch, bytes):
                ch = ch.decode()
            if ch != self.channel:
                continue
            for raw, sweep_id, doppler_hz, snr_db, range_km in zip(
                    rec["raw_delay_sample"], rec["sweep_id"], rec["doppler_hz"],
                    rec["snr_db"], rec["range_km"]):
                if int(sweep_id) not in tmm:
                    continue
                if self.min_snr_db is not None and snr_db < self.min_snr_db:
                    continue
                if self.range_window is not None and not (
                        self.range_window[0] <= range_km <= self.range_window[1]):
                    continue
                raw = float(raw)
                if self.refine_reader is not None:
                    raw = self.refine(int(key), int(sweep_id), raw, float(doppler_hz))
                delay = raw - tmm[int(sweep_id)]["tx0"] - self.channel_delay_samples
                out.setdefault(int(key), []).append((delay, float(doppler_hz), float(snr_db)))
        return out

    def rows(self, i0, i1):
        """(pulse key, delay, Doppler, SNR dB) of every echo in [i0, i1)."""
        return [(k, d[0], d[1], d[2]) for k, v in self._detections(i0, i1).items() for d in v]

    def refine(self, key, sweep_id, raw, doppler_hz):
        m = tmm[sweep_id]
        n_read = int(raw) + (m["tx1"] - m["tx0"]) + 32
        try:
            echo = self.reader_rf(key, n_read, self.channel)
            tx = self.reader_rf(key, m["tx1"], "tx-h")
        except Exception:
            return raw
        return refine_delay(echo, tx, raw, doppler_hz, m["tx0"], m["tx1"])[0]

    def reader_rf(self, key, length, channel):
        return self.refine_reader.read_vector_1d(key, length, channel)


class TrackedDetections:
    """
    Detections linked into passes, with a column for every pulse of a pass,
    detected or not.

    The single-pulse detector finds a pulse's echo only above about 19 dB
    matched-filter SNR (end-to-end injection, 13 October), so the faint
    pulses at the start and end of every pass get no column and bias the
    plasma ACF. A satellite moves smoothly, so its delay and Doppler on the
    pulses in between and around the detections follow from the detected
    ones:

      1. the detections of [i0 - margin, i1 + margin), in time order, are
         linked into passes: a detection joins the pass whose predicted
         delay (from its last detection, moved by its Doppler) is within
         link_delay samples, whose Doppler is within link_doppler_hz, and
         whose last detection is at most link_gap_s earlier;
      2. a pass with at least min_detections detections (default one: at the
         detector's false-alarm probability of 1e-10 a single detection is a
         satellite, and a satellite stays in the beam for seconds, Memo 13)
         gets delay(t) as a
         quadratic, Doppler(t) as a line, and SNR(t) in dB as a parabola (a
         beam crossing), fitted to its detections;
      3. its span is extended beyond the first and last detection to where
         the fitted SNR falls to extend_snr_db, at most max_extend_s; when
         the SNR cannot be fitted (fewer than five detections, or a parabola
         that does not open downward), by default_extend_s (2 s: passes last
         3 to over 10 s, Memo 13);
      4. every pulse in the span gets a column at the fitted delay and
         Doppler, except detected pulses, which keep their own; columns of
         one pulse closer than merge_delay samples and merge_doppler_hz are
         merged (one satellite reached twice, by two fragments of a pass;
         two such columns would be indistinguishable).

    Too long a span costs a little precision (a column without an echo fits
    about zero); too short a span leaves the faint pulses' bias, so the
    defaults lean long.

    base       an object with rows(i0, i1) -> [(key, delay, Doppler, SNR dB)]
               (CatalogueDetections, or a test's own)
    pulse_keys a function (i0, i1) -> the pulse keys a column may be put on
    """

    def __init__(self, base, pulse_keys, link_delay=10.0, link_doppler_hz=1500.0, link_gap_s=1.5,
                 min_detections=1, extend_snr_db=0.0, max_extend_s=3.0, default_extend_s=2.0,
                 margin_s=10.0, merge_delay=2.0, merge_doppler_hz=500.0, sr=1e6):
        self.base = base
        self.pulse_keys = pulse_keys
        self.link_delay = link_delay
        self.link_doppler_hz = link_doppler_hz
        self.link_gap_s = link_gap_s
        self.min_detections = min_detections
        self.extend_snr_db = extend_snr_db
        self.max_extend_s = max_extend_s
        self.default_extend_s = default_extend_s
        self.margin = int(margin_s * sr)
        self.merge_delay = merge_delay
        self.merge_doppler_hz = merge_doppler_hz
        self.sr = sr
        self.last_passes = []

    def passes(self, rows):
        """Link detections (key, delay, Doppler, SNR) into passes: lists of rows."""
        lam = 299792458.0 / 440.2e6
        tracks = []
        for r in sorted(rows, key=lambda x: x[0]):
            k, d, f, _ = r
            best, best_cost = None, None
            for tr in tracks:
                k0, d0, f0, _ = tr[-1]
                dt = (k - k0) / self.sr
                if dt > self.link_gap_s or dt < 0:
                    continue
                # delay rate in samples per second from the Doppler: -f lambda / c
                pred = d0 - f0 * lam / 299792458.0 * self.sr * dt
                dd, df = abs(d - pred), abs(f - f0)
                if dd <= self.link_delay and df <= self.link_doppler_hz:
                    cost = dd / self.link_delay + df / self.link_doppler_hz
                    if best is None or cost < best_cost:
                        best, best_cost = tr, cost
            if best is None:
                tracks.append([r])
            else:
                best.append(r)
        return tracks

    def span(self, tr):
        """Fitted delay(t), Doppler(t) and the time span (keys) of a pass."""
        k = n.array([r[0] for r in tr], dtype=float)
        t = (k - k[0]) / self.sr
        d = n.array([r[1] for r in tr])
        f = n.array([r[2] for r in tr])
        snr = n.array([r[3] for r in tr])
        deg_d = 2 if len(tr) >= 5 and t[-1] > 0.5 else 1
        pd = n.polyfit(t, d, deg_d) if t[-1] > 0 else n.array([d.mean()])
        pf = n.polyfit(t, f, 1) if t[-1] > 0 else n.array([f.mean()])
        ext0 = ext1 = self.default_extend_s
        if len(tr) >= 5 and t[-1] > 0.5:
            ps = n.polyfit(t, snr, 2)
            if ps[0] < 0:
                # times where the parabola equals extend_snr_db
                roots = n.roots([ps[0], ps[1], ps[2] - self.extend_snr_db])
                roots = n.real(roots[n.isreal(roots)])
                if len(roots) == 2:
                    ext0 = float(n.clip(0.0 - roots.min(), 0.0, self.max_extend_s))
                    ext1 = float(n.clip(roots.max() - t[-1], 0.0, self.max_extend_s))
        k0 = k[0] - ext0 * self.sr
        k1 = k[-1] + ext1 * self.sr
        return pd, pf, k[0], k0, k1

    def for_period(self, i0, i1):
        rows = self.base.rows(i0 - self.margin, i1 + self.margin)
        detected = {}
        for r in rows:
            detected.setdefault(int(r[0]), []).append((r[1], r[2]))
        out = {k: list(v) for k, v in detected.items() if i0 <= k < i1}
        keys = n.array(sorted(self.pulse_keys(i0, i1)), dtype=n.int64)
        self.last_passes = []
        for tr in self.passes(rows):
            if len(tr) < self.min_detections:
                continue
            pd, pf, kref, k0, k1 = self.span(tr)
            own = set(int(r[0]) for r in tr)
            sel = keys[(keys >= k0) & (keys <= k1)]
            added = 0
            for k in sel:
                if int(k) in own:
                    continue
                t = (k - kref) / self.sr
                out.setdefault(int(k), []).append((float(n.polyval(pd, t)), float(n.polyval(pf, t))))
                added += 1
            self.last_passes.append(dict(detections=len(tr), added=added, start=int(k0), stop=int(k1)))
        for k in out:
            kept = []
            for d, f in out[k]:
                if all(abs(d - d2) >= self.merge_delay or abs(f - f2) >= self.merge_doppler_hz for d2, f2 in kept):
                    kept.append((d, f))
            out[k] = kept
        return out


class RandomControlColumns:
    """
    The pure-plasma control (Memo 46, option F): satellite columns at random
    pulses, delays and Dopplers where there is no echo, added to those of a
    base detector.

    The claim tested is that a column where no satellite is (a false
    detection, or the extrapolated end of a pass) costs precision and not
    bias: the inversion fits its amplitude to about zero, takes from the
    plasma the information of that pulse's rows at that range, and says so
    in the error bar. Columns on pure plasma show whether that holds.

    Per period, per_period columns are drawn uniformly over the period's
    pulses (pulse_keys), delays in [delay_min, delay_max] samples (the round
    trip in microseconds: 150-800 km is 1000-5337) and Dopplers in
    [-doppler_max_hz, doppler_max_hz] (the detector's search band, where a
    false alarm would fall). A draw is kept only if it stays clear of every
    column the base gives, on its own pulse and on the pulses within
    `neighbours` of it (the clutter-subtracted products carry the echoes of
    the pulse a pulse is differenced against), by at least avoid_delay
    samples (the pulse length) in delay and avoid_doppler_hz in Doppler.
    The draws are seeded by the seed and the period, so a run and its
    analysis see the same columns.

    base        an object with for_period(i0, i1), or None for no others
    pulse_keys  a function (i0, i1) -> the pulse keys a column may be put on
    """

    def __init__(self, base, pulse_keys, per_period, seed=1, delay_min=1000.0, delay_max=5337.0,
                 doppler_max_hz=60e3, avoid_delay=480.0, avoid_doppler_hz=2000.0, neighbours=3):
        self.base = base
        self.pulse_keys = pulse_keys
        self.per_period = int(per_period)
        self.seed = int(seed)
        self.delay_min, self.delay_max = float(delay_min), float(delay_max)
        self.doppler_max_hz = float(doppler_max_hz)
        self.avoid_delay = float(avoid_delay)
        self.avoid_doppler_hz = float(avoid_doppler_hz)
        self.neighbours = int(neighbours)

    def random_columns(self, i0, i1, base_cols):
        """{key: [(delay, Doppler), ...]} of the random columns of [i0, i1)."""
        keys = n.array(sorted(self.pulse_keys(i0, i1)), dtype=n.int64)
        if len(keys) == 0 or self.per_period <= 0:
            return {}
        rng = n.random.default_rng([self.seed, int(i0 // 1000000)])
        out = {}
        n_made, n_tried = 0, 0
        while n_made < self.per_period and n_tried < 50 * self.per_period:
            n_tried += 1
            pi = int(rng.integers(len(keys)))
            d = float(rng.uniform(self.delay_min, self.delay_max))
            f = float(rng.uniform(-self.doppler_max_hz, self.doppler_max_hz))
            near = keys[max(pi - self.neighbours, 0):pi + self.neighbours + 1]
            clash = False
            for k in near:
                for d2, f2 in base_cols.get(int(k), []):
                    if abs(d - d2) < self.avoid_delay or abs(f - f2) < self.avoid_doppler_hz:
                        clash = True
                for d2, f2 in out.get(int(k), []):
                    # two random columns of one pulse kept apart as well
                    if int(k) == int(keys[pi]) and abs(d - d2) < self.avoid_delay:
                        clash = True
            if clash:
                continue
            out.setdefault(int(keys[pi]), []).append((d, f))
            n_made += 1
        return out

    def for_period(self, i0, i1):
        base_cols = self.base.for_period(i0, i1) if self.base is not None else {}
        rand = self.random_columns(i0, i1, base_cols)
        out = {k: list(v) for k, v in base_cols.items()}
        for k, v in rand.items():
            out.setdefault(k, []).extend(v)
        return out


def zenith_pulse_keys(data_dir):
    """A pulse_keys function for TrackedDetections: the pulses of [i0, i1)
    whose code is in the timing table."""
    reader = DigitalMetadataReader("%s/metadata/id_metadata" % data_dir)

    def keys(i0, i1):
        sid = reader.read(int(i0), int(i1), "sweepid")
        return [int(k) for k, v in sid.items() if int(v) in tmm]
    return keys
