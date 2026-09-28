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
    CatalogueDetections    detections from the single pulse catalogue
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
    to +- 4 samples, +- 0.6 km in range. Nothing here refines them.

    for_period() returns {pulse key: [(delay_samples, doppler_hz), ...]}, the
    form outlier_lpi.lpi_files() takes. The optional windows restrict it to one
    pass: t_window in unix seconds, range_window in km as the catalogue states
    range, min_snr_db on the catalogue's matched filter SNR.
    """

    def __init__(self, metadata_dir, channel, channel_delay_samples,
                 t_window=None, range_window=None, min_snr_db=None):
        self.reader = DigitalMetadataReader(metadata_dir)
        self.channel = channel
        self.channel_delay_samples = float(channel_delay_samples)
        self.t_window = t_window
        self.range_window = range_window
        self.min_snr_db = min_snr_db

    def for_period(self, i0, i1):
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
                delay = float(raw) - tmm[int(sweep_id)]["tx0"] - self.channel_delay_samples
                out.setdefault(int(key), []).append((delay, float(doppler_hz)))
        return out
