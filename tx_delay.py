#!/usr/bin/env python3
"""
Signal processing channel delay between the transmit sample channel (tx-h) and
the echo channels (zenith-l, misa-l).

tx-h and the echo channels are digitised by separate receiver chains, each with
its own filter and decimation group delay, so a transmit waveform read from
tx-h is not aligned in time with the same pulse as it appears in an echo
channel. The offset is around 10-12 microseconds, far too long to be cable
length, and it differs per channel and per experiment.

avg_range_doppler_spec.py used to correct this with a hardcoded n.roll(z_tx,11),
eyeballed from the leakthrough and quantised to whole samples (1 us, 150 m in
range). The functions here replace that with a measurement:

    estimate_channel_delay()  matched filters tx-h against the transmit pulse
                              leaking into the echo channel, over many pulses
    fractional_shift()        shifts a waveform by a non integer number of
                              samples, so the correction is not quantised

The estimate deliberately uses only the phase coded pulses. An uncoded long
pulse autocorrelates to a broad triangle with no sharp peak to locate, and
gives a delay biased by ~0.2 us relative to the coded pulses in the same
recording. The delay is a property of the receiver chain, not of the waveform,
so the coded pulse estimate is the one to apply to every mode.
"""

import numpy as n
from digital_rf import DigitalRFReader, DigitalMetadataReader

import millstone_radar_state as mrs

# channel carrying the transmit samples
TX_CHANNEL = "tx-h"

# transmit pulse timing in samples from the start of the interpulse period,
# mirrors the tables in outlier_lpi.py and avg_range_doppler_spec.py
tmm = {}
tmm[300] = {"noise0": 7800, "noise1": 8371, "tx0": 76, "tx1": 645, "gc": 1000, "last_echo": 7700, "e_gc": 800}
for i in range(1, 33):
    tmm[i] = {"noise0": 8400, "noise1": 8850, "tx0": 76, "tx1": 624, "gc": 1000, "last_echo": 8200, "e_gc": 800}
tmm[800] = {"noise0": 30176, "noise1": 32033, "tx0": 69, "tx1": 2171, "gc": 3721, "last_echo": 30000, "e_gc": 3721}

# sweepids of the phase coded pulses, the ones with a sharp correlation peak
CODED_SWEEPIDS = list(range(1, 33))

# fallback if the delay cannot be measured, the value the code used to hardcode
DEFAULT_DELAY_US = 11.0


def matched_filter_delay(z_tx, z_echo, oversample=100, max_lag=100.0):
    """
    Matched filter of the isolated transmit pulses of two channels:

        xc(tau) = sum_t z_tx(t) conj(z_echo(t - tau))

    i.e. z_tx convolved with the time reversed complex conjugate of z_echo. The
    cross correlation spectrum is zero padded before the inverse transform,
    which band limited interpolates xc onto a 1/oversample sample grid, so the
    delay is resolved to a fraction of a sample.

    Returns (delay in samples, peak amplitude). A positive delay means the echo
    channel lags the transmit channel.
    """
    n_s = len(z_tx)
    nfft = int(2 ** n.ceil(n.log2(2 * n_s)))
    C = n.fft.fft(z_tx, nfft) * n.conj(n.fft.fft(z_echo, nfft))

    nover = nfft * oversample
    Ci = n.zeros(nover, dtype=n.complex128)
    Ci[0:(nfft // 2)] = C[0:(nfft // 2)]
    Ci[(nover - nfft // 2):] = C[(nfft // 2):]
    xc = n.fft.ifft(Ci) * oversample

    lag = n.arange(nover) / float(oversample)
    lag[lag > nfft / 2] -= nfft

    # only look for the peak at physically plausible lags
    ok = n.abs(lag) <= max_lag
    pi = n.argmax(n.abs(xc[ok]))

    # xc peaks at tau = -delay, see the definition above
    return -lag[ok][pi], n.abs(xc[ok][pi])


def fractional_shift(z, d):
    """
    Delay z by d samples, d need not be an integer.

    Applies a phase ramp in the frequency domain, which is band limited
    (sinc) interpolation of the samples. Apply this to the raw waveform, before
    any hard windowing: shifting an array that has been zeroed outside a window
    rings at the window edges.
    """
    if d == 0.0:
        return z
    N = len(z)
    f = n.fft.fftfreq(N)
    return n.array(n.fft.ifft(n.fft.fft(z) * n.exp(-2j * n.pi * f * d)), dtype=z.dtype)


def channel_pulse_ok(channel, key, tx_ant, rx_ant, zpm, mpm, min_tx_pwr):
    """Was this pulse transmitted and received on the antenna of this channel?"""
    if channel in ["zenith-l", "zenith-l2"]:
        return (tx_ant(key) <= -0.99) and (rx_ant(key) <= -0.99) and (zpm(key / 1e6) >= min_tx_pwr)
    if channel == "misa-l":
        return (tx_ant(key) >= 0.99) and (rx_ant(key) >= 0.99) and (mpm(key / 1e6) >= min_tx_pwr)
    return True


def estimate_channel_delay(dirname, channel, n_pulses=100, oversample=100,
                           max_lag_us=100.0, min_tx_pwr=400e3, t0_unix=None, verbose=True,
                           zpm=None, mpm=None, tx_ant=None, rx_ant=None):
    """
    Measure the delay of an echo channel relative to tx-h, in microseconds.

    Matched filters the isolated transmit pulse of up to n_pulses phase coded
    pulses and returns the median, which is robust against the occasional pulse
    ruined by interference. Returns (delay_us, spread_us, n_used); the delay is
    None when it could not be measured.

    The transmit power and antenna selection models are read from the metadata
    when not given. They take a while to build, so pass them in when the caller
    already has them.
    """
    d_il = DigitalRFReader("%s/rf_data/" % (dirname))
    if TX_CHANNEL not in d_il.get_channels() or channel not in d_il.get_channels():
        if verbose:
            print("cannot measure tx delay: need both %s and %s" % (TX_CHANNEL, channel))
        return None, None, 0

    sr = float(d_il.get_properties(channel)["samples_per_second"])

    id_read = DigitalMetadataReader("%s/metadata/id_metadata" % (dirname))
    idsr = int(id_read.get_samples_per_second())
    idb = id_read.get_bounds()
    i0 = idb[0] if t0_unix is None else int(t0_unix * idsr)

    if zpm is None or mpm is None:
        zpm, mpm = mrs.get_tx_power_model("%s/metadata/powermeter" % (dirname))
    if tx_ant is None or rx_ant is None:
        tx_ant, rx_ant = mrs.get_antenna_select("%s/metadata/antenna_control_metadata" % (dirname))

    # enough pulses to find n_pulses coded ones on the right antenna
    sid = id_read.read(i0, min(i0 + 20 * n_pulses * 10000, idb[1]), "sweepid")
    keys = sorted(sid.keys())

    def measure(keys, screen):
        out = []
        for key in keys:
            if len(out) >= n_pulses:
                break
            if sid[key] not in CODED_SWEEPIDS:
                continue
            if screen and not channel_pulse_ok(channel, key, tx_ant, rx_ant, zpm, mpm, min_tx_pwr):
                continue
            t = tmm[sid[key]]
            try:
                z_tx = d_il.read_vector_c81d(int(key) + t["tx0"], t["tx1"] - t["tx0"], TX_CHANNEL)
                z_echo = d_il.read_vector_c81d(int(key) + t["tx0"], t["tx1"] - t["tx0"], channel)
            except Exception:
                continue
            d_samples, amp = matched_filter_delay(z_tx, z_echo, oversample=oversample,
                                                  max_lag=max_lag_us * sr / 1e6)
            out.append(d_samples / sr * 1e6)
        return out

    # prefer pulses transmitted on this channel's own antenna, they have the
    # strongest leakthrough
    delays = measure(keys, screen=True)

    # the antenna control metadata and the power meter disagree in some
    # recordings, which rejects every pulse. the delay is a property of the
    # receiver chain, not of which antenna transmitted, so fall back to any
    # coded pulse: cross antenna leakthrough is weaker but still well above
    # the noise.
    if len(delays) < 10:
        fallback = measure(keys, screen=False)
        if len(fallback) > len(delays):
            if verbose:
                print("%s: only %d pulses pass the antenna and power check, "
                      "using %d unscreened coded pulses instead"
                      % (channel, len(delays), len(fallback)))
            delays = fallback

    if len(delays) == 0:
        if verbose:
            print("cannot measure tx delay for %s: no usable coded pulses found" % (channel))
        return None, None, 0

    delays = n.array(delays)
    delay_us = float(n.median(delays))
    # median absolute deviation, scaled to be comparable to a standard deviation
    spread_us = float(1.4826 * n.median(n.abs(delays - delay_us)))

    # a pulse ruined by interference can correlate anywhere within max_lag, so
    # drop the outliers before reporting
    if spread_us > 0:
        ok = n.abs(delays - delay_us) < 5.0 * spread_us
        if n.sum(ok) >= 10:
            delays = delays[ok]
            delay_us = float(n.median(delays))
            spread_us = float(1.4826 * n.median(n.abs(delays - delay_us)))

    if verbose:
        print("%s delay relative to %s: %1.3f us (spread %1.3f us, %d coded pulses)"
              % (channel, TX_CHANNEL, delay_us, spread_us, len(delays)))
    return delay_us, spread_us, len(delays)
