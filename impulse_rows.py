"""
Removal of the lagged products that impulsive interference touches.

Sparks on power-line hardware near the radar emit broadband impulses of about
a microsecond, locked to the 60 Hz mains (memo 21); in the day zenith-l also
carries a second, unlocked population of single-sample impulses (memos 26,
33). At 1 MHz sampling each is one or two samples tens to hundreds of times
the noise power, anywhere in the echo window.

The method (Juha's recommendation, Henrik 25 October 2026) leaves the voltage
as it is and removes from the lag profile inversion every lagged product
(row of m_j, with the same rows of the theory matrix A_j and of any satellite
column) whose samples touch an impulse. The rows that remain contain no
impulse power, so the least-squares estimate from them is unbiased; the price
is the removed rows' information. Earlier, impulse_blanking.BlankingReader
(branch main-before-review-2026-10-03) set the impulse samples to zero
instead, which leaves the zeroed samples in the products.

Detection is BlankingReader's test: in each pulse's echo window [gc,
last_echo) the power of every raw sample is compared with the running median
of the `window` samples around it, and samples above `thresh` times that
median are impulses, together with `guard` samples either side. The running
median follows the echo's range profile, so extended echoes, satellites
included, are not detected. For noise alone a sample exceeds 30 times the
median power in about one case in 2e9 (memo 21).

Which rows an impulse touches. The echo is low-pass filtered (outlier_lpi
fft_lpf, 2L taps centred on the sample, L = filter_len): the filtered sample
t contains the raw samples t-L+1 ... t+L, so an impulse at raw sample s
reaches the filtered samples s-L ... s+L-1. The lagged product of lag tau at
t is z(t) z*(t+tau), so it is touched when t or t+tau is. A row is the sum
of the products of `rg` consecutive t (simple_decimator), and it is removed
when any of them is touched. The clutter-subtracted products use the
difference of two pulses of the same code, so the impulses of both pulses
count; as the plain and the clutter-subtracted products share one design
matrix and one solve, a row is removed from both.
"""

import numpy as n
from scipy.ndimage import median_filter, binary_dilation

# detection, as impulse_blanking.BlankingReader (memo 21)
THRESH = 30.0
GUARD = 2
WINDOW = 65


def detect(z, a, b, thresh=THRESH, guard=GUARD, window=WINDOW):
    """
    Impulse samples of one pulse's raw voltage z (complex, not filtered):
    a boolean array of len(z), True at the samples of [a, b) whose power
    exceeds thresh times the running median of the `window` samples around
    it, and at `guard` samples either side of each.
    """
    hot = n.zeros(len(z), dtype=bool)
    b = min(b, len(z))
    if b - a < window:
        return hot
    p = n.abs(z[a:b]) ** 2
    # mirrored at the window's ends: BlankingReader's "nearest" repeats the
    # end sample half a window long, and when its power is low the median
    # there is too, which gave false impulses at the window's edges (in about
    # one noise-only pulse in five)
    h = p > thresh * median_filter(p, size=window, mode="mirror")
    if h.any() and guard:
        h = binary_dilation(h, iterations=guard)
    hot[a:b] = h
    return hot


def widen(hot, filter_len):
    """
    The filtered samples an impulse reaches: True at t when any raw impulse
    sample lies in t-L+1 ... t+L (L = filter_len), the span of outlier_lpi's
    fft_lpf (2L taps, centred).
    """
    L = int(filter_len)
    if not hot.any():
        return n.zeros(len(hot), dtype=bool)
    c = n.convolve(hot.astype(n.int64), n.ones(2 * L, dtype=n.int64))
    return c[L:L + len(hot)] > 0


def hit_rows(touched, lag, decim, m0, m1):
    """
    Rows m0 ... m1-1 of one lag whose lagged products touch an impulse:
    True where any product z(t) z*(t+lag) summed into the row has t or t+lag
    among the touched filtered samples. The products are summed into rows by
    decim (outlier_lpi.simple_decimator) as the measurement is.
    """
    N = len(touched)
    t = touched[0:(N - lag)] | touched[lag:N]
    return decim.decimate(t.astype(n.float64))[m0:m1] > 0


def masked_running_sigma(p, keep, width=10):
    """
    The per-row noise level of the variance weights, a running mean of |m|^2
    over `width` pulses (outlier_lpi: a circular moving average along the
    pulses, then the square root, shifted by width/2), taken over the kept
    rows only, so that the removed impulse power does not lower the weight of
    the clean rows near it. p: |m|^2, shape (pulses, rows); keep: bool, the
    same shape. Where none of the `width` pulses is kept the plain running
    mean is used (those rows are removed anyway, unless a neighbour's
    impulse caused it).
    """
    from numpy.fft import fft, ifft
    wf = n.repeat(1.0 / width, width)
    WF = fft(wf, p.shape[0])[:, None]
    k = keep.astype(n.float64)
    num = ifft(WF * fft(p * k, axis=0), axis=0).real
    den = ifft(WF * fft(k, axis=0), axis=0).real
    plain = ifft(WF * fft(p, axis=0), axis=0).real
    # den is a multiple of 1/width; below half of that no row was kept
    m = n.where(den > 0.5 / width, num / n.maximum(den, 1e-300), plain)
    return n.roll(n.sqrt(n.maximum(m, 0.0)), -(width // 2), axis=0)
