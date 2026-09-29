"""
Blanking of impulsive interference in the raw voltage.

Sparks on power-line hardware near the radar emit broadband impulses of about
a microsecond, locked to the 60 Hz mains (memo 21).  At 1 MHz sampling each
is one or two samples tens to hundreds of times the noise power.  They fall
anywhere in the interpulse period, so nothing tied to the radar's timing
removes them.

BlankingReader wraps a reader with read_vector_1d (RawReader or
DigitalRFReader).  In every vector it returns it compares each sample's power
with the running median of the surrounding `window` samples and sets samples
above `thresh` times that median to zero, together with `guard` samples on
either side.  The running median follows the echo's range profile, so extended
echoes, satellites included, are not blanked.  Only the echo window
[gc, last_echo) is blanked.  The transmitted pulse and the clutter before gc
are never touched, as the inversion takes its transmit template from them.

The noise-injection calibration must not use blanked samples.  Its background
window is the last 500 samples of the echo window, its injection window lies
beyond, and during the injection the noise is some 8 times higher, so an
impulse that stands out in the background drowns in the injection: blanking
would remove impulses from one window and not the other, and bias alpha.
Unblanked, the impulses are in both and cancel in alpha.  lpi_files therefore
takes the calibration windows from the unwrapped reader (blank_impulses).

For power over the running median, noise alone exceeds 30 in about one sample
in 2e9; sums of blanked and read samples are kept in n_blanked and n_samples.
"""

import numpy as n
from scipy.ndimage import median_filter, binary_dilation


class BlankingReader:
    def __init__(self, reader, dirname, thresh=30.0, guard=2, window=65):
        from raw_reader import pulse_index
        from radar_timing import TMM
        self.reader = reader
        ix = pulse_index(dirname)
        self.keys, self.sids = ix["key"], ix["sweepid"]
        self.tmm = TMM
        self.thresh, self.guard, self.window = thresh, guard, window
        self.n_blanked = 0
        self.n_samples = 0

    def _segments(self, key, length):
        i = n.searchsorted(self.keys, key)
        if i >= len(self.keys) or self.keys[i] != key or int(self.sids[i]) not in self.tmm:
            return []
        t = self.tmm[int(self.sids[i])]
        a, b = t["gc"], min(t["last_echo"], length)
        return [(a, b)] if a < b else []

    def read_vector_1d(self, start_sample, vector_length, channel_name, sub_channel=0):
        v = self.reader.read_vector_1d(start_sample, vector_length, channel_name)
        for a, b in self._segments(int(start_sample), len(v)):
            p = n.abs(v[a:b])**2
            if len(p) < self.window:
                continue
            hot = p > self.thresh * median_filter(p, size=self.window, mode="nearest")
            if hot.any():
                if self.guard:
                    hot = binary_dilation(hot, iterations=self.guard)
                v[a:b][hot] = 0
                self.n_blanked += int(hot.sum())
            self.n_samples += len(p)
        return v
