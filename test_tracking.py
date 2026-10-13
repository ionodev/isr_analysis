import numpy as n

from satellite_columns import TrackedDetections

LAM = 299792458.0 / 440.2e6


class Rows:
    def __init__(self, rows):
        self.r = rows

    def rows(self, i0, i1):
        return [x for x in self.r if i0 <= x[0] < i1]


def a_pass(t0_s, dur_s, delay0, doppler, peak_db, keys, threshold_db=19.0):
    """Every pulse of the pass, and the detected ones (SNR above threshold)."""
    rate = -doppler * LAM / 299792458.0 * 1e6        # samples per second
    tc = t0_s + dur_s / 2
    truth, det = {}, []
    for k in keys:
        t = k / 1e6
        if t0_s <= t < t0_s + dur_s:
            x = (t - tc) / (dur_s / 2)
            snr = peak_db - 12 * x ** 2
            d = delay0 + rate * (t - tc)
            truth[k] = d
            if snr >= threshold_db:
                det.append((k, d, doppler, snr))
    return truth, det


KEYS = list(range(0, 10_000_000, 10_000))            # a pulse every 10 ms for 10 s


def test_undetected_pulses_of_a_pass_get_columns_at_the_right_delay():
    truth, det = a_pass(3.0, 4.0, 3000.0, -500.0, 23.0, KEYS)
    td = TrackedDetections(Rows(det), lambda i0, i1: [k for k in KEYS if i0 <= k < i1], margin_s=0)
    out = td.for_period(0, 10_000_000)
    assert len(det) < 0.7 * len(truth)               # the detector misses the flanks
    covered = [k for k in truth if k in out]
    assert len(covered) >= 0.95 * len(truth)
    err = [abs(out[k][0][0] - truth[k]) for k in covered]
    assert max(err) < 0.5                             # samples


def test_detections_with_gaps_form_one_pass():
    truth, det = a_pass(3.0, 4.0, 3000.0, -500.0, 23.0, KEYS)
    det = [d for i, d in enumerate(det) if (i // 10) % 3 != 1]   # holes of 10 pulses
    td = TrackedDetections(Rows(det), lambda i0, i1: [k for k in KEYS if i0 <= k < i1], margin_s=0)
    assert len(td.passes(det)) == 1


def test_two_passes_at_different_dopplers_stay_apart():
    _, a = a_pass(3.0, 4.0, 3000.0, -500.0, 23.0, KEYS)
    _, b = a_pass(3.0, 4.0, 3003.0, 4000.0, 23.0, KEYS)
    td = TrackedDetections(Rows(a + b), lambda i0, i1: KEYS, margin_s=0)
    tr = td.passes(a + b)
    assert len(tr) == 2
    for t in tr:
        assert len(set(round(r[2]) for r in t)) == 1


def test_isolated_detections_get_no_track():
    det = [(1_000_000, 3000.0, -500.0, 20.0), (6_000_000, 2000.0, 300.0, 20.0)]
    td = TrackedDetections(Rows(det), lambda i0, i1: KEYS, margin_s=0)
    out = td.for_period(0, 10_000_000)
    assert sorted(out) == [1_000_000, 6_000_000]
