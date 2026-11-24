import numpy as n

from satellite_columns import TrackedDetections

KEYS = list(range(0, 20_000_000, 8910))


class Rows:
    def __init__(self, rows):
        self.r = rows

    def rows(self, i0, i1):
        return [x for x in self.r if i0 <= x[0] < i1]


def keys(i0, i1):
    return [k for k in KEYS if i0 <= k < i1]


def run(rows, **kw):
    tr = TrackedDetections(Rows(rows), keys, min_weak_pulses=3, **kw)
    return tr.for_period(0, 10_000_000)


def weak_pass(k0, n_det, delay=3000.0, f=0.0, step=20):
    return [(KEYS[k0 + i * step], delay, f, 12.0, False) for i in range(n_det)]


def test_three_linked_weak_detections_make_a_pass():
    out = run(weak_pass(300, 3))
    assert len(out) > 3          # the pass got its tracked columns


def test_two_weak_detections_are_dropped():
    assert run(weak_pass(300, 2)) == {}


def test_weak_detections_on_one_pulse_count_once():
    k = KEYS[300]
    rows = [(k, 3000.0, 0.0, 12.0, False), (k, 3003.0, 100.0, 12.0, False), (KEYS[320], 3000.0, 0.0, 12.0, False)]
    assert run(rows) == {}


def test_a_weak_detection_joins_a_strong_pass():
    rows = [(KEYS[300], 3000.0, 0.0, 25.0, True), (KEYS[320], 3000.0, 0.0, 12.0, False)]
    out = run(rows)
    assert KEYS[320] in out and KEYS[300] in out


def test_a_single_strong_detection_is_kept():
    out = run([(KEYS[300], 3000.0, 0.0, 25.0, True)])
    assert KEYS[300] in out


def test_unlinked_weak_detections_are_dropped():
    # three weak detections far apart in delay: three passes of one each
    rows = [(KEYS[300], 1500.0, 0.0, 12.0, False), (KEYS[320], 3000.0, 0.0, 12.0, False),
            (KEYS[340], 4500.0, 0.0, 12.0, False)]
    assert run(rows) == {}


def test_without_confirmation_rows_with_the_mark_behave_as_before():
    rows4 = [(KEYS[300], 3000.0, 0.0, 25.0)]
    rows5 = [(KEYS[300], 3000.0, 0.0, 25.0, True)]
    a = TrackedDetections(Rows(rows4), keys).for_period(0, 10_000_000)
    b = TrackedDetections(Rows(rows5), keys).for_period(0, 10_000_000)
    assert a == b


def test_cfar_alpha_is_the_detectors_threshold():
    import single_pulse_satellite_cfar as cfar
    rng = n.random.default_rng(0)
    L = 60
    power = rng.exponential(size=(1200, 120)).astype(n.float32)
    freqs = n.arange(120, dtype=float)
    for pfa in (1e-10, 1e-6, 1e-4):
        _, alpha = cfar.cfar_peaks(power, freqs, L, pfa, 32, training_padding=cfar.ceil_div(128, 8))
        assert cfar.cfar_alpha(pfa, L, 8) == float(alpha)
