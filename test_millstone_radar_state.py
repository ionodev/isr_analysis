import numpy as np

from millstone_radar_state import ANTENNA_SWITCH_GUARD_S as G, _step, _unknown_at_switches

S = 1_000_000  # microseconds per second


def antenna(t, v, names, guard_s):
    """The step function get_antenna_select builds from one field's events."""
    t = np.asarray(t, dtype=np.int64)
    tt, vv = _unknown_at_switches(t, v, names, guard_s)
    return _step(tt, vv, v[0])


def test_unknown_from_guard_to_opening_event():
    # cycle on MISA opens at 10 s and closes at 40 s; zenith opens at 50 s
    t = [10 * S, 40 * S, 50 * S]
    v = [1.0, 1.0, -1.0]
    names = {10 * S: b"misa_cycle", 40 * S: b"", 50 * S: b"zenith_cycle"}
    f = antenna(t, v, names, G)
    assert f(5 * S) == 1.0
    assert f((40 - G - 0.1) * S) == 1.0
    assert f((40 - G) * S) == 0.0
    assert f(40 * S) == 0.0
    assert f(49.9 * S) == 0.0
    assert f(50 * S) == -1.0
    assert f(100 * S) == -1.0


def test_first_event_closes_a_cycle():
    # the recording starts inside a MISA cycle: its first event is a closing
    # event (as in eclipse2024, +52.3 s), and zenith opens at +60.6 s. Before
    # the guard the antenna is MISA, not unknown.
    t = [int(52.3 * S), int(60.6 * S)]
    v = [1.0, -1.0]
    names = {t[0]: b"", t[1]: b"zenith_cycle"}
    f = antenna(t, v, names, G)
    assert f(0) == 1.0
    assert f((52.3 - G - 0.1) * S) == 1.0
    assert f((52.3 - G + 0.1) * S) == 0.0
    assert f(55 * S) == 0.0
    assert f(60.6 * S) == -1.0


def test_last_event_is_held():
    t = [10 * S, 40 * S, 50 * S]
    v = [1.0, 1.0, -1.0]
    names = {10 * S: b"a", 40 * S: b"", 50 * S: b"b"}
    assert antenna(t, v, names, G)(1e6 * S) == -1.0


def test_guard_does_not_reach_back_past_the_previous_event():
    # a very short cycle: the guard would start before the opening
    # event at 10 s, so it is clamped to just after it
    t = [10 * S, 11 * S, 20 * S]
    v = [1.0, 1.0, -1.0]
    names = {10 * S: b"a", 11 * S: b"", 20 * S: b"b"}
    f = antenna(t, v, names, G)
    assert f(10 * S) == 1.0
    assert f(10 * S + 1) == 0.0
    assert f(19 * S) == 0.0


def test_no_change_without_closing_event():
    # an antenna change between two named events is left as recorded
    t = [10 * S, 20 * S]
    v = [1.0, -1.0]
    names = {10 * S: b"a", 20 * S: b"b"}
    f = antenna(t, v, names, G)
    assert f(19.9 * S) == 1.0
    assert f(20 * S) == -1.0


def test_values_are_steps():
    rng = np.random.default_rng(20261001)
    t, v, names = [], [], {}
    now, ant = 0, 1.0
    for _ in range(50):
        now += int(rng.uniform(5, 60) * S)
        t.append(now)
        v.append(ant)
        names[now] = b"cycle"
        now += int(rng.uniform(30, 300) * S)
        t.append(now)
        v.append(ant)
        names[now] = b""
        if rng.random() < 0.5:
            ant = -ant
    tt, vv = _unknown_at_switches(np.array(t, dtype=np.int64), v, names, G)
    f = _step(tt, vv, v[0])
    x = np.linspace(-10 * S, now + 10 * S, 200_000)
    assert set(np.unique(f(x))) <= {-1.0, 0.0, 1.0}
    assert np.all(np.diff(tt) >= 0)
