import numpy as n

from satellite_columns import RandomControlColumns


class Base:
    def __init__(self, cols):
        self.cols = cols

    def for_period(self, i0, i1):
        return {k: list(v) for k, v in self.cols.items() if i0 <= k < i1}


def keys(i0, i1):
    return list(range(i0, i1, 8910))


def test_random_columns_count_range_and_clearance():
    base = Base({1000000 + 8910 * 5: [(2000.0, 1000.0)], 1000000 + 8910 * 6: [(3000.0, -5000.0)]})
    rc = RandomControlColumns(base, keys, per_period=200, seed=3)
    out = rc.for_period(1000000, 11000000)
    rand = rc.random_columns(1000000, 11000000, base.for_period(1000000, 11000000))
    assert sum(len(v) for v in rand.values()) == 200
    # the base's columns are all kept
    for k, v in base.cols.items():
        for c in v:
            assert c in out[k]
    ks = n.array(keys(1000000, 11000000))
    for k, v in rand.items():
        for d, f in v:
            assert 1000.0 <= d <= 5337.0 and abs(f) <= 60e3
            pi = int(n.flatnonzero(ks == k)[0])
            for k2 in ks[max(pi - 3, 0):pi + 4]:
                for d2, f2 in base.cols.get(int(k2), []):
                    assert abs(d - d2) >= 480.0 and abs(f - f2) >= 2000.0


def test_random_columns_are_reproducible_and_differ_by_period():
    rc = RandomControlColumns(None, keys, per_period=50, seed=1)
    a = rc.for_period(1000000, 11000000)
    b = rc.for_period(1000000, 11000000)
    c = rc.for_period(11000000, 21000000)
    assert a == b
    assert sorted(v for vs in a.values() for v in vs) != sorted(v for vs in c.values() for v in vs)
