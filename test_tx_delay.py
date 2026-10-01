import numpy as np
import pytest

import tx_delay

S = 1_000_000          # metadata samples per second
IPP = 10_000           # a pulse every 10 ms
OWN_US, OTHER_US = 10.8, 9.6


class FakeMeta:
    def __init__(self, path, length_s, reads):
        self.length, self.reads = int(length_s * S), reads

    def get_samples_per_second(self):
        return S

    def get_bounds(self):
        return (0, self.length)

    def read(self, i0, i1, field):
        self.reads.append((i1 - i0) / S)
        return {k: 1 for k in range(i0, i1, IPP)}   # sweep 1: coded


def run(monkeypatch, own_from_s, length_s=3000.0, max_search_s=2400.0):
    """Own-antenna pulses from own_from_s on (None: never)."""
    reads = []
    own = (lambda k: own_from_s is not None and k >= own_from_s * S)

    class FakeRF:
        def __init__(self, path):
            pass

        def get_channels(self):
            return ["tx-h", "zenith-l"]

        def get_properties(self, ch):
            return {"samples_per_second": 1e6}

        def read_vector_1d(self, start, length, ch):
            return np.full(length, 1.0 if own(start - 76) else 0.0)

    monkeypatch.setattr(tx_delay, "DigitalRFReader", FakeRF)
    monkeypatch.setattr(tx_delay, "DigitalMetadataReader",
                        lambda path: FakeMeta(path, length_s, reads))
    monkeypatch.setattr(tx_delay, "matched_filter_delay",
                        lambda z_tx, z_echo, oversample, max_lag:
                        ((OWN_US if z_echo[0] == 1.0 else OTHER_US), 1.0))
    ant = lambda k: -1.0 if own(k) else 1.0    # zenith-l: own antenna is -1
    pwr = lambda t: 1e6
    d, spread, n_used = tx_delay.estimate_channel_delay(
        "x", "zenith-l", n_pulses=100, verbose=False, zpm=pwr, mpm=pwr,
        tx_ant=ant, rx_ant=ant, max_search_s=max_search_s)
    return d, n_used, reads


def test_own_pulses_in_the_first_window(monkeypatch):
    d, n_used, reads = run(monkeypatch, own_from_s=0)
    assert reads == [20.0]
    assert d == pytest.approx(OWN_US) and n_used == 100


def test_window_grows_until_it_holds_own_pulses(monkeypatch):
    # on main the 20 s window held no own pulse and fell back to the other
    # antenna's delay
    d, n_used, reads = run(monkeypatch, own_from_s=100)
    assert reads == [20.0, 40.0, 80.0, 160.0]
    assert d == pytest.approx(OWN_US) and n_used == 100


def test_cap_then_fallback(monkeypatch):
    d, n_used, reads = run(monkeypatch, own_from_s=None)
    assert reads == [20.0, 40.0, 80.0, 160.0, 320.0, 640.0, 1280.0, 2400.0]
    assert d == pytest.approx(OTHER_US)


def test_recording_shorter_than_the_window(monkeypatch):
    d, n_used, reads = run(monkeypatch, own_from_s=None, length_s=30.0)
    assert reads == [20.0, 30.0]
    assert d == pytest.approx(OTHER_US)
