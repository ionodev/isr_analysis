import numpy as np

from satellite_columns import refine_delay
from tx_delay import fractional_shift


def test_refine_delay_recovers_subsample_delay():
    rng = np.random.default_rng(20260930)
    tx0, tx1, n = 76, 624, 4000
    # a binary phase code of 10 us bauds, smoothed as a receiver would
    bauds = np.repeat(rng.choice([-1.0, 1.0], (tx1 - tx0) // 10 + 1), 10)[: tx1 - tx0]
    tx = np.zeros(n, dtype=np.complex128)
    tx[tx0:tx1] = bauds
    tx = np.convolve(tx, np.ones(3) / 3, "same")
    t = np.arange(n)
    for _ in range(20):
        true = 2000.0 + rng.uniform(0, 8)        # delay of the copy of tx[tx0:tx1]
        doppler = rng.uniform(-3000, 3000)
        echo = fractional_shift(tx, true - tx0) * np.exp(2j * np.pi * doppler * t / 1e6)
        echo = echo + 0.05 * (rng.normal(size=n) + 1j * rng.normal(size=n))
        grid = 8 * np.round(true / 8)            # the catalogue's grid
        got, gain = refine_delay(echo, tx, grid, doppler, tx0, tx1)
        assert abs(got - true) < 0.15
        assert gain >= 1.0
