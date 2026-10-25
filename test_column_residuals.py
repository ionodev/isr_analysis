import numpy as n
from scipy import sparse

from outlier_lpi import column_residuals


def setup(rng, shape_echo):
    n_rows, n_pl = 400, 5
    Ap = rng.normal(size=(n_rows, n_pl)) + 0j
    col = n.zeros(n_rows, dtype=complex)
    col[100:110] = n.exp(1j * n.linspace(0, 1, 10))            # the column's shape
    echo = n.zeros(n_rows, dtype=complex)
    echo[100:110] = shape_echo
    AA = sparse.csc_matrix(n.column_stack([Ap, col]))
    m = Ap @ rng.normal(size=n_pl) + 50.0 * echo + (rng.normal(size=n_rows) + 1j * rng.normal(size=n_rows)) / n.sqrt(2)
    Sinv = sparse.identity(n_rows, format="csc")
    B = Sinv @ AA
    x = n.linalg.lstsq(B.toarray(), m, rcond=None)[0]
    return column_residuals(AA, Sinv, x, m, n_pl)[0]


def test_a_described_echo_leaves_noise():
    rng = n.random.default_rng(1)
    q = [setup(rng, n.exp(1j * n.linspace(0, 1, 10))) for _ in range(200)]
    assert 0.8 < n.mean(q) < 1.0          # (10 - 1) / 10 of the noise


def test_an_echo_of_another_shape_is_left_in_the_rows():
    rng = n.random.default_rng(2)
    q = setup(rng, n.exp(1j * n.linspace(0, 6, 10)))
    assert q > 20
