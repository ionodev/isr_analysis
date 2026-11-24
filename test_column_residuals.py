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


def test_residual_flag_marks_the_gates_of_a_column_left_with_power():
    import fit_quality as fq
    r = n.arange(100.0, 800.0, 9.0)
    t0 = fq.T0_RECORDING + 10 * n.array([100.0, 106.0])
    F = dict(t0=t0, t1=t0 + 60, r=r)
    for p in ("ne", "dne", "Te", "dTe", "Ti", "dTi", "vi", "dvi"):
        F[p] = n.ones((2, len(r)))
    F["Te"] = F["Te"] * 1500.0
    F["Ti"] = F["Ti"] * 1000.0
    F["dne"] = F["dne"] * 0.1
    cols = n.array([[102, 450.0, 5.0],        # left with power, in window 0
                    [103, 300.0, 0.9],        # removed
                    [108, 600.0, 3.0]])       # left, in window 1
    f = fq.flags(F, columns=cols)
    res = (f & fq.RESIDUAL) > 0
    assert res[0, n.abs(r - 450) <= 9].all() and res[0].sum() == (n.abs(r - 450) <= 9).sum()
    assert res[1, n.abs(r - 600) <= 9].all() and res[1].sum() == (n.abs(r - 600) <= 9).sum()
    # information only: a gate is still good
    assert fq.good(f).all()
