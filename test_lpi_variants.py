import numpy as n
import pytest

from outlier_lpi import row_sigma, design_block


def test_variance_weights_are_the_local_estimate():
    loc = n.array([[1.0, 2.0], [3.0, 4.0]])
    assert row_sigma(loc, 2.5, "variance") is loc


def test_equal_weights_are_one_common_level():
    loc = n.array([[1.0, 2.0], [3.0, 4.0]])
    out = row_sigma(loc, 2.5, "equal")
    assert out.shape == loc.shape
    assert n.all(out == 2.5)
    # the local estimate is left as it was
    assert loc[1, 1] == 4.0


def test_unknown_weighting_is_refused():
    with pytest.raises(ValueError):
        row_sigma(n.ones(3), 1.0, "median")


def test_background_column_is_a_column_of_ones():
    amb = n.arange(6, dtype=n.complex64).reshape(3, 2)
    with_bg = design_block(amb, True)
    assert with_bg.shape == (3, 3)
    assert n.all(with_bg[:, :2] == amb)
    assert n.all(with_bg[:, 2] == 1)


def test_no_background_column_leaves_the_ambiguity_rows():
    amb = n.arange(6, dtype=n.complex64).reshape(3, 2)
    assert design_block(amb, False) is amb


def test_equal_weights_give_the_ordinary_least_squares_estimate():
    # with one common sigma the weighted estimate equals the unweighted one
    rng = n.random.default_rng(1)
    A = rng.normal(size=(40, 4))
    m = rng.normal(size=40)
    sig = row_sigma(rng.uniform(0.5, 3, size=40), 1.7, "equal")
    B = A / sig[:, None]
    xw = n.linalg.solve(B.T @ B, B.T @ (m / sig))
    xo = n.linalg.lstsq(A, m, rcond=None)[0]
    assert n.allclose(xw, xo)


def test_ols_covariance_is_the_sandwich():
    from scipy import sparse
    from outlier_lpi import ols_covariance, invert_normal
    rng = n.random.default_rng(2)
    A = rng.normal(size=(60, 4)) + 1j * rng.normal(size=(60, 4))
    A[:, 2] = 0                                  # one unconstrained unknown
    s = rng.uniform(0.5, 4, size=60)             # true noise level of each row
    s0 = n.median(s)
    B = sparse.csc_matrix(A / s0)
    Sigma, unc = invert_normal((n.conj(B.T) @ B).toarray())
    C = ols_covariance(Sigma, B, (s / s0) ** 2, unc)
    # direct: (A^H A)^-1 A^H diag(s^2) A (A^H A)^-1 on the constrained unknowns
    k = [0, 1, 3]
    Ak = A[:, k]
    P = n.linalg.inv(n.conj(Ak.T) @ Ak)
    D = P @ n.conj(Ak.T) @ n.diag(s ** 2) @ Ak @ P
    assert n.allclose(C[n.ix_(k, k)], D)
    assert n.isnan(C[2, 2]) and n.isfinite(C[n.ix_(k, k)]).all()
    # with every row at the common level it is the nominal covariance
    assert n.allclose(ols_covariance(Sigma, B, n.ones(60), unc)[n.ix_(k, k)], Sigma[n.ix_(k, k)])
