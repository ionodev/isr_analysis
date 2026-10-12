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
