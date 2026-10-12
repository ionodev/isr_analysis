import numpy as np

from outlier_lpi import invert_normal


def design(rng, zero_cols=()):
    B = rng.normal(size=(60, 12)) + 1j * rng.normal(size=(60, 12))
    B = B.astype(np.complex64)
    B[:, list(zero_cols)] = 0
    return B


def test_all_constrained_is_the_plain_inverse():
    B = design(np.random.default_rng(1))
    ATA = B.conj().T @ B
    Sigma, unconstrained = invert_normal(ATA)
    assert len(unconstrained) == 0
    assert np.array_equal(Sigma, np.linalg.inv(ATA))


def test_zero_columns_are_left_nan_and_the_rest_is_least_squares():
    rng = np.random.default_rng(2)
    B = design(rng, zero_cols=(2, 5))
    m = (rng.normal(size=60) + 1j * rng.normal(size=60)).astype(np.complex64)
    ATA = B.conj().T @ B
    Sigma, unconstrained = invert_normal(ATA)
    assert unconstrained.tolist() == [2, 5]
    xhat = Sigma @ (B.conj().T @ m)
    keep = [i for i in range(12) if i not in (2, 5)]
    ref = np.linalg.lstsq(B[:, keep], m, rcond=None)[0]
    assert np.allclose(xhat[keep], ref, atol=1e-4)
    assert np.all(np.isnan(xhat[[2, 5]]))
    d = np.diag(Sigma)
    assert np.all(np.isnan(d[[2, 5]]))
    assert np.all(np.isfinite(d[keep]))
    # NaN only on the diagonal of the unconstrained unknowns
    assert np.sum(np.isnan(Sigma)) == 2


def test_nan_is_not_taken_for_unconstrained():
    B = design(np.random.default_rng(3))
    ATA = B.conj().T @ B
    ATA[4, 4] = np.nan
    Sigma, unconstrained = invert_normal(ATA)
    assert len(unconstrained) == 0
