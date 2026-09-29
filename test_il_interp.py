"""
Tests of the multilinear interpolation in il_interp.ilint.getspec.

The table is four dimensional, (ne, ion fraction, te/ti, ti), so every
evaluation mixes the 16 surrounding table entries. The weight of an entry must
be the product of one weight per axis, the 16 weights must sum to one, and at a
grid point the table entry itself must come back. The last test compares the
interpolated spectra with spectra computed directly by isr_spec.py.

Run with: python -m pytest test_il_interp.py
"""
import os

import numpy as np
import pytest

import il_interp
import isr_spec


def _synthetic_table(n_freq=5, seed=1):
    """
    An ilint on the production grid whose table holds a multilinear function,

        L[i,j,k,l,m] = (1 + a_m ne_i)(1 + b_m fr_j)(1 + c_m r_k)(1 + d_m ti_l).

    A product of functions that are each linear in one parameter is reproduced
    exactly by multilinear interpolation, at any point inside the grid. A
    corner with a wrong weight or a wrong index breaks that.
    """
    rng = np.random.default_rng(seed)
    il = il_interp.ilint.__new__(il_interp.ilint)
    il.ne = 10 ** np.linspace(8, 12, num=13)
    il.mol_fracs = np.linspace(0, 1, num=10)
    il.te_ti_ratio = np.linspace(1, 3, num=6)
    il.tis = np.linspace(100, 4000, num=20)
    for name, ax in (("mol_fracs", il.mol_fracs),
                     ("te_ti_ratio", il.te_ti_ratio),
                     ("tis", il.tis)):
        setattr(il, name + "0", ax.min())
        setattr(il, name + "N", len(ax))
    il.dmol_fracs = np.diff(il.mol_fracs)[0]
    il.dte_ti_ratio = np.diff(il.te_ti_ratio)[0]
    il.dtis = np.diff(il.tis)[0]
    il.radar_freq = 440.2e6
    il.radar_wavelength = 299792458.0 / il.radar_freq

    a, b, c, d = rng.uniform(0.5, 2.0, size=(4, n_freq))
    fne = 1 + a[None, :] * il.ne[:, None] / 1e12
    ffr = 1 + b[None, :] * il.mol_fracs[:, None]
    fr = 1 + c[None, :] * il.te_ti_ratio[:, None]
    fti = 1 + d[None, :] * il.tis[:, None] / 1e3
    il.S = (fne[:, None, None, None, :] * ffr[None, :, None, None, :]
            * fr[None, None, :, None, :] * fti[None, None, None, :, :])
    il.A = il.S.astype(np.complex64)

    def exact(ne, fr_, r, ti):
        return ((1 + a * ne / 1e12) * (1 + b * fr_) * (1 + c * r)
                * (1 + d * ti / 1e3))
    return il, exact


def _unnormalized(il, ne, te, ti, fr):
    """getspec output with the ion-line power scaling divided back out."""
    s = il.getspec(ne=np.array([ne]), te=np.array([te]), ti=np.array([ti]),
                   ion1_frac=np.array([fr]), normalize=False)[0]
    alpha2 = (il_interp.sc.epsilon_0 * il_interp.sc.k * te / (ne * il_interp.sc.e**2)
              * (4 * np.pi / il.radar_wavelength) ** 2)
    return s / (ne / ((1 + alpha2) * (1 + alpha2 + te / ti)))


def _indices(il, ne, te, ti, fr):
    """The fractional table indices getspec works with, for one point."""
    ne = np.array([ne])
    idx = np.searchsorted(il.ne, ne[0])
    ne_idxl = np.array([max(idx - 1, 0)])
    ne_idxh = np.array([min(idx, len(il.ne) - 1)])
    fri = np.array([(fr - il.mol_fracs0) / il.dmol_fracs])
    ri = np.array([(te / ti - il.te_ti_ratio0) / il.dte_ti_ratio])
    tii = np.array([(ti - il.tis0) / il.dtis])
    return ne, ne_idxl, ne_idxh, fri, ri, tii


# points strictly inside a cell on every axis, with n_e both on a grid node
# (the value the fitters pass) and between nodes
MID_CELL = [
    # ne, te/ti, ti, ion fraction
    (1e12, 1.55, 1000.0, 0.5),
    (1e12, 2.93, 3700.0, 0.97),
    (1e11, 1.13, 420.0, 0.03),
    (3e11, 2.05, 1650.0, 0.61),
    (2.5e9, 1.71, 777.0, 0.29),
]


@pytest.mark.parametrize("ne,r,ti,fr", MID_CELL)
def test_corner_weights_mid_cell(ne, r, ti, fr):
    il, _ = _synthetic_table()
    ne_a, ne_idxl, ne_idxh, fri, ri, tii = _indices(il, ne, r * ti, ti, fr)
    idx, c = il.corner_weights(0, ne_a, ne_idxl, ne_idxh, fri, ri, tii)

    # expected axis weights, low/high corner of each axis
    lo = ne_idxl[0]
    hi = ne_idxh[0]
    t = (ne - il.ne[lo]) / (il.ne[hi] - il.ne[lo])
    axes = [
        ((lo, hi), (1 - t, t)),
        ((int(np.floor(fri[0])), int(np.floor(fri[0])) + 1), (1 - fri[0] % 1, fri[0] % 1)),
        ((int(np.floor(ri[0])), int(np.floor(ri[0])) + 1), (1 - ri[0] % 1, ri[0] % 1)),
        ((int(np.floor(tii[0])), int(np.floor(tii[0])) + 1), (1 - tii[0] % 1, tii[0] % 1)),
    ]
    assert np.sum(c) == pytest.approx(1.0, abs=1e-12)
    for k in range(16):
        bits = [(k >> 3) & 1, (k >> 2) & 1, (k >> 1) & 1, k & 1]
        want_idx = [axes[a][0][bits[a]] for a in range(4)]
        want_w = np.prod([axes[a][1][bits[a]] for a in range(4)])
        assert list(idx[k]) == want_idx, "corner %d" % k
        assert c[k] == pytest.approx(want_w, abs=1e-12), "corner %d" % k


@pytest.mark.parametrize("ne,r,ti,fr", MID_CELL)
def test_multilinear_function_reproduced(ne, r, ti, fr):
    il, exact = _synthetic_table()
    got = _unnormalized(il, ne, r * ti, ti, fr)
    np.testing.assert_allclose(got, exact(ne, fr, r, ti), rtol=1e-6)


def test_grid_points_synthetic():
    il, _ = _synthetic_table()
    # interior nodes on every axis but ne, where the end nodes are fine too;
    # the other axes clamp their ends 1e-4 of a cell inwards
    for i, j, k, l in [(12, 1, 1, 1), (9, 8, 4, 18), (0, 4, 2, 7), (5, 5, 3, 10)]:
        got = _unnormalized(il, il.ne[i], il.te_ti_ratio[k] * il.tis[l],
                            il.tis[l], il.mol_fracs[j])
        np.testing.assert_allclose(got, il.S[i, j, k, l, :], rtol=1e-6)


TABLE_DIR = il_interp._DEFAULT_TABLE_DIR
needs_tables = pytest.mark.skipif(
    not all(os.path.exists(os.path.join(TABLE_DIR, "ion_line_interpolate_%d_%d_440.2.h5" % m))
            for m in [(32, 16), (16, 1)]),
    reason="interpolation tables not generated")


@pytest.fixture(scope="module")
def tables():
    return {m: il_interp.ilint(ion_mass1=m[0], ion_mass2=m[1], verbose=False)
            for m in [(32, 16), (16, 1)]}


@needs_tables
@pytest.mark.parametrize("masses", [(32, 16), (16, 1)])
def test_grid_points_real_table(tables, masses):
    il = tables[masses]
    for i, j, k, l in [(12, 1, 1, 1), (9, 8, 4, 18), (12, 5, 3, 4), (6, 4, 2, 10)]:
        ne = il.ne[i]
        ti = il.tis[l]
        te = il.te_ti_ratio[k] * ti
        fr = il.mol_fracs[j]
        args = dict(ne=np.array([ne]), te=np.array([te]), ti=np.array([ti]),
                    ion1_frac=np.array([fr]))
        spec = il.getspec(acf=False, **args)[0]
        want = il.S[i, j, k, l, :] / np.max(il.S[i, j, k, l, :])
        np.testing.assert_allclose(spec, want, rtol=1e-5, atol=1e-6)
        acf = il.getspec(acf=True, **args)[0]
        want = il.A[i, j, k, l, :] / il.A[i, j, k, l, 0].real
        np.testing.assert_allclose(acf, want, rtol=1e-5, atol=1e-6)


def _direct(il, masses, ne, te, ti, fr):
    """A spectrum computed the way isr_spec.il_table fills the table."""
    fr = min(max(fr, 1e-4), 1 - 1e-4)
    plpar = {"t_i": [ti, ti], "n_e": ne, "t_e": te, "m_i": list(masses),
             "freq": il.radar_freq, "B": 45000e-9, "alpha": 90,
             "ion_fractions": [fr, 1 - fr]}
    s = isr_spec.isr_spectrum(2 * np.pi * il.doppler_hz, plpar=plpar,
                              n_points=1e3, ni_points=1e3)
    s = s - s[-1]
    return s / np.sum(s)


# relative L2 difference, interpolated against direct. at these points the
# faulty weights gave 2-24 per cent; the fixed interpolation stays below 2.5
DIRECT = [
    # ne, te/ti, ti, ion fraction
    (1e12, 1.6, 1000.0, 0.5),
    (1e12, 2.4, 1250.0, 0.95),
    (1e12, 2.8, 1600.0, 0.999),
    (1e11, 1.5, 700.0, 0.7),
]


@needs_tables
@pytest.mark.parametrize("masses", [(32, 16), (16, 1)])
def test_against_isr_spec(tables, masses):
    il = tables[masses]
    # a node first: there the table and isr_spec must agree to float precision
    ne, r, ti, fr = il.ne[12], il.te_ti_ratio[3], il.tis[5], il.mol_fracs[4]
    d = _direct(il, masses, ne, r * ti, ti, fr)
    s = _unnormalized(il, ne, r * ti, ti, fr)
    assert np.linalg.norm(s - d) / np.linalg.norm(d) < 1e-4

    for ne, r, ti, fr in DIRECT:
        d = _direct(il, masses, ne, r * ti, ti, fr)
        s = _unnormalized(il, ne, r * ti, ti, fr)
        err = np.linalg.norm(s - d) / np.linalg.norm(d)
        assert err < 0.025, "ne=%g te/ti=%g ti=%g fr=%g err=%.3f" % (ne, r, ti, fr, err)
