"""
Tests of the topside fit, fit_lpi.fit_acf_ts, on noise-free ACFs made with
its own model.

The fit's fifth parameter x4 gives the O+ fraction 1-1/x4. The returned
model and the Jacobian must be evaluated at that fraction, and a fit on the
Te/Ti bound of 3, the edge of the table, must still give a covariance.

Run with: python -m pytest test_fit_lpi_ts.py
"""
import os

import numpy as np
import pytest

import il_interp
import fit_lpi

needs_tables = pytest.mark.skipif(
    not all(os.path.exists(os.path.join(il_interp._DEFAULT_TABLE_DIR, "ion_line_interpolate_%d_%d_440.2.h5" % m))
            for m in [(32, 16), (16, 1)]),
    reason="interpolation tables not generated")

# the lags of the lpi_30 products
LAGS = np.arange(1, 46) * 10e-6


def _fit(te_ti, ti, o_frac, vi=50.0, zl=1e4):
    fit_lpi._init_tables(440.2e6)
    acf = zl * fit_lpi.model_acf(te_ti * ti, ti, o_frac, vi, LAGS, hplus=True)
    var = np.full(len(LAGS), (0.01 * zl) ** 2)
    guess = np.array([np.nan, np.nan, np.nan, np.nan])
    return acf, fit_lpi.fit_acf_ts(acf, LAGS, 900.0, var, guess=guess)


# topside O+ fractions. well below 0.9 the three starting points of the fit
# can end in a wrong minimum even without noise (Te/Ti and composition trade
# off), which is a property of the minimisation, not tested here
@needs_tables
@pytest.mark.parametrize("te_ti,ti,o_frac", [(2.0, 1500.0, 0.95), (1.5, 1000.0, 0.99)])
def test_composition_used_in_model(te_ti, ti, o_frac):
    acf, (res, model, sigmas, Sigma, fitted) = _fit(te_ti, ti, o_frac)
    # a noise-free ACF is fitted to far better than its 1 % noise level, and
    # the returned model is the fitted one, at the fitted composition
    assert np.max(np.abs(model - acf)) < 1e-3 * np.abs(acf[0])
    assert fitted == pytest.approx(o_frac, abs=0.01)
    assert res[0] == pytest.approx(te_ti, abs=0.02)
    assert res[1] == pytest.approx(ti, rel=0.01)
    assert np.all(np.isfinite(Sigma))


@needs_tables
def test_covariance_on_te_ti_bound():
    # Te/Ti beyond the table: the fit ends on its bound of 3
    acf, (res, model, sigmas, Sigma, o_frac) = _fit(3.5, 1200.0, 0.95)
    assert res[0] == pytest.approx(3.0, abs=1e-3)
    assert np.all(np.isfinite(Sigma))
    assert np.all(np.isfinite(sigmas)) and np.all(sigmas > 0)
    # a forward step from the bound stays on the clamped table edge and makes
    # the Te/Ti column of the Jacobian (nearly) zero: sigma(Te/Ti) came out
    # ~1e5 here, and on real data the inversion failed outright
    assert sigmas[0] < 1.0
