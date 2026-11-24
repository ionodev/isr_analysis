import numpy as n

import plasma_line_density_scale as pls


def test_uncorrected_density_is_f_over_8_98_squared():
    assert abs(pls.plasma_density(8.98e6) - 1e12) < 1e3


def test_correction_inverts_the_dispersion_relation():
    ne, Te = 6e11, 2500.0
    fp2 = (8.98 * n.sqrt(ne)) ** 2
    f = n.sqrt(fp2 + pls.line_offset_hz2(Te))
    assert abs(pls.plasma_density(f, Te) / ne - 1) < 1e-9
    # the line reads 2-3 % high in density at 7 MHz and 2500 K
    assert 0.02 < pls.plasma_density(f) / ne - 1 < 0.03


def test_peak_of_a_log_parabola():
    r = n.arange(150.0, 600.0, 9.0)
    lr, lr0 = n.log10(r), n.log10(270.0)
    ne = 10 ** (3.0 - 40.0 * (lr - lr0) ** 2)
    pk, hm = pls.peak(r, ne, n.ones(len(r), bool))
    assert abs(pk / 1e3 - 1) < 1e-6 and abs(hm - 270.0) < 1e-3


def test_peak_needs_a_maximum():
    r = n.arange(150.0, 600.0, 9.0)
    pk, hm = pls.peak(r, r ** 2, n.ones(len(r), bool))   # rises to the 500 km limit
    assert n.isnan(pk) and n.isnan(hm)


def test_rejection_drops_outliers_and_keeps_the_scatter():
    rng = n.random.default_rng(3)
    x = rng.normal(0.0, 0.02, 1000)
    x[:5] = [0.3, -0.4, 0.5, 1.0, -0.2]
    keep = pls.reject_outliers(x)
    assert not keep[:5].any()
    assert keep[5:].mean() > 0.99


def test_nearest_respects_the_limit():
    t_ref = n.array([0.0, 60.0, 120.0])
    assert list(pls.nearest(t_ref, n.array([-100.0, 25.0, 95.0, 400.0]), 60.0)) == [-1, 0, 2, -1]
