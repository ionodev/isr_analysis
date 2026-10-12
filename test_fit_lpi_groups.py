"""
Tests of how fit_lpi groups the LPI files into fit windows, and when an
incremental run (reanalyze false) fits a window again.

Run with: python -m pytest test_fit_lpi_groups.py
"""
import h5py

import fit_lpi


def test_groups():
    g = fit_lpi.group_by_time
    assert g([], 60) == []
    assert g([0.0], 60) == [[0]]
    # 10 s files in 60 s windows; the short last window is kept
    assert g([0, 10, 20, 30, 40, 50, 60, 70], 60) == [[0, 1, 2, 3, 4, 5], [6, 7]]
    # a file exactly max_dt after the window start opens the next window
    assert g([0, 60], 60) == [[0], [1]]
    # a gap: the next window starts at the next file, not on a fixed grid
    assert g([0, 10, 200, 210, 300], 60) == [[0, 1], [2, 3], [4]]


def _pp(path, n_files=None):
    with h5py.File(path, "w") as h:
        h["ne"] = [1.0]
        if n_files is not None:
            h["n_lpi_files"] = n_files


def test_refit_decision(tmp_path):
    f = str(tmp_path / "pp-0.h5")
    # no fit yet
    assert not fit_lpi.fit_is_current(f, 1)
    # a fit written before the count was stored
    _pp(f)
    assert not fit_lpi.fit_is_current(f, 1)
    # the last window of an earlier run, which has since gained files
    _pp(f, 2)
    assert not fit_lpi.fit_is_current(f, 6)
    assert fit_lpi.fit_is_current(f, 2)
    assert fit_lpi.fit_is_current(f, 1)
    # an unreadable file is fitted again
    with open(f, "w") as o:
        o.write("not hdf5")
    assert not fit_lpi.fit_is_current(f, 1)
