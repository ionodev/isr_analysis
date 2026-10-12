"""
Tests of the space object cut of fit_lpi.fit_lpifiles (space_object_flagged):
off by default, and the two versions that can be switched on.

Run with: python -m pytest test_space_object_cut.py
"""
import inspect

import pytest

import fit_lpi
import fit_ionline
from fit_lpi import space_object_flagged


def test_off_by_default():
    for module in (fit_lpi, fit_ionline):
        default = inspect.signature(module.fit_lpifiles).parameters["space_object_cut"].default
        assert default is None
    # off flags nothing, however narrow the width
    assert not space_object_flagged(0.0, 0.0, None)


def test_absolute_is_juhas_test():
    # width < 300 m/s and its uncertainty < 100 m/s
    assert space_object_flagged(299.0, 99.0, "absolute")
    assert not space_object_flagged(301.0, 10.0, "absolute")
    assert not space_object_flagged(100.0, 101.0, "absolute")


def test_significance_is_three_sigma():
    # width + 3 sigma < 300 m/s
    assert space_object_flagged(150.0, 49.0, "significance")
    assert not space_object_flagged(150.0, 51.0, "significance")
    # a narrow width with a large uncertainty is flagged by neither
    assert not space_object_flagged(250.0, 90.0, "significance")


def test_unknown_cut_is_refused():
    with pytest.raises(ValueError):
        space_object_flagged(100.0, 10.0, "two sigma")
    with pytest.raises(ValueError):
        fit_lpi.fit_lpifiles(dirn="/nonexistent", space_object_cut="two sigma")
    with pytest.raises(ValueError):
        fit_ionline.fit_lpifiles(dirn="/nonexistent", space_object_cut="two sigma")
