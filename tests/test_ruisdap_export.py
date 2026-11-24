"""ruisdap_export.export_day on two synthetic fit windows and three
inversions: the arrays, their layout and the density scaling."""
import json

import h5py
import numpy as np
import pytest

import ruisdap_export as rex

DATE = "2024-04-08"
T0 = 1712534400.0 + 18 * 3600      # 18 UTC on DATE
RGS = np.array([200.0, 300.0, 400.0])
LAGS = np.arange(1, 4) * 1e-5


def _write_run(run):
    d = run / "lpi_60" / "zenith-l"
    d.mkdir(parents=True)
    for i in range(2):
        t = T0 + 300 * i
        with h5py.File(d / ("pp-%d.h5" % t), "w") as h:
            h["ne"], h["dne"] = np.array([1e3, 2e3, 3e3]) * (i + 1), np.full(3, 10.0)
            h["Ti"], h["dTi"] = np.full(3, 1000.0), np.full(3, 50.0)
            h["Te"], h["dTe"] = np.full(3, 2000.0), np.full(3, 100.0)
            h["vi"], h["dvi"] = np.full(3, 5.0), np.full(3, 2.0)
            h["rgs"], h["t0"], h["t1"] = RGS, t, t + 200
            h["P_tx"], h["T_sys"], h["el"], h["az"] = 1e6, 150.0, 88.0, 0.0
    for i in range(3):
        t = T0 + 10 * i
        with h5py.File(d / ("lpi-%d.h5" % t), "w") as h:
            h["rgs_km"], h["lags"], h["t0"] = RGS, LAGS, t
            h["acfs_e"] = np.ones((3, 3), complex) * 2.0
            h["acfs_g"] = np.ones((3, 3), complex) * 6.0
            h["P_tx"], h["alpha"] = 1e6, 0.5


def _plasma(tmp_path):
    m = dict(date=DATE, antenna="ZENITH", plasma_times=[T0], plasma_frequency_hz=[0.0, 1.0],
             plasma_records=[dict(start=T0, end=T0 + 60)], collapse_range_km=[185, 560],
             collapse_definition="test")
    p = tmp_path / "plasma.json"
    p.write_text(json.dumps(m))
    return p


def test_export_day(tmp_path):
    run, out = tmp_path / "run", tmp_path / "export"
    _write_run(run)
    rex.export_day(run, out, DATE, _plasma(tmp_path), "/plasma.f32", "test run",
                   magic_constant=2.0, quality=False)
    m = json.loads((out / "manifest.json").read_text())
    assert m["schema"] == "ruisdap-millstone-quicklook-v1" and len(m["fits"]) == 2
    assert m["calibration"]["gain"] == 2.0 and m["plasma_url"] == "/plasma.f32"
    # [parameter][range][time]: Ne = C * proxy
    p = np.fromfile(out / "parameters.f32", dtype="<f4").reshape(4, 3, 2)
    assert np.allclose(p[0, :, 0], [2e3, 4e3, 6e3]) and np.allclose(p[0, :, 1], [4e3, 8e3, 12e3])
    assert np.allclose(p[1], 1000.0) and np.allclose(p[2], 2000.0) and np.allclose(p[3], 5.0)
    # raw Ne = 2 C Re ACF(10 us) / alpha * r^2 / P_tx; the clutter-cancelled
    # ACF / 2 at and below 300 km
    raw = np.fromfile(out / "raw-ne.f32", dtype="<f4").reshape(3, 3)
    expect = 2 * 2.0 * np.array([3.0, 3.0, 2.0]) / 0.5 * RGS ** 2 / 1e6
    assert np.allclose(raw[:, 0], expect)


def test_export_refuses_other_day_and_existing_folder(tmp_path):
    run = tmp_path / "run"
    _write_run(run)
    plasma = _plasma(tmp_path)
    with pytest.raises(ValueError):
        rex.export_day(run, tmp_path / "a", "2024-04-07", plasma, "/p", "x", magic_constant=1.0, quality=False)
    (tmp_path / "b").mkdir()
    with pytest.raises(FileExistsError):
        rex.export_day(run, tmp_path / "b", DATE, plasma, "/p", "x", magic_constant=1.0, quality=False)
