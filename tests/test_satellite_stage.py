import run_analysis


class FakeCatalogue:
    def __init__(self, path, channel, delay, min_snr_db=None):
        self.args = (path, channel, delay, min_snr_db)


def test_no_satellite_section_means_no_stage():
    assert run_analysis.satellite_detections("/data", {"channel": "zenith-l"}) is None


def test_disabled_section_means_no_stage():
    cfg = {"channel": "zenith-l", "satellites": {"catalogue": "/cat", "enabled": False}}
    assert run_analysis.satellite_detections("/data", cfg) is None


def test_given_delay_is_used(monkeypatch):
    import satellite_columns
    monkeypatch.setattr(satellite_columns, "CatalogueDetections", FakeCatalogue)
    cfg = {"channel": "misa-l", "satellites": {"catalogue": "/cat", "channel_delay_us": 12.34, "min_snr_db": 10}}
    det = run_analysis.satellite_detections("/data", cfg)
    assert det.args == ("/cat", "misa-l", 12.34, 10)


def test_missing_delay_is_measured(monkeypatch):
    import satellite_columns, tx_delay
    monkeypatch.setattr(satellite_columns, "CatalogueDetections", FakeCatalogue)
    monkeypatch.setattr(tx_delay, "estimate_channel_delay", lambda d, ch, verbose=False: (10.62, 0.05, 100))
    cfg = {"channel": "zenith-l", "satellites": {"catalogue": "/cat"}}
    assert run_analysis.satellite_detections("/data", cfg).args == ("/cat", "zenith-l", 10.62, None)


def test_barrier_returns_with_one_process():
    run_analysis.barrier()


def test_steps_wait_for_each_other(monkeypatch):
    calls = []
    monkeypatch.setattr(run_analysis, "run_lpi", lambda *a, **k: calls.append("lpi"))
    monkeypatch.setattr(run_analysis, "run_fit_lpi", lambda *a, **k: calls.append("fit_lpi"))
    monkeypatch.setattr(run_analysis, "barrier", lambda: calls.append("barrier"))
    monkeypatch.setattr(run_analysis, "load_config", lambda p: {
        "data_dir": "/d", "steps": {"lpi": {"enabled": True, "range_gate_us": 60},
                                    "fit_lpi": {"enabled": True}}})
    monkeypatch.setattr("sys.argv", ["run_analysis.py", "cfg.json"])
    run_analysis.main()
    assert calls == ["lpi", "barrier", "fit_lpi"]
