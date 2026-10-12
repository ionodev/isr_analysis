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
