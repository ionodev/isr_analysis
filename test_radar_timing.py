from radar_timing import TMM

# the diode's switch-on and the end of its overshoot in mode 300 (memos 27, 35)
DIODE_ON_300 = 7813
DIODE_SETTLED_300 = 7822


def test_mode300_injection_window_after_the_diode_has_settled():
    m = TMM[300]
    assert m["noise0"] > DIODE_SETTLED_300
    assert m["noise0"] < m["noise1"] <= m["read_length"]
    assert m["noise1"] - m["noise0"] >= 500


def test_windows_inside_the_read():
    for sweep, m in TMM.items():
        assert m["last_echo"] < m["noise0"] < m["noise1"] <= m["read_length"], sweep
