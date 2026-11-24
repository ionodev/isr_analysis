from radar_timing import TMM

# the diode's switch-on in mode 300, the 50 % point, 7813-7816 over 7-9 April
# (memo 27; gate B review, 2026-10-02)
DIODE_ON_300 = 7815
# the raw power settles by about 7822, but the range-Doppler path's 60 kHz
# low-pass filter smears the switch-on burst: 35 us after switch-on, the
# residual on misa-l is down to 0.2 % (0.85 % at 7830)
FILTER_SETTLE_US = 35


def test_mode300_injection_window_after_the_filtered_burst():
    m = TMM[300]
    assert m["noise0"] >= DIODE_ON_300 + FILTER_SETTLE_US
    assert m["noise0"] < m["noise1"] <= m["read_length"]
    assert m["noise1"] - m["noise0"] >= 500


def test_windows_inside_the_read():
    for sweep, m in TMM.items():
        assert m["last_echo"] < m["noise0"] < m["noise1"] <= m["read_length"], sweep
