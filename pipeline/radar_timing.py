#!/usr/bin/env python3
"""
Transmit and receive gate timing for the Millstone Hill experiment modes.

One table, imported by outlier_lpi.py, avg_range_doppler_spec.py and
tx_delay.py. It used to be copied into each of them, with the copies drifting:
two carried the coded modes, two carried the horizon scanning mode, one carried
the read lengths. They agreed on every value they shared, which is luck rather
than design.

All indices are samples from the start of the interpulse period at the
1 MHz recording rate, so they are also microseconds.

    tx0, tx1     transmit gate. Wider than the pulse itself: the pulse measured
                 at half maximum on tx-h runs 103 to 581 for mode 300 and 103 to
                 2100 for mode 800, so 479 and 1998 microseconds inside gates of
                 569 and 2102. The difference is guard time.
    gc           end of the ground clutter region
    e_gc         start of the usable echo, earlier than gc for the lag profile
                 inversion which can work into the clutter
    last_echo    last sample of usable echo before the noise injection
    noise0/1     noise injection gate, used to calibrate the system temperature
    read_length  samples to read for one interpulse period of this mode

Sweepids 1 to 32 are the 16-baud alternating codes, 300 the uncoded long pulse,
both on the 8910 microsecond interpulse period; 800 is the low elevation
horizon scan on a 34600 microsecond period.

These values belong to one experiment configuration. They are constants here
because nothing in the recording states them, but a recording that did carry
them should be preferred over this table.
"""

# noise injection temperature, May 24th 2022 value
T_INJECTION = 1172.0

TMM = {}

# uncoded long pulse
# The diode switches on at 7813-7816 us in mode 300 (50 % point, over 7-9
# April) with a 2 us incoherent noise burst at 7816-7817: single-sample peaks
# 108-187 times the settled level on misa-l and 10-13 times on zenith-l,
# down to 1.12 by 7822 (gate B review, 2026-10-02). The receiver's 60 kHz
# low-pass filter in the range-Doppler path smears the burst: on misa-l it
# leaves +0.85 +/- 0.11 % of the injection power at 7830 and +0.19 +/- 0.11 %
# at 7850, so the window starts at 7850 (the LPI path is clean at both). The
# old window from 7800 took in the burst and 13-16 background samples, and
# biased misa-l's T_sys low by a factor of 1.24-1.38 over 8 April, varying
# through the day (1.31-1.38 before 15 UT, 1.24-1.29 after). The fix changes
# zenith-l's T_sys by -0.8 %.
# The satellite catalogue has its own mode table with noise0 7800, kept on
# purpose: its quiet window ends before the diode. Do not synchronise.
TMM[300] = {"noise0": 7850, "noise1": 8371, "tx0": 76, "tx1": 645,
            "gc": 1000, "last_echo": 7700, "e_gc": 800, "read_length": 10000}

# 16-baud alternating codes
for _sweepid in range(1, 33):
    TMM[_sweepid] = {"noise0": 8400, "noise1": 8850, "tx0": 76, "tx1": 624,
                     "gc": 1000, "last_echo": 8200, "e_gc": 800, "read_length": 10000}

# low elevation horizon scan
TMM[800] = {"noise0": 30176, "noise1": 32033, "tx0": 69, "tx1": 2171,
            "gc": 3721, "last_echo": 30000, "e_gc": 3721, "read_length": 40000}

# sweepids of the phase coded pulses, which have a sharp correlation peak
CODED_SWEEPIDS = list(range(1, 33))

# transmit pulse length in microseconds, measured from tx-h at half maximum
TX_PULSE_LENGTH_US = {300: 479, 800: 1998}
