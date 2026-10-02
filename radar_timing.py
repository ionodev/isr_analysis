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
TMM[300] = {"noise0": 7800, "noise1": 8371, "tx0": 76, "tx1": 645,
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
