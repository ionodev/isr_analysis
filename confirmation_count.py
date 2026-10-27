#!/usr/bin/env python3
"""
False alarms after the tracking confirmation: on an interval, the columns
the tracker gives with a catalogue made at a lower false-alarm probability
(weak detections confirmed by min_weak_pulses), against those it gives with
the default catalogue (P_FA 1e-10). Per period: detections, strong ones,
detections kept after confirmation, passes, and columns.

    python3 confirmation_count.py <low_pfa_catalogue> <default_catalogue> <first_period> <last_period>
"""
import sys

import numpy as n

import satellite_columns as satcol

DATA = "/mnt/data/juha/millstone_hill/isr/eclipse2024/usrp-rx0-r_20240407T100000_20240409T110000"
T0 = 1712484000


def count(det, keys, i0, i1, **kw):
    tr = satcol.TrackedDetections(det, keys, **kw)
    out = tr.for_period(i0, i1)
    return sum(len(v) for v in out.values()), len(tr.last_passes)


def main():
    low, default, a, b = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    keys = satcol.zenith_pulse_keys(DATA)
    dl = satcol.CatalogueDetections(low, "zenith-l", 10.62, strong_pfa=1e-10)
    dd = satcol.CatalogueDetections(default, "zenith-l", 10.62)
    tot = n.zeros(6)
    for ai in range(a, b + 1):
        i0 = int((T0 + 10 * ai) * 1e6)
        i1 = i0 + 10_000_000
        rows = dl.rows(i0, i1)
        strong = sum(1 for r in rows if r[4])
        c_def, p_def = count(dd, keys, i0, i1)
        c_low, p_low = count(dl, keys, i0, i1, min_weak_pulses=3)
        tot += [len(rows), strong, c_def, p_def, c_low, p_low]
        print("period %d: detections %d (strong %d); default: %d columns in %d passes; "
              "confirmed: %d columns in %d passes" % (ai, len(rows), strong, c_def, p_def, c_low, p_low), flush=True)
    m = tot / (b - a + 1)
    print("per period: detections %.1f, strong %.1f; default %.1f columns (%.1f passes); confirmed %.1f columns "
          "(%.1f passes); extra columns %.1f" % (m[0], m[1], m[2], m[3], m[4], m[5], m[4] - m[2]))


if __name__ == "__main__":
    main()
