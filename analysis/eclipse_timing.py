#!/usr/bin/env python3
"""
The 8 April 2024 eclipse at Millstone Hill: start, maximum and end of the
partial eclipse, and the obscuration (the fraction of the Sun's disk the
Moon covers) against time, at the ground and at a height above the site.

The Sun's and the Moon's apparent positions and distances are computed with
astropy (its built-in ephemeris) for an observer at the site, at the given
height straight above it; the obscuration is the area of overlap of the two
disks (radii from the bodies' radii and distances) over the Sun's disk.
The Sun's altitude is printed too: at heights where the Sun is above the
horizon the eclipse geometry at that height is what matters for the
ionosphere's photoionisation there.

    python3 eclipse_timing.py [height_km ...]        (default 0 300)
    python3 eclipse_timing.py --json out.json [height_km ...]
"""
import json
import sys

import numpy as n
import astropy.units as u
from astropy.coordinates import EarthLocation, get_body, AltAz
from astropy.time import Time

# Millstone Hill (MIT Haystack Observatory, Westford, Massachusetts)
LAT, LON, ALT_M = 42.619, -71.491, 146.0
R_SUN_KM = 695700.0
R_MOON_KM = 1737.4


def overlap_fraction(rs, rm, d):
    """Area of overlap of disks of angular radii rs (Sun) and rm (Moon) whose
    centres are d apart, over the Sun's area (all in the same angle unit)."""
    rs, rm, d = n.broadcast_arrays(n.asarray(rs, float), n.asarray(rm, float), n.asarray(d, float))
    out = n.zeros(rs.shape)
    full = d <= n.abs(rm - rs)
    out[full] = n.minimum(1.0, (n.minimum(rm, rs)[full] / rs[full]) ** 2)
    part = (d < rs + rm) & ~full
    r1, r2, dd = rs[part], rm[part], d[part]
    a1 = r1 ** 2 * n.arccos(n.clip((dd ** 2 + r1 ** 2 - r2 ** 2) / (2 * dd * r1), -1, 1))
    a2 = r2 ** 2 * n.arccos(n.clip((dd ** 2 + r2 ** 2 - r1 ** 2) / (2 * dd * r2), -1, 1))
    a3 = 0.5 * n.sqrt(n.maximum((-dd + r1 + r2) * (dd + r1 - r2) * (dd - r1 + r2) * (dd + r1 + r2), 0))
    out[part] = (a1 + a2 - a3) / (n.pi * r1 ** 2)
    return out


def obscuration(times, height_km):
    loc = EarthLocation.from_geodetic(LON * u.deg, LAT * u.deg, ALT_M * u.m + height_km * u.km)
    sun = get_body("sun", times, loc)
    moon = get_body("moon", times, loc)
    d = sun.separation(moon).to(u.rad).value
    rs = n.arcsin(R_SUN_KM / sun.distance.to(u.km).value)
    rm = n.arcsin(R_MOON_KM / moon.distance.to(u.km).value)
    alt = sun.transform_to(AltAz(obstime=times, location=loc)).alt.deg
    return overlap_fraction(rs, rm, d), alt


def contacts(height_km, step_s=10.0):
    t = Time("2024-04-08T17:00:00", scale="utc") + n.arange(0, 5 * 3600, step_s) * u.s
    f, alt = obscuration(t, height_km)
    on = n.where(f > 0)[0]
    i0, i1, im = on[0], on[-1], int(n.argmax(f))
    return dict(height_km=height_km, start=t[i0].isot[:19], maximum=t[im].isot[:19], end=t[i1].isot[:19],
                max_obscuration=float(f[im]), sun_altitude_at_max_deg=float(alt[im]),
                times=[x[:19] for x in t.isot], obscuration=f.tolist())


def main():
    args = sys.argv[1:]
    out = None
    if args and args[0] == "--json":
        out, args = args[1], args[2:]
    heights = [float(h) for h in args] or [0.0, 300.0]
    res = [contacts(h) for h in heights]
    for r in res:
        print("%5.0f km: start %s, maximum %s (obscuration %.3f, Sun at %.1f deg), end %s UTC" % (
            r["height_km"], r["start"], r["maximum"], r["max_obscuration"], r["sun_altitude_at_max_deg"], r["end"]))
    if out:
        json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
