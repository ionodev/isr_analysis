#!/usr/bin/env python3
"""
Sky brightness temperature seen by the Millstone Hill antennas, from PyGDSM.

The system temperature measured by the noise injection calibration is

    T_sys = T_rx + T_spill + T_atm + eta * T_sky

where only T_sky moves on a sidereal timescale. This module supplies T_sky: the
antenna beam weighted average of a global sky model, evaluated at the pointing
and time of each integration period, so that the measured T_sys can be regressed
against it.

Model
-----
GSM2008 (de Oliveira-Costa et al. 2008) with the default `haslam` basemap, which
is locked to the Haslam 408 MHz survey and carries its 1 degree resolution. At
440 MHz that is a 8 per cent extrapolation in frequency from the anchor map,
which is as close as any all sky model gets to this radar's band. GSM2016 is
available through the same interface as a cross check.

The CMB is included. PyGDSM leaves it out by default, because its models
describe galactic and extragalactic emission, but an antenna does not know the
difference: 2.725 K enters the measured T_sys like any other sky brightness, and
omitting it would push a constant into the fitted receiver temperature.

Beam weighting
--------------
At 440.2 MHz the wavelength is 0.681 m, so the diffraction limited beams are

    zenith, 68 m dish    1.22 lambda / D = 0.70 degrees
    MISA,   46 m dish    1.22 lambda / D = 1.03 degrees

Both are at or below the 1 degree resolution of the model map, which is already
the sky convolved with the Haslam survey beam. Convolving that product again
with the full antenna beam would smooth twice. What is wanted is a map whose
total resolution equals the antenna beam, so the map is smoothed by the
quadrature difference

    fwhm_extra = sqrt(fwhm_beam^2 - fwhm_map^2)

which is 0.25 degrees for MISA and imaginary for the zenith antenna. In the
latter case no smoothing is applied and the model is left slightly coarser than
the real beam: it cannot be sharpened without deconvolving Haslam. The
consequence is that the model dilutes compact bright structure, most of all the
galactic plane in Cygnus, so a fit against it should be expected to want a beam
efficiency slightly above unity rather than below.

Only the main beam is modelled. Far sidelobe and ground pickup are constant on a
sidereal timescale for a fixed pointing and are therefore absorbed by the
constant term of the fit, not by this model.
"""

import numpy as n
import healpy as hp

import millstone_radar_state as mrs

# Millstone Hill, from millstone_radar_state.py
RADAR_LAT = mrs.radar_lat
RADAR_LON = mrs.radar_lon
RADAR_HGT = mrs.radar_hgt

# radar centre frequency; the configs carry 440.2 MHz as radar_freq_hz
DEFAULT_FREQ_MHZ = 440.2

# aperture diameters, metres
DISH_DIAMETER_M = {"zenith": 68.0, "misa": 46.0}

# angular resolution of the GSM2008 haslam basemap, degrees
GSM2008_FWHM_DEG = 1.0

C_M_S = 299792458.0


def beam_fwhm_deg(antenna, freq_mhz=DEFAULT_FREQ_MHZ):
    """Diffraction limited half power beamwidth of one of the two antennas."""
    lam = C_M_S / (freq_mhz * 1e6)
    return n.degrees(1.22 * lam / DISH_DIAMETER_M[antenna])


class SkyModel:
    """
    A global sky model at one frequency, smoothed to one antenna's resolution.

    The map is generated once and reused, because generating it is far more
    expensive than sampling it and the frequency does not change within a run.
    """

    def __init__(self, freq_mhz=DEFAULT_FREQ_MHZ, model="gsm2008",
                 antenna=None, include_cmb=True, verbose=True):
        self.freq_mhz = freq_mhz
        self.model = model
        self.antenna = antenna
        self.include_cmb = include_cmb

        if model == "gsm2008":
            from pygdsm import GlobalSkyModel
            gsm = GlobalSkyModel(freq_unit="MHz", basemap="haslam",
                                 include_cmb=include_cmb)
            self.map_fwhm_deg = GSM2008_FWHM_DEG
        elif model == "gsm2016":
            from pygdsm import GlobalSkyModel2016
            gsm = GlobalSkyModel2016(freq_unit="MHz", resolution="hi",
                                     include_cmb=include_cmb)
            # 48 arcmin below 10 GHz, as documented in pygdsm.gsm16
            self.map_fwhm_deg = 48.0 / 60.0
        else:
            raise ValueError("unknown sky model %s" % (model))

        m = n.array(gsm.generate(freq_mhz)).ravel()
        self.nside = hp.npix2nside(len(m))

        self.beam_fwhm_deg = None
        self.smoothed_fwhm_deg = 0.0
        if antenna is not None:
            self.beam_fwhm_deg = beam_fwhm_deg(antenna, freq_mhz)
            extra2 = self.beam_fwhm_deg ** 2.0 - self.map_fwhm_deg ** 2.0
            if extra2 > 0:
                self.smoothed_fwhm_deg = n.sqrt(extra2)
                m = hp.smoothing(m, fwhm=n.radians(self.smoothed_fwhm_deg))

        self.map = m
        if verbose:
            print(self.describe())

    def describe(self):
        s = "%s at %.1f MHz, nside %d, map resolution %.2f deg, CMB %s" % (
            self.model, self.freq_mhz, self.nside, self.map_fwhm_deg,
            "included" if self.include_cmb else "excluded")
        if self.antenna is not None:
            s += "\n  %s beam %.2f deg, extra smoothing %.2f deg, effective %.2f deg" % (
                self.antenna, self.beam_fwhm_deg, self.smoothed_fwhm_deg,
                max(self.beam_fwhm_deg, self.map_fwhm_deg))
        return s

    def at_galactic(self, gl_deg, gb_deg):
        """
        Sky temperature at galactic coordinates, bilinearly interpolated.

        Nearest pixel sampling would quantise the model onto the 6.9 arcmin
        nside 512 grid and put a sawtooth into the model time series as the beam
        drifts across pixel boundaries.
        """
        return hp.get_interp_val(self.map, gl_deg, gb_deg, lonlat=True)

    def at_horizontal(self, t_unix, az_deg, el_deg):
        """Sky temperature for pointings given as time, azimuth and elevation."""
        gl, gb = horizontal_to_galactic(t_unix, az_deg, el_deg)
        return self.at_galactic(gl, gb)


def horizontal_to_galactic(t_unix, az_deg, el_deg):
    """
    Galactic longitude and latitude of a horizontal pointing at Millstone Hill.

    Azimuth is degrees east of north, as the antenna control metadata records it,
    and may be negative. Elevation is degrees above the horizon.
    """
    from astropy.coordinates import EarthLocation, AltAz, Galactic
    from astropy.time import Time
    import astropy.units as u

    loc = EarthLocation(lat=RADAR_LAT * u.deg, lon=RADAR_LON * u.deg,
                        height=RADAR_HGT * u.m)
    frame = AltAz(az=n.asarray(az_deg) * u.deg,
                  alt=n.asarray(el_deg) * u.deg,
                  obstime=Time(n.asarray(t_unix), format="unix"),
                  location=loc)
    g = frame.transform_to(Galactic())
    return n.asarray(g.l.deg), n.asarray(g.b.deg)


def local_sidereal_time_h(t_unix):
    """Apparent local sidereal time at the radar, hours. Used only for plotting."""
    from astropy.coordinates import EarthLocation
    from astropy.time import Time
    import astropy.units as u

    loc = EarthLocation(lat=RADAR_LAT * u.deg, lon=RADAR_LON * u.deg,
                        height=RADAR_HGT * u.m)
    t = Time(n.asarray(t_unix), format="unix", location=loc)
    return n.asarray(t.sidereal_time("apparent").hour)


if __name__ == "__main__":
    # what the two antennas can see, and how much of a signal there is to find
    for antenna in ["zenith", "misa"]:
        sm = SkyModel(antenna=antenna)
        t = n.linspace(1712484000, 1712484000 + 48 * 3600, num=2000)
        el = 90.0 if antenna == "zenith" else 45.0
        az = 0.0 if antenna == "zenith" else -90.0
        ts = sm.at_horizontal(t, n.full(len(t), az), n.full(len(t), el))
        print("  T_sky over 48 h: min %.1f  median %.1f  max %.1f K (swing %.1f K)\n"
              % (ts.min(), n.median(ts), ts.max(), ts.max() - ts.min()))
