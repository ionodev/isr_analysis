"""
Export one day of fits into the RUISDAP Millstone viewer's quicklook format.

The viewer (github.com/jvierine/ruisdap) shows an analysis as one folder per
day and antenna:

    parameters.f32   Ne, Ti, Te, velocity, [parameter][range][time]
    raw-ne.f32       raw Ne from the run's own 10 microsecond lag, [range][time]
    manifest.json    the viewer's manifest; fits carry their 1-sigma errors
    quicklook.h5     the same arrays with units and source hashes

Density. The fit's "ne" is the proxy p = P0 (1 + Te/Ti) r_km^2 / P_tx (P0 the
zero-lag power with the receiver gain removed, in K); Ne = C p, with C the
run's density constant (fit_quality.density_constant by default).

Raw Ne is 2 C Re ACF(10 us) / alpha * r_km^2 / P_tx per 10-second inversion
(alpha the receiver gain, Te/Ti = 1), with the ground-clutter-cancelled ACF
below raw_clutter_cutoff_km.

The plasma-line panel. The viewer shows every analysis of a day next to the
day's plasma lines, which come from RUISDAP's own export of that day: its
manifest's plasma fields are copied and its spectrum file is referenced by
URL, not copied. So a day can be exported only when RUISDAP has exported its
plasma lines.

run_analysis.py runs this as the step "ruisdap_export"
(docs/CONFIGURATION.md); tools/export_ruisdap.py exports a finished run.
"""
import datetime as dt
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

import fit_quality as fq


def _finite(x):
    if isinstance(x, np.ndarray):
        return _finite(x.tolist())
    if isinstance(x, (list, tuple)):
        return [_finite(v) for v in x]
    if isinstance(x, dict):
        return {k: _finite(v) for k, v in x.items()}
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    if isinstance(x, np.generic):
        return x.item()
    return x


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _echo_rows(catalogue, channel, p0, p1):
    """Rows (period index, range km, SNR dB) of the catalogued echoes of the
    channel in the integration periods p0..p1 (fit_quality's numbering)."""
    from digital_rf import DigitalMetadataReader
    rd, rows = DigitalMetadataReader(str(catalogue)), []
    i0 = int((fq.T0_RECORDING + fq.PERIOD_S * p0) * 1e6)
    for k, rec in rd.read(i0, int((fq.T0_RECORDING + fq.PERIOD_S * p1) * 1e6),
                          columns=['channel', 'range_km', 'snr_db']).items():
        ch = rec['channel']
        if (ch.decode() if isinstance(ch, bytes) else str(ch)) != channel:
            continue
        ai = (k / 1e6 - fq.T0_RECORDING) // fq.PERIOD_S
        rows += [(ai, float(r), float(q)) for r, q in zip(rec['range_km'], rec['snr_db'])]
    return np.array(rows) if rows else None


def export_day(run, output, date, plasma_manifest, plasma_url, label, postfix='lpi_60', channel='zenith-l',
               antenna='ZENITH', magic_constant=None, calibration_method=None, source_commit='',
               quality=True, echo_catalogue=None, echo_rows=None, raw_clutter_cutoff_km=300.0):
    """
    Write the fits of one UTC day of a run into a new export folder.

    run: the run's output_dir (holds <postfix>/<channel>/pp-*.h5 and lpi-*.h5);
    output: the new export folder (must not exist); date: "YYYY-MM-DD";
    plasma_manifest, plasma_url: the day's RUISDAP export's manifest.json and
    the URL of its plasma spectrum file; label: what the run is, for the
    viewer; magic_constant: C (None: fit_quality.density_constant of the day);
    quality: write fit_quality's flags per fit and gate, with SATELLITE from
    echo_catalogue (a catalogue path) or echo_rows (a .npy of its rows) if
    given. Returns the output folder.
    """
    run, output, plasma_manifest = Path(run), Path(output), Path(plasma_manifest)
    if output.exists():
        raise FileExistsError('choose a new export folder: %s exists' % output)
    start = dt.datetime.fromisoformat(date).replace(tzinfo=dt.timezone.utc).timestamp()
    end = start + 86400
    folder = run / postfix / channel
    C = fq.density_constant(start) if magic_constant is None else magic_constant
    if calibration_method is None:
        calibration_method = ('isr_analysis fit_quality.py: the plasma-line density constant of the day'
                              if magic_constant is None else 'given by the user')
    semantics = 'plasma_line_calibrated_density' if magic_constant is None else 'calibrated_density'

    fits, sources = [], []
    for path in sorted(folder.glob('pp-*.h5')):
        t = float(path.stem.split('-')[1])
        if not start <= t < end:
            continue
        with h5py.File(path, 'r') as h:
            if not all(k in h for k in ['ne', 'Ti', 'Te', 'vi', 'rgs', 't0', 't1', 'P_tx']):
                continue
            proxy = h['ne'][()]
            params = np.stack([proxy * C, h['Ti'][()], h['Te'][()], h['vi'][()]], axis=1)
            err = np.stack([h[k][()] for k in ['dne', 'dTi', 'dTe', 'dvi']], axis=1)
            uerr = err.copy()
            err[:, 0] *= C
            digest = _sha(path)
            fits.append(dict(
                id=hashlib.sha256((digest + ':/').encode()).hexdigest()[:24],
                start_unix_s=float(h['t0'][()]), end_unix_s=float(h['t1'][()]),
                range_m=h['rgs'][()] * 1000.0, parameters=params, standard_error=err,
                uncalibrated_density_proxy=proxy, uncalibrated_standard_error=uerr,
                transmitter_w=float(h['P_tx'][()]), system_temperature_k=float(h['T_sys'][()]),
                elevation_deg=float(h['el'][()]), azimuth_deg=float(h['az'][()]),
                space_object_count=h['space_object_count'][()] if 'space_object_count' in h else [],
                n_lpi_files=int(h['n_lpi_files'][()]) if 'n_lpi_files' in h else None))
            sources.append(dict(path=str(path), sha256=digest))
    if not fits:
        raise ValueError('no fits of %s in %s' % (date, folder))
    fits.sort(key=lambda f: f['start_unix_s'])
    ranges = fits[0]['range_m']
    if any(not np.array_equal(f['range_m'], ranges) for f in fits):
        raise ValueError('fit range grids differ')

    raw, rawranges = [], None
    for path in sorted(folder.glob('lpi-*.h5')):
        t = float(path.stem.split('-')[1])
        if not start <= t < end:
            continue
        with h5py.File(path, 'r') as h:
            r = h['rgs_km'][()]
            lags = h['lags'][()]
            lag = int(np.argmin(abs(lags - 1e-5)))
            if abs(lags[lag] - 1e-5) > 1e-9:
                continue
            acf = h['acfs_e'][:, lag].copy()
            cut = np.where(r > raw_clutter_cutoff_km)[0][0]
            acf[:cut] = h['acfs_g'][:cut, lag] / 2.0
            power, alpha = float(h['P_tx'][()]), float(h['alpha'][()])
            values = 2 * C * acf.real / alpha * r * r / power if power > 0 and alpha > 0 else np.full(len(r), np.nan)
            time = float(h['t0'][()]) if 't0' in h else float(h['i0'][()])
        if rawranges is None:
            rawranges = r * 1000.0
        elif not np.array_equal(rawranges, r * 1000.0):
            raise ValueError('raw range grids differ')
        raw.append((time, values))
    raw.sort(key=lambda x: x[0])
    rawtimes = np.array([x[0] for x in raw])
    rawpixels = np.array([x[1] for x in raw]).T

    base = json.loads(plasma_manifest.read_text())
    if base['date'] != date or base['antenna'] != antenna:
        raise ValueError('the plasma manifest is for another day or antenna')

    quality_flags = None
    if quality:
        # FAILED (a parameter or its error not finite), BOUND (the fit at a
        # limit of its search), NOISY (dne/ne or dTe/Te above 0.5), and, for
        # information, SATELLITE (a catalogued echo of at least 20 dB within
        # one gate during the fit window). The viewer hides the first three
        # on request.
        F = fq.load_fits(str(run), (start, end), constant=1.0, subdir=str(Path(postfix) / channel))
        p0 = int((min(f['start_unix_s'] for f in fits) - fq.T0_RECORDING) // fq.PERIOD_S)
        p1 = int((max(f['end_unix_s'] for f in fits) - fq.T0_RECORDING) // fq.PERIOD_S) + 1
        rows = None
        if echo_rows is not None:
            rows = np.load(echo_rows)
            rows = rows[(rows[:, 0] >= p0) & (rows[:, 0] < p1)]
        elif echo_catalogue is not None:
            rows = _echo_rows(echo_catalogue, channel, p0, p1)
        flags = fq.flags(F, rows if rows is not None and len(rows) else None)
        by_t0 = {float(t): i for i, t in enumerate(F['t0'])}
        for f in fits:
            f['quality_flags'] = flags[by_t0[f['start_unix_s']]].astype(int)
        quality_flags = dict(
            bits=dict(FAILED=int(fq.FAILED), BOUND=int(fq.BOUND), NOISY=int(fq.NOISY), SATELLITE=int(fq.SATELLITE)),
            hide=['FAILED', 'BOUND', 'NOISY'],
            definition='isr_analysis fit_quality.py: FAILED, a parameter or its error is not finite; BOUND, the fit '
                       'ends at a limit of its search (Te/Ti 0.99 or 3, Ti 150 or 3000 K, |vi| 1500 m/s); NOISY, '
                       'dne/ne or dTe/Te above 0.5; SATELLITE (information only), a catalogued echo of at least '
                       '20 dB within one gate during the fit window. A gate is good when it has none of the first '
                       'three.',
            source='fit_quality.py; echoes from %s' % (echo_rows or echo_catalogue or 'none'))

    times = np.array([(f['start_unix_s'] + f['end_unix_s']) / 2 for f in fits])
    parameters = np.array([f['parameters'] for f in fits])
    output.mkdir(parents=True)
    np.asarray(parameters.transpose(2, 1, 0), dtype='<f4').tofile(output / 'parameters.f32')
    np.asarray(rawpixels, dtype='<f4').tofile(output / 'raw-ne.f32')
    manifest = dict(
        schema='ruisdap-millstone-quicklook-v1', station='Millstone Hill', antenna=antenna, date=date,
        radar_frequency_hz=440.2e6, fit_times=times, range_m=ranges, parameters_url='parameters.f32',
        parameter_names=['Ne / m⁻³', 'Ti / K', 'Te / K', 'Velocity / m s⁻¹'], fits=fits,
        raw_times=rawtimes, raw_range_m=rawranges, raw_url='raw-ne.f32',
        raw_definition='Raw Ne = 2 C Re ACF(10 us) / receiver gain * range_km^2 / transmitter_W from each '
                       '10-second isr_analysis inversion; Te/Ti = 1; ground-clutter-cancelled ACF below %g km'
                       % raw_clutter_cutoff_km,
        raw_units='m^-3',
        plasma_times=base['plasma_times'], plasma_frequency_hz=base['plasma_frequency_hz'],
        plasma_url=plasma_url, plasma_records=base['plasma_records'],
        collapse_range_km=base['collapse_range_km'], collapse_definition=base['collapse_definition'],
        density_semantics=semantics,
        frequency_convention=base.get('frequency_convention', 'signed receiver offset; no sign reversal'),
        frequency_sign_verified=base.get('frequency_sign_verified', False),
        velocity_convention='isr_analysis fit_lpi.py convention; sign not independently verified',
        calibration=dict(gain=C, revision=1, magic_constant_units='W m^-3 K^-1 km^-2', points=[],
                         method=calibration_method),
        isr_analysis=dict(label=label, run=str(run), postfix=postfix, channel=channel,
                          commit=source_commit, fit_count=len(fits), raw_count=len(raw)))
    if quality_flags is not None:
        manifest['quality_flags'] = quality_flags
    (output / 'manifest.json').write_text(json.dumps(_finite(manifest), allow_nan=False))
    with h5py.File(output / 'quicklook.h5', 'w') as h:
        h.attrs['schema'] = manifest['schema']
        h.attrs['metadata_json'] = json.dumps(
            _finite({k: v for k, v in manifest.items() if k not in ['fits', 'plasma_records']}), allow_nan=False)
        h.attrs['sources_json'] = json.dumps(sources)
        for k, v, unit in [('fit_time_unix_s', times, 'UTC Unix s'), ('range_m', ranges, 'm'),
                           ('parameters', parameters, 'm^-3,K,K,m s^-1'),
                           ('standard_error', np.array([f['standard_error'] for f in fits]), 'm^-3,K,K,m s^-1'),
                           ('raw_time_unix_s', rawtimes, 'UTC Unix s'), ('raw_range_m', rawranges, 'm'),
                           ('raw_ne', rawpixels, 'm^-3')]:
            h.create_dataset(k, data=v, compression='gzip')
            h[k].attrs['units'] = unit
    for f in output.iterdir():
        f.chmod(0o644)
    print('Export complete: %s %s, %d fits, %d raw profiles -> %s' % (date, antenna, len(fits), len(raw), output))
    return output
