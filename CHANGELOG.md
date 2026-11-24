# Changelog

## 0.1.0-preview (unreleased)

The first preview release of this fork of
[jvierine/isr_analysis](https://github.com/jvierine/isr_analysis).

### Added
- Satellite echoes as unknowns in the lag-profile inversion: a single-pulse
  matched-filter CFAR detector (`pipeline/single_pulse_satellite_cfar.py`),
  pass tracking with confirmation of weak detections, and one unknown per
  echo (`pipeline/satellite_columns.py`; `steps.lpi.satellites`).
- Removal of the lagged products touched by power-line impulses
  (`pipeline/impulse_rows.py`; `steps.lpi.impulse_rows`).
- Leave-one-out variance weights, a longer weighting window and a row cap
  (`weights_exclude_self`, `weights_width`, `weights_row_cap`).
- Fit options: per-band range averaging with correlated gate errors, the
  clutter-subtracted ACF below a cutoff, a wider theory table
  (`tools/make_table.py`), restarts and best-start selection.
- Quality flags and density scaling of the fits (`pipeline/fit_quality.py`);
  the density scale from the plasma lines (`calibration/`).
- Example configurations, `requirements.txt`, a user guide and a
  configuration reference; tests on GitHub Actions.

### Changed
- Layout: only `run_analysis.py` in the top folder; the pipeline in
  `pipeline/`, other scripts in `calibration/`, `analysis/`, `plotting/`,
  `tools/`, tests in `tests/`.
- The pipeline modules no longer contain machine-specific data paths; the
  data folder is a required argument, and running a module directly points
  to `run_analysis.py`.

### Fixed
- Corrections in the inversion, the fits and the range averaging, among
  them the range-Doppler outlier blanking at the last gate, the median
  averaging of the long pulse, the corner weights of the ion-line
  interpolation and the transmitter delay per antenna.
