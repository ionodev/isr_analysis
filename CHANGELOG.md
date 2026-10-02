# Changelog

All changes to `main` of this fork (`ionodev/isr_analysis`), newest first,
and the fixes waiting on branches. The memos named are in the project's
documents (`isr_documents`, available to project contributors).

**How to add an entry.** Every merge into `main` adds its entry here, in the
same merge (REVIEW_PROCESS.md, gate A step 4). Under the date of the merge
and the right heading, write one bullet with:
- what changed and why, in a sentence or two;
- the effect on the products, if any;
- the memo, the review record, and the commit or merge.

The headings are **Fixed** (a defect corrected), **Added** (new capability or
tool), **Changed** (behaviour or defaults changed on purpose) and **Process**
(review, documentation, repository). A product-changing fix goes under
"Waiting for gate B" while it is on its branch, and moves to its merge date
when it is merged.

## Waiting for gate B (on branches, not merged)

Each has passed gate A (code review, tests, benchmark; record in
`review/records/2026-10-02-<branch>.md` on the branch). All change the
products, so all need gate B: a scientific review, then the merge. Henrik's
approvals are from 2 October 2026; the calibration fix also needs Juha's.

### Fixed

- **`lpi-unconstrained-gates`** (c2df47d): when the outlier rejection
  empties a range gate at some lag, the normal matrix is singular and the
  whole lag was lost (36 of 78 MISA periods in a test hour lose a median of
  11 of 45 lags). The inversion now solves for the constrained unknowns and
  leaves only the unconstrained ones NaN (`invert_normal`). Memo 29, item 17.
  Approved by Henrik.
- **`mode300-injection-window`** (1184dd9): mode 300's noise-injection window
  started 13 µs before the noise diode switches on and caught its overshoot,
  so MISA's mode-300 T_sys was about 26 % too low and its range-Doppler
  densities 1.29-1.37 times too low. The window now starts at 7830 µs.
  Memo 27; the site's own T_sys log of 8 April agrees with the fix. A
  calibration change: approved by Henrik, needs Juha's approval and a full
  verification.
- **`antenna-switch-gap`** (28e658f): the antenna metadata record a change
  of antenna 2.3-20.8 s (median 8.6 s) after the radar makes it, so each
  channel took in a median of 190-244 pulses per change that were sent on
  the other antenna or not at all. The antenna is now unknown from 2.5 s
  before each cycle's closing event until the next opening event, and those
  pulses are skipped. Memo 30. Must be merged together with
  `lpi-unconstrained-gates`. Approved by Henrik.
- **`catalogue-noise-window`** (387d933): the satellite catalogue's noise
  window contained the noise diode in the coded modes and mode 800, so its
  SNRs were 1.3 and 0.7 dB too low. The window now ends at the last echo.
  Memo 28. Approved by Henrik; the stored catalogue needs a rerun or a
  correction (open).
- **`tx-delay-own-antenna`** (979c1b7): the channel-delay estimator looked
  only 20 s ahead for pulses on the receiver's own antenna before falling
  back to the other antenna's, whose leakage arrives 1.2-2.3 µs away. It now
  searches up to 40 minutes first. Memo 29. Approved by Henrik.
- **`fit-lpi-range-avg`** (b093f5c): the range average gave the lowest
  gates infinite variance and averaged 2r_a gates, off centre by half a
  gate, instead of 2r_a+1. Item 16b. Approved by Henrik.
- **`fit-lpi-last-group`** (26b726e): `fit_lpifiles` never fitted the last
  group of files of a run (the last `max_dt`). Memo 35. Waiting for
  Henrik's decision.

All seven together are on `integration-all-fixes` (3fbc1ef); their combined
benchmark is the sum of the single branches, except at one antenna switch,
where `antenna-switch-gap` needs `lpi-unconstrained-gates`.

## 2026-10-01

### Process
- The review process: two gates before anything is merged into `main`
  (REVIEW_PROCESS.md), the regression benchmark (`review/regression.py`,
  `review/benchmark.json`), review records, a pull-request template, and the
  independent code-review and scientific-review agents (`.claude/agents/`).
  Pull request #1 (merge 51b9a0c); its record 8206c84.
- The scientific reviewer's paths updated to the new documents layout. Pull
  request #2 (merge 8ce21f2).

## 2026-09-30

### Fixed
- `plasma_line_clicker` chooses its gates by range (185-560 km), not by the
  hard-coded indices 50:150, and the spectra can be built without the
  interactive clicker. Item 13 (241190d).
- `outlier_lpi` skips, and counts, pulses of modes whose echo window exceeds
  the 10 000 samples it reads: a mode 800 pulse crashed the integration
  period on misa-l. Memo 20 (171ecfb).

### Added
- `satellite_columns`: optional refinement of the catalogue's delays below
  its 8-sample grid (`refine_delay`, off by default): 3-5 times more precise
  for the coded pulses; the plasma is unchanged. Memo 29 (9c3817b).
- `validate_satellite_hours`: the satellite columns over whole hours on both
  receivers. Memo 29 (2b968f5, merge b3c4879).
- `validate_impulses_rd`: power-line impulses in the mode 300 range-Doppler
  path; they add noise, not bias. Memo 26 (ba37540).
- `millstone_radar_state`: MISA's pointing per pulse, and the zenith beam's
  pointing from satellite matches (70429cd).

### Process
- README: a short overview of the project, the data and the pipeline
  (a43af93; edited by Henrik, 0e93f51).

## 2026-09-29

### Fixed
- The noise calibration's DC offset is estimated over the integration
  period instead of per window. The old estimate biased T_sys low by a
  signal-dependent amount (6 % at a 160 K baseline, 15 % at the Cygnus A
  peak at 18 kHz), and with it the injection gain α and the density scale.
  Memo 18 (1964d01; test 2bed00e). Products made before this need the
  correction in `compare_tsys_sky.py`, or regeneration.
- The range-Doppler outlier blanking now reaches the last range gate
  (ebe3be7).
- `avg_type="median"` in the range-Doppler averaging works (385095b).
- The corner weights of the ion-line interpolation (bcaaaee).
- `fit_acf_ts` evaluates the model and Jacobian at the fitted composition,
  and keeps it (b58f1f9).

### Added
- Impulse blanking of power-line sparks (`blank_impulses`), off by default:
  impulses locked to the 60 Hz mains inflate the ACF error variance by up to
  32 %. Memo 21 (88940ed, merge 9fc6540).
- `raw_reader`: whole-second reads and a cached pulse index, 30-40 times
  faster than reading pulse by pulse. Memo 19 (10a46c3).
- `validate_rfi_lpi` and `validate_rfi_rd800`: the effect of the outside
  transmitter's line interference on the LPI and on the mode 800
  range-Doppler path. Memo 20 (a284fbb, e635bba).

## 2026-09-28

### Fixed
- The lag profile inversion's noise weights are computed in double
  precision. In single precision, roundoff near a bright satellite reached
  every pulse in the echo's rows and made the ACF wrong by up to 40σ.
  Memo 17 (f846e1c, merge 5301cc2).

### Added
- Satellite echoes as extra unknowns ("columns") in the lag profile
  inversion, with the injection and real-transit tests. Memos 14-16
  (a0ae342, 6f5ff1a, merge 5301cc2).
- The system temperature compared against a 440 MHz sky model
  (`compare_tsys_sky.py`). Memo 7 (9ab0f82, merge decea5e).

### Process
- Merged from `jvierine/isr_analysis`: lazy MPI loading of the ion-line
  tables (9c74398, merge c671ba3).

## 2026-09-17

### Fixed
- The `lpi` step could not run at all (e981091).
- The antenna selection is interpolated as the step function it is; linear
  interpolation swept through zero at every switch and dropped pulses
  (0330e30).
- `acfs_var` is read as the complex variance it is, and `var_scale` is kept
  (e2cfc6b).
- Moved off `read_vector_c81d`, which digital_rf 3 removes (c75fd76).

### Changed
- The space object cut is a significance, not an absolute number (61c79aa),
  and its blanking width is derived instead of hard-coded (66cee68).
- One transmit timing table instead of four copies (9a01db0); the code says
  what it does when a pulse code is not in the table (e943c6e).
- The design matrix of the lag profile inversion is scaled once (7657505).
- The jammer notch is configurable, and the driver block guarded (5bed66c).

### Added
- Output metadata that describes the run (9a901e8), and where an
  integration period ends (ef0420e).
- A note of what the fourth-moment term is worth (3c01c02).

## 2026-09-16

### Fixed
- The plasma parameter uncertainties repaired and completed (95aba03), and
  carried into the revived `fit_lp.py` (0ebc5db).

### Changed
- The delay between `tx-h` and each echo channel is measured instead of
  hard-coded (ffe0e35).

### Added
- Merged from `jvierine/isr_analysis`: the MPI single-pulse satellite CFAR
  detector and its metadata output (merge 6cd1897).

### Process
- Experiment configs, cached tables, plots and working notes are ignored by
  git (9ee62f5).

## Before 2026-09-16

The pipeline as written by Juha Vierinen in `jvierine/isr_analysis` (from
2023), including the driver script with a JSON configuration, the MPI table
generation and the eclipse mode. See `git log 254cca8`.
