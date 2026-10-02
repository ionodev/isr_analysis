# Changelog

All changes to `main` of this fork (`ionodev/isr_analysis`), newest first,
and the fixes waiting on branches. The memos named are in the project's
documents (`isr_documents`, available to project contributors).

**How to add an entry.** Every merge into `main` adds its entry here, in the
same merge (REVIEW_PROCESS.md, gate A step 4). Under the date of the merge
and the right heading, write one bullet with:
- what changed and why, in a sentence or two;
- the effect on the products, if any;
- the memo, the review record, and the last commit of the change itself,
  before the commit that adds the entry (the entry cannot name the merge
  commit; find that with `git log --merges --first-parent main`).

The headings are **Fixed** (a defect corrected), **Added** (new capability or
tool), **Changed** (behaviour or defaults changed on purpose) and **Process**
(review, documentation, repository). Documentation-only merges get a
Process entry too. A product-changing fix goes under "Waiting for gate B"
while it is on its branch, and moves to its merge date, in a commit on the
branch, before it is merged. Branches that edit this file at the same time
conflict here; resolve by keeping all the entries.

## Waiting for gate B (on branches, not merged)

Each has passed gate A (code review, tests, benchmark; record in
`review/records/2026-10-02-<branch>.md` on the branch). All change the
products, so all need gate B: a scientific review, then the merge. Henrik's
approvals are from 2 October 2026 (Henrik decides on the code changes, and
which of them go to Juha for review).

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
  In the inversion, MISA's T_sys moves by 5.7-8.1 %; zenith-l by about
  1 % or less. Memo 27. The site's own T_sys log of 8 April preliminarily
  agrees (site/ours for MISA's range-Doppler T_sys: 1.27). Open: start the
  window at 7830 or at about 7850 µs. A calibration change: approved by
  Henrik; as a calibration change it gets a full verification in gate B.
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
  searches up to 40 minutes first. On the benchmark, zenith-l's delay moves
  from 11.88 to 10.62 µs, shifting every zenith-l range-Doppler spectrum by
  1.3 samples (195 m); single gates changing by up to 12-23 σ, and the
  0.1 µs against an independent 10.72 µs, are not yet explained (for gate B).
  Memo 29. Approved by Henrik, on a condition the fix meets.
- **`fit-lpi-range-avg`** (b093f5c): the range average gave the lowest
  gates infinite variance and averaged 2r_a gates, off centre by half a
  gate, instead of 2r_a+1. In a test before the benchmark, fit parameters
  moved by a median of 0.2-0.7 σ below 700 km. Item 16b. Approved by Henrik.
- **`fit-lpi-last-group`** (26b726e): `fit_lpifiles` never fitted the last
  group of files of a run (the last `max_dt`). Memo 35. Waiting for
  Henrik's decision.

All seven together are on `integration-all-fixes` (3fbc1ef); their combined
benchmark is the sum of the single branches, except at one antenna switch,
where `antenna-switch-gap` needs `lpi-unconstrained-gates`. **Merge order:**
`lpi-unconstrained-gates`, `antenna-switch-gap`, `fit-lpi-range-avg` and
`fit-lpi-last-group` go in together, since each of the other three alone
leaves fits with no lags (on its own the new misa-l last-group fit is all
NaN); `tx-delay-own-antenna` belongs with `antenna-switch-gap`.

## 2026-10-02

### Process
- This changelog, and the rule that every merge adds its entry
  (REVIEW_PROCESS.md, gate A step 4; the record and pull-request templates).
  Branch `changelog`; the process change needs Henrik's approval.

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
  its 8-sample grid (`refine_delay`, used when `CatalogueDetections` is given
  a `refine_reader`; off by default): 3-5 times more precise
  for the coded pulses; the plasma is unchanged. Memo 29 (9c3817b).
- `validate_satellite_hours`: the satellite columns over whole hours on both
  receivers. Memo 29 (2b968f5, merge b3c4879).
- `validate_impulses_rd`: power-line impulses in the mode 300 range-Doppler
  path; they add noise, not bias. Memo 26 (ba37540).
- `millstone_radar_state`: MISA's pointing per pulse, and the zenith beam's
  pointing from satellite matches. Memo 7 (70429cd).

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
- The range-Doppler outlier blanking now reaches the last range gate. In
  mode 800 the outside line interference survived in that gate at up to
  5e4 σ, and its re-detection blanked the 16 gates below it over and over
  (ebe3be7).
- `avg_type="median"` in the range-Doppler averaging works. The median is
  not an estimate of the mean power: in mode 300 noise its level is 0.688 of
  the outlier average (ln 2 = 0.693 for exponentially distributed power), so
  absolute levels differ while ratios such as the SNR do not (385095b).
- The corner weights of the ion-line interpolation: the interpolated
  spectra were wrong by up to 24 %, now below 2.5 % (bcaaaee).
- `fit_acf_ts` evaluates the model and Jacobian at the fitted composition,
  and keeps it: in the commit's test, 84 of 84 topside gates are fitted,
  where 50-59 were (b58f1f9).

### Added
- Impulse blanking of power-line sparks (`blank_impulses`), off by default:
  impulses locked to the 60 Hz mains inflate the ACF error variance by up to
  32 %. Memo 21 (88940ed, merge 9fc6540).
- `raw_reader`: whole-second reads and a cached pulse index, 25-30 times
  faster than reading pulse by pulse in one process, up to about 40 times
  with 4-16 readers. `outlier_lpi.lpi_files` now reads through it by
  default (`fast_read=True`); the output is bit-identical, and the inversion
  gains little, as reading is about 1 % of its time. Memo 19 (10a46c3).
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
  inversion, with the injection and real-transit tests. Memos 12, 14-16
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
  interpolation swept through zero at every switch and dropped pulses (93.1 %
  of sampled pulses accepted instead of 96.3 %; about 3 % more pulses on
  reanalysis, all on misa-l) (0330e30).
- `acfs_var` is read as the complex variance it is, and `var_scale` is
  kept; every uncertainty was overstated by sqrt(2) (e2cfc6b).
- Moved off `read_vector_c81d`, which digital_rf 3 removes (c75fd76).

### Changed
- The space object cut is a significance, not an absolute number (61c79aa),
  and its blanking width is derived from the gate and pulse length instead
  of hard-coded at ±17 gates: mode 300 ±16, mode 800 ±67 (66cee68).
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
- The plasma parameter uncertainties repaired and completed: all of them
  shrink, dvi by 43 %, and the fitted parameters do not change (95aba03),
  and
  carried into the revived `fit_lp.py` (0ebc5db).

### Changed
- The delay between `tx-h` and each echo channel is measured instead of
  hard-coded: 12.34 µs (misa-l) and 11.88 µs (zenith-l) in the eclipse
  recording (ffe0e35).

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
