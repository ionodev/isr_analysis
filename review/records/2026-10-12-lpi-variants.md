# Review record: `lpi-variants`

- **Gates:** A only (output-neutral)
- **Author:** Claude session of 12 October 2026, for Henrik
- **Reviewer:** `isr-code-reviewer` subagent (Sonnet 5.5), two rounds
- **Commits:** `fc0eb49..7e5d5a8`, based on `main` at `7f44991`

## What the change does

The weeks 7-9 comparison runs five LPI variants on the same data (Henrik,
12 October): (1) equal weights (ordinary least squares), (2) variance-weighted,
(3) variance-weighted with the background column, (4) variance-weighted with
the satellite columns, (5) with both. `outlier_lpi.lpi_files` gets two options:
`weighting` (`"variance"`, the default as before, or `"equal"`: one common
noise level for every lagged product, with the error covariance from each
product's own level, `ols_covariance`) and `background_column` (`True`, the
default as before, or `False`: no column of ones). `run_analysis.py` reads
these, `outlier_rejection` and `periods` from the lpi configuration, and a new
`satellites` section (catalogue, channel delay given or measured, minimum SNR)
makes it model the catalogued satellite echoes: the satellite handling as a
pipeline stage (Henrik, Q4 and 12 October).

## Tests

`python3 -m pytest -q`: 61 passed. New: `test_lpi_variants.py` (the two
weightings, the background column on and off, equal weights give the OLS
estimate, the OLS covariance against the direct formula with an unconstrained
unknown) and `test_satellite_stage.py` (no section, disabled section, given
delay, measured delay).

End to end (not part of the tests): all five variants on period 3070
(STARLINK-31553, zenith-l) with the catalogue of 2 October; variant 5 through
`run_analysis.py` gives the same file as the direct call (895 satellite
columns, channel delay 10.62 us).

## Benchmark

`python3 ~/isr_project/isr_analysis/review/regression.py lpi-variants`: exit
status 2. Report:
`~/isr_project/regression/reports/lpi-variants-7e5d5a8_vs_main-7f44991.md`.

80 output files compared; **every output is bit-identical**. The status is 2
only because the tool expects `lpi_60/zenith-l/lpi-1712594080.h5` (period
11008, made only of wrong-antenna pulses), which `main` rightly no longer
makes since pull request #4 (`antenna-switch-gap`); the same gap as in
`2026-10-12-apply-review-decisions.md`. 0 failed jobs. The tool's expected
list needs updating (a change to `review/`, which needs Henrik's approval).

## Summary for Henrik's review (gate B)

Not applicable: output-neutral.

## Scientific review (only if Henrik asked for one)

Not requested.

## Code review

First round (on fc0eb49, d579402): passes after fixes.
- F1 (fix): with equal weights, `acfs_var` was the nominal covariance
  `(B^H B)^-1`, a lower bound, too low most at bright gates, short lags and
  satellites, so variant 1's error bars and its weights in `fit_lpi` would not
  be comparable with the others. Fixed in 07963e6 (`ols_covariance`, the
  sandwich covariance); on period 3070 the estimate is unchanged and the lag-1
  sigma goes from a flat 357 to 1.6e3 far from the echo and up to 1e8 at it.
- F2 (fix): no record or changelog entry. Changelog in 7e5d5a8; this record.
- F3 (caveat): the variants are not separated on disk, and the skip test
  looks only at the png: give each variant its own `output_dir` (the run
  configurations in `~/isr_project/variants/configs/` do).
- F4 (caveat): `outlier_rejection` is not tied to the satellite stage; each
  variant's configuration sets it explicitly (off in all five, so that they
  differ only in weighting and columns).
- F5 (caveat): with `channel_delay_us` null every MPI rank measures the delay
  from the start of the recording; the configurations give it (10.62 us for
  zenith-l).
- F6 (caveat): the tests cover the helpers, not the slicing of the unknowns
  inside `lpi_files`; the reviewer checked the indexing by reading, in all four
  combinations.
- F7 (drop): equal mode does not exclude rows with a NaN local sigma; cannot
  occur (no NaN before rejection).

Second round (on 07963e6, 7e5d5a8): passes. F1 and F2 fixed; F3-F6 accepted as
caveats, F7 dropped.

## Changelog

Under 2026-10-12, Added (7e5d5a8): the LPI variants, and the satellite stage
of the lpi step. Moved to the merge date if that differs.

## Decision

- Gate A: passed on 2026-10-12, waiting for Henrik's review of the pull request.
