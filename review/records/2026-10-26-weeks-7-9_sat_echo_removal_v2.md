# Review record: `weeks-7-9_sat_echo_removal_v2`

- **Gates:** A (every new option is off by default; runs without them are
  unchanged). Using the lower threshold or outlier rejection for products is
  Henrik's decision (gate B if it becomes a default).
- **Author:** Claude (Opus 5.5), autonomous run of 25-26 October 2026
  (queue `documents_logs/overnight/plan-2026-10-25-night.md`, items 3-7).
- **Reviewer:** `isr-code-reviewer` subagent (Sonnet 5.5), 26 October 2026.
- **Commits:** `5a6a298..5c88bac`, based on `weeks-7-9_sat_echo_removal_v1`
  at `31025ca`.

## What the change does

Options F, A and C of Memo 46 and an outlier-rejection test (Memo 48):
`RandomControlColumns` (pure-plasma control); `cfar_alpha` and `--channels`
in the detector; strong/weak rows in `CatalogueDetections(strong_pfa)` and
the tracking confirmation `TrackedDetections(min_weak_pulses)`;
`column_residuals` (`sat_resid`) in the LPI and the RESIDUAL flag in
`fit_quality`; with outlier rejection and satellite columns, what each test
removes (`outlier_*_fraction`, `sat_rows_rejected`); end-to-end test at a
given P_FA (`realistic_confirmed`); analysis scripts (`random_control.py`,
`column_checks.py`, `confirmation_count.py`, `interval_measures.py`,
`rejection_measures.py`, `e2e_threshold_figure.py`).

## Tests

`OMP_NUM_THREADS=1 python3 -m pytest -q`: 84 passed, nothing skipped (tables
linked in). New: `test_random_columns.py`, `test_confirmation.py` (incl.
`cfar_alpha` equal to the detector's threshold), `test_column_residuals.py`
(incl. the RESIDUAL flag).

## Benchmark

`python3 ~/isr_project/isr_analysis/review/regression.py
weeks-7-9_sat_echo_removal_v2 --base weeks-7-9_sat_echo_removal_v1 --jobs 16`
(main's copy of the tool), started 26 October about 03:35 UTC on `68cc8fb`;
still running at the end of the night (the machine was overloaded): its
result is in `/mnt/data/henrik/sat_v2_2026-10/benchmark_report.txt` and
`~/isr_project/regression/reports/`. The commit after it (`5c88bac`) changes
`run_analysis.py` (not covered by the benchmark; shown neutral by
`test_satellite_stage.py`, which it fixes) and docstrings only.
Independent evidence of neutrality: the control base run (satv2 code, tracked
columns, the quiet interval) reproduced the v1 tracked output bit for bit in
all 18 periods (`acfs_e`, `acfs_g`, `acfs_var`, `sat_amp_e`, `sat_amp_var`).

## Scientific review

Not requested.

## Code review

Verdict: passes after one fix.
- Must-fix: `run_analysis.py` passed `strong_pfa` to `CatalogueDetections`
  always, which broke the stub in `test_satellite_stage.py` (2 tests failed).
  Fixed in `5c88bac` (passed only when configured); reviewer confirmed.
- Notes, answered: `sat_resid` is residual over the row's local power, which
  includes the echo, so an echo left on a long pass is not flagged
  (docstrings corrected in `5c88bac`; Memo 48 shows the flag does not
  discriminate); `--null` documented but absent (removed); a column with a
  significant echo but no usable long lags counts as coherent (edge case,
  accepted); a weak detection linked to a strong pass gets a column (by
  design; precision only, by the control); random columns clear in both
  delay and Doppler (docstring); the control tests the nearest gate only
  (a limitation in Memo 48); `sat_rows_rejected` counts rows already NaN as
  not rejected (slightly low; accepted).

## Changelog

To be added with the delivery pull request of weeks 7-9, under "Added"
(Memo 48, this record).

## Decision

- Gate A: tests and code review passed on 2026-10-26; benchmark pending.
  Not merged (Henrik's decision for this run: no merges, no pull requests);
  the branch is pushed for his review.
