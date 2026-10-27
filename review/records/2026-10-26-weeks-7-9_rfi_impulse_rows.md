# Review record: `weeks-7-9_rfi_impulse_rows`

- **Gates:** A (output-neutral with the switch off, its default). Turning
  the switch on for a product would be gate B (Henrik).
- **Author:** Claude (Opus 5.5), autonomous run of 25-26 October 2026
  (queue `documents_logs/overnight/plan-2026-10-25-night.md`, item 1).
- **Reviewer:** `isr-code-reviewer` subagent (Sonnet 5.5), 26 October 2026.
- **Commits:** `7173173` (the change) and `64cc604` (the measure script
  `compare_impulse_rows.py`), based on `weeks-7-9` at `b64be39`.

## What the change does

RFI mitigation as Juha recommended (Henrik, 25 October): the lagged products
that impulsive interference touches are removed from the lag profile
inversion. `impulse_rows.py` detects impulses in each pulse's raw echo window
(power over the running median of 65 samples above 30, two guard samples:
the test of Memo 21's blanking, with the median mirrored at the window ends),
widens them by the filter span (+-filter_len) and marks the rows of each lag
whose products at t or t+lag touch them. `outlier_lpi.lpi_files(
impulse_rows=False)` removes those rows (of the pulse and of its
clutter-subtraction partner) from the measurement, the theory matrix and
any satellite column, takes the variance weights' running mean over the kept
rows only, and writes the removed fraction per lag. `run_analysis.py` maps
`"impulse_rows"`. Memo 47 has the method and the test on spark minutes.

## Tests

`python3 -m pytest -q`: 68 passed, nothing skipped (tables linked in). New:
`test_impulse_rows.py` (detection finds an impulse and not an extended echo;
the widening equals the real filter's reach; kept rows carry no trace of a
synthetic impulse and only rows near it are removed; the lag-10 ACF of a
correlated process is unbiased after removal and biased without; the masked
running mean equals the original with nothing masked).

## Benchmark

`python3 ~/isr_project/isr_analysis/review/regression.py
weeks-7-9_rfi_impulse_rows --base weeks-7-9 --jobs 16` (main's copy of the
tool; base is the period branch, as the branch contains it, not `main`):
exit status 0, "Every output is bit-identical": 80 output files compared,
80 bit-identical (run 2026-10-26 04:45 UTC, environment key 23efb11a80e8).
Report: `~/isr_project/regression/reports/weeks-7-9_rfi_impulse_rows-7173173_vs_weeks-7-9-b64be39.md`.

`64cc604` adds an analysis script that no pipeline stage imports; the outputs
cannot depend on it, so the benchmark was not repeated for it.

## Scientific review

Not requested.

## Code review

Verdict: passes, nothing must be fixed. Notes, answered:
- Stale partner pulse when a pulse's code matches neither neighbour (an
  existing behaviour of `lpi_files`, not changed here): with the switch on,
  impulses in that stale partner are not looked for. Accepted: rare, and
  only that pulse.
- Possible false detections at a steep power edge such as the ground clutter
  at the start of the window (seen on a synthetic 10-sample e-fold): checked
  on real periods with the switch on, the retained fraction at the first
  usable delay (1020 us) equals that of the rest of the window (0.92 against
  0.91 and 0.74 against 0.69), so no excess removal at the clutter edge.
- `impulse_sample_fraction` counts the union of the two pulses' impulse
  samples over both pulses' samples, so it is slightly low; diagnostic only.
- No integration test runs `lpi_files` with the switch on; covered by the
  runs of Memo 47 (71 periods).
- `compare_impulse_rows.py` was uncommitted at review time; committed in
  `64cc604`.

## Changelog

To be added with the delivery pull request of weeks 7-9, under "Added":
the `impulse_rows` switch (off by default), Memo 47, this record.

## Decision

- Gate A: tests and code review passed on 2026-10-26; benchmark as above.
  Not merged (Henrik's decision for this run: no merges, no pull requests);
  the branch is pushed for his review.
