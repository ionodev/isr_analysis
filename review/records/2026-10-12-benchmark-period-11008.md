# Review record: `benchmark-period-11008`

- **Gates:** A only (review tooling; needs Henrik's approval through this pull request)
- **Author:** Claude session of 12 October 2026, at Henrik's request
- **Reviewer:** `isr-code-reviewer` subagent (Sonnet 5.5)
- **Commits:** `d65606e..46c5445`, based on `main` at `7f44991`

## What the change does

Since pull request #4 (`antenna-switch-gap`) the pipeline rightly makes no
LPI file for zenith-l period 11008, which holds only wrong-antenna pulses
(Memo 30). The benchmark still expected one, so every run against `main`
ended with exit status 2 (invalid), even with every output bit-identical
(records `2026-10-12-apply-review-decisions.md`, `2026-10-12-lpi-variants.md`).
`review/benchmark.json` now lists the period under `lpi_no_output`, and
`expected()` in `review/regression.py` skips it. The period is still run; a
file appearing for it on one side shows as a difference. Also corrects the
label of period 10812: it lies in the Starlink hour and is not quiet (352
catalogued echoes, catalogue of 2 October).

## Tests

No pipeline code changes. The check is a run of the changed tool (below).

## Benchmark

`python3 review/regression.py benchmark-period-11008 --allow-other-tool`
(this branch's tool, as the process allows for a change to the tool itself):
80 output files compared, **every output bit-identical**, no missing file
(exit status 0). Report:
`~/isr_project/regression/reports/benchmark-period-11008-d65606e_vs_main-7f44991.md`.
The environment key changes with `benchmark.json` (23efb11a80e8), so cached
runs are redone once.

## Summary for Henrik's review (gate B)

Not applicable: no product changes.

## Scientific review (only if Henrik asked for one)

Not requested.

## Code review

Passes, no fixes. Caveats: the CHANGELOG cites the code commit (a commit
cannot cite itself); the fit stretch check sees 29 instead of 30 files, far
above its minimum of 10. Checked: the string keys match so the skip fires,
11008 is still run, nothing else reads the new keys, the JSON is valid.

## Changelog

Under 2026-10-12, Process (46c5445).

## Decision

- Gate A: passed on 2026-10-12, waiting for Henrik's approval of the pull request.
