# Review record: `lpi-unconstrained-gates`

- **Gates:** A and B (product-changing)
- **Author:** Claude (Opus 5.5) sessions for Henrik: the fix on 2026-09-30, the gate A review fixes on 2026-10-01
- **Reviewer:** `isr-code-reviewer` (Sonnet 5.5, a fresh subagent), 2026-10-01; report in `documents/documents_logs/reviews/2026-10-01-gate-a-code-reviews/lpi-unconstrained-gates.md`
- **Commits:** `a999223`, `33022a1` (with a merge of main 8206c84: `0877912`); based on `main` at `8206c84`. Main is now `8ce21f2`, which differs only in `.claude/agents/`; merge main before the pull request.

## What the change does

With the default outlier rejection, the rejection can empty the clutter gates of a lag; that range gate is then an unknown no measurement constrains, the normal matrix is singular, and the whole lag was lost (36 of 78 MISA periods in a test hour, a median of 11 of 45 lags; Memo 29). `invert_normal` inverts for the constrained unknowns and leaves the unconstrained ones NaN; with every unknown constrained it is the plain inverse, as before.

## Tests

`python3 -m pytest -q` (with the ion-line tables linked into `data/`): 26 passed. New: `test_outlier_lpi.py` (the solve against least squares on the constrained columns; NaN exactly at the zero columns; bit-identical to `inv` when all are constrained; a NaN in the matrix is not taken for unconstrained).

## Benchmark

`33022a1`: exit 1, 74 of 79 files bit-identical; report `~/isr_project/regression/reports/33022a1-33022a1_vs_8206c84-8206c84.md`. The same differences as `0877912` before the review fixes (`lpi-unconstrained-gates-0877912_vs_main-8206c84.md`), so the refactor is output-neutral. All 5 changed files are misa-l: three LPI files get back every lost lag (8, 11 and 12 lags; every value main had is unchanged, only NaN became finite), and two fits change because of that. Period 2's fit (`pp-1712484020`) had no n_e on main; on the branch it fits at 107 of 108 gates, but 4 gates sit at the Te bound (9000 K), 7 at Ti (3000 K) and 9 at |v_i| = 1500 m/s, mostly at the top: finite is not the same as usable (for gate B).

Together with the other six on branch `integration-all-fixes` (3fbc1ef, report `~/isr_project/regression/reports/3fbc1ef-3fbc1ef_vs_8206c84-8206c84.md`): exit 1, 0 of 81 files bit-identical, as expected. Every combined change is the sum of the single branches' changes, except one interaction: misa-l period 10875 (`pp-1712592750`), which contains an antenna switch. `antenna-switch-gap` alone loses that fit (108 of 108 gates); with `lpi-unconstrained-gates` it is recovered, and n_e changes by a median of 35 % (fewer pulses, the wrong-antenna ones removed). That is why those two, `fit-lpi-range-avg` and `fit-lpi-last-group` must be merged together.

## Scientific review (gate B, and memos with results)

Not done yet. Gate B waits for the supervisor's decision on the fix branches (TODO.tex, Q8), and is an L task (Henrik's go-ahead first).

## Full verification (when flagged, and always for calibration)

Not decided yet: the scientific review says whether it is needed.

## Code review

- **No test** of the solve: fixed in `33022a1` (helper `invert_normal`, `test_outlier_lpi.py`).
- **No review record:** this file.
- **A NaN diagonal counted as unconstrained** (`NaN > 0` is false), which would solve around NaN-contaminated rows silently: fixed in `33022a1`, only an exactly zero diagonal counts.
- **The benchmark summary's "no finite value at all" for period 2** was wrong (main held the fit's start values at one gate): corrected in the summary.
- **Log volume** (one line per affected lag and period): noted, left.

## Decision

- Gate A: passed on 2026-10-02 (code review answered, tests pass, every benchmark change explained or handed to gate B). Not merged: the change is product-changing.
- Gate B: waiting for the supervisor's decision on the fix branches (Q8), then the scientific review, then Henrik's approval (and Juha's, for calibration or physics).
