# Review record: `fit-lpi-range-avg`

- **Gates:** A and B (product-changing)
- **Author:** Claude (Opus 5.5) sessions for Henrik: the fix on 2026-09-30, the gate A review fixes on 2026-10-01
- **Reviewer:** `isr-code-reviewer` (Sonnet 5.5, a fresh subagent), 2026-10-01; report in `documents/documents_logs/reviews/2026-10-01-gate-a-code-reviews/fit-lpi-range-avg.md`
- **Commits:** `446d239` (with a merge of main 8206c84: `f52508e`); based on `main` at `8206c84`. Main is now `8ce21f2`, which differs only in `.claude/agents/`; merge main before the pull request.

## What the change does

In `fit_lpi`'s range average the variance slice lacked `max(0, ...)`, so the lowest `ra` gates got infinite variance, and both slices averaged 2 ra gates, half a gate off centre, instead of the stated 2 ra + 1. Both now use `[max(0, ri-ra), min(n, ri+ra+1))`.

## Tests

`python3 -m pytest -q` (with the ion-line tables linked into `data/`): 23 passed (`test_fit_lpi_ts.py` with the tables linked). No new test: the slice sits inside `fit_lpifiles` (finding below).

## Benchmark

`f52508e`: exit 1; report `~/isr_project/regression/reports/fit-lpi-range-avg-f52508e_vs_8206c84-8206c84.md`. The 8 `fit_lpi` files change and the other 71 are bit-identical. Each fit gains its lowest gate (finite where it was infinite variance), and parameters move with the recentred window (zenith-l n_e by a median relative 0.16-0.20 in these 60-s fits, error bars 1.3-7.7 % smaller). The empty fit of misa-l period 2 changes v_i from +100 to -100 m/s: explained, see below. The commit message's numbers (85 km gate, 0.2-0.7 σ) come from an earlier test against ba37540, not from the benchmark.

Together with the other six on branch `integration-all-fixes` (3fbc1ef, report `~/isr_project/regression/reports/3fbc1ef-3fbc1ef_vs_8206c84-8206c84.md`): exit 1, 0 of 81 files bit-identical, as expected. Every combined change is the sum of the single branches' changes, except one interaction: misa-l period 10875 (`pp-1712592750`), which contains an antenna switch. `antenna-switch-gap` alone loses that fit (108 of 108 gates); with `lpi-unconstrained-gates` it is recovered, and n_e changes by a median of 35 % (fewer pulses, the wrong-antenna ones removed). That is why those two, `fit-lpi-range-avg` and `fit-lpi-last-group` must be merged together.

## Scientific review (gate B, and memos with results)

Not done yet. Gate B waits for the supervisor's decision on the fix branches (TODO.tex, Q8), and is an L task (Henrik's go-ahead first).

## Full verification (when flagged, and always for calibration)

Not decided yet: the scientific review says whether it is needed.

## Code review

- **No record:** this file.
- **`fit_ionline.py` has the same bug** (a verbatim copy, lines 710-722): it looks unused; TODO.tex.
- **The ±100 m/s:** a pre-existing bug in `fit_acf`. With lag 0 lost the fit returns its start guess, and `fit_acf` flips the caller's guess in place, so a gate that raises leaves the next flat gate at -100 m/s. Not caused by this branch; TODO.tex (fix: copy the guess, NaN without lag 0). `lpi-unconstrained-gates` removes the case here.
- **The variance of the mean** (1/Σ(1/v)) does not match the r²-weighted ACF mean, and a NaN gate still counts in the ACF mean's weights: pre-existing; decision open (TODO.tex).
- **No unit test** of the slice: it would need the averaging pulled into its own function; left for now.
- **Edge windows** hold fewer than 2 ra + 1 gates while `range_avg_window_km` stays nominal; the default `range_avg` does not average below 300 km.
- **Merge order:** `fit-lpi-correlated-gates` has the same slice edit; no conflict.

## Decision

- Gate A: passed on 2026-10-02 (code review answered, tests pass, every benchmark change explained or handed to gate B). Not merged: the change is product-changing.
- Gate B: waiting for the supervisor's decision on the fix branches (Q8), then the scientific review, then Henrik's approval (and Juha's, for calibration or physics).
