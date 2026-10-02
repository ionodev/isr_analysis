# Review record: `tx-delay-own-antenna`

- **Gates:** A and B (product-changing)
- **Author:** Claude (Opus 5.5) sessions for Henrik: the fix on 2026-09-30, the gate A review fixes on 2026-10-01
- **Reviewer:** `isr-code-reviewer` (Sonnet 5.5, a fresh subagent), 2026-10-01; report in `documents/documents_logs/reviews/2026-10-01-gate-a-code-reviews/tx-delay-own-antenna.md`
- **Commits:** `784f00b`, `c908db7` (with a merge of main 8206c84: `e43b5fd`); based on `main` at `8206c84`. Main is now `8ce21f2`, which differs only in `.claude/agents/`; merge main before the pull request.

## What the change does

The channel-delay estimator looked 20 s ahead for coded pulses on the receiver's own antenna and otherwise used the other antenna's, whose leakage arrives 1.2 µs (zenith-l) or 2.3 µs (misa-l) away (Memo 29). The window now doubles, up to `max_search_s` = 2400 s, until it holds `n_pulses` own-antenna pulses; only then does it fall back.

## Tests

`python3 -m pytest -q` (with the ion-line tables linked into `data/`): 27 passed. New: `test_tx_delay.py` with fake readers (own pulses in the first window; a window that has to grow, which fails on main; the cap and the fallback; a recording shorter than the window).

## Benchmark

`c908db7`: exit 1, 64 of 79 files bit-identical; report `~/isr_project/regression/reports/c908db7-c908db7_vs_8206c84-8206c84.md`; the same differences as `e43b5fd` before the review fixes. zenith-l's delay moves from 11.88 to 10.62 µs (the recording starts in MISA's cycle, so main fell back to MISA's pulses), which shifts every zenith-l range-Doppler spectrum by 1.3 samples (195 m) and moves the debris mask by 2 gates. misa-l is bit-identical (its first 20 s hold own pulses). **Not explained:** single gates change by up to 12-23 σ; and the pipeline's 10.62 µs differs from an independent 10.72 µs (`verify_combined_txdelay.txt`). Both are handed to gate B.

Together with the other six on branch `integration-all-fixes` (3fbc1ef, report `~/isr_project/regression/reports/3fbc1ef-3fbc1ef_vs_8206c84-8206c84.md`): exit 1, 0 of 81 files bit-identical, as expected. Every combined change is the sum of the single branches' changes, except one interaction: misa-l period 10875 (`pp-1712592750`), which contains an antenna switch. `antenna-switch-gap` alone loses that fit (108 of 108 gates); with `lpi-unconstrained-gates` it is recovered, and n_e changes by a median of 35 % (fewer pulses, the wrong-antenna ones removed). That is why those two, `fit-lpi-range-avg` and `fit-lpi-last-group` must be merged together.

## Scientific review (gate B, and memos with results)

Not done yet. Gate B waits for the supervisor's decision on the fix branches (TODO.tex, Q8), and is an L task (Henrik's go-ahead first).

## Full verification (when flagged, and always for calibration)

Not decided yet: the scientific review says whether it is needed (it has the two open points above).

## Code review

- **No test, no record:** `test_tx_delay.py` in `c908db7`; this file.
- **Wrong-antenna pulses at the end of an own cycle** still pass the antenna test (the metadata lag): this is what `antenna-switch-gap` removes, so the two belong together. Probably the two exceptions of Memo 29.
- **Benchmark coverage:** the cap, the fallback, mode 800 and other recordings are not exercised by the benchmark; the unit test now covers the first two.
- **Units:** the first window was a microsecond count, the cap in seconds: `c908db7` uses seconds (`n_pulses/5` s, 20 s for 100) and clips it to the cap.
- **Cost:** each doubling rescreened every pulse: `c908db7` screens only the new ones and stops at `n_pulses`.
- **Callers** take the new default silently, and stored validation results (Memo 29's hour run) will not reproduce exactly: noted.
- **Provenance** of "118 of 120" in the commit message: from Memo 29's run, not the benchmark.

## Decision

- Gate A: passed on 2026-10-02 (code review answered, tests pass, every benchmark change explained or handed to gate B). Not merged: the change is product-changing.
- Gate B: waiting for the supervisor's decision on the fix branches (Q8), then the scientific review, then Henrik's approval (and Juha's, for calibration or physics).
