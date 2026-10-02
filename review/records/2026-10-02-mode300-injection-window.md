# Review record: `mode300-injection-window`

- **Gates:** A and B (product-changing)
- **Author:** Claude (Opus 5.5) sessions for Henrik: the fix on 2026-09-30, the gate A review fixes on 2026-10-01
- **Reviewer:** `isr-code-reviewer` (Sonnet 5.5, a fresh subagent), 2026-10-01; report in `documents/documents_logs/reviews/2026-10-01-gate-a-code-reviews/mode300-injection-window.md`
- **Commits:** `8e9184e`, `881bde0` (with a merge of main 8206c84: `a150107`); based on `main` at `8206c84`. Main is now `8ce21f2`, which differs only in `.claude/agents/`; merge main before the pull request.

## What the change does

Mode 300's noise-injection window started at 7800 µs, 13 µs before the noise diode switches on (about 7813 µs), so it took in the diode's switch-on overshoot: MISA's mode-300 T_sys was 26 % too low and its range-Doppler densities too low by 1.29-1.37 (Memo 27). The window now starts at 7830 µs (`radar_timing.TMM[300]["noise0"]`).

## Tests

`python3 -m pytest -q` (with the ion-line tables linked into `data/`): 25 passed. New: `test_radar_timing.py` (the window starts after the diode has settled, about 7822 µs, holds at least 500 samples, and every mode's windows lie inside the read; fails on main's 7800).

## Benchmark

`881bde0`: exit 1, 0 of 79 files bit-identical; report `~/isr_project/regression/reports/881bde0-881bde0_vs_8206c84-8206c84.md`; the same differences as `a150107` before the review fixes (a comment and a test only). Every file changes because α and T_sys change everywhere: misa-l LPI T_sys +5.7 to +8.1 % and α -5.5 to -7.5 %; misa-l range-Doppler T_sys and spectra +26 to +30 %, α -21 to -23 %; zenith-l up to 1.4 % (range-Doppler) and 0.55 % (LPI). The ACFs themselves are unchanged. Two effects need the record: a few poorly constrained zenith fit gates move between fit solutions under a 0.1-0.4 % change of α (Sigma up to 1.1e3 relative, medians 1e-8 to 4e-3), and one debris detection in `range_doppler_300_outlier/zenith-l/pp-1712592120` is lost (12 to 11), because `fit_lp`'s debris test runs on α-scaled power and that one sat at the threshold.

Together with the other six on branch `integration-all-fixes` (3fbc1ef, report `~/isr_project/regression/reports/3fbc1ef-3fbc1ef_vs_8206c84-8206c84.md`): exit 1, 0 of 81 files bit-identical, as expected. Every combined change is the sum of the single branches' changes, except one interaction: misa-l period 10875 (`pp-1712592750`), which contains an antenna switch. `antenna-switch-gap` alone loses that fit (108 of 108 gates); with `lpi-unconstrained-gates` it is recovered, and n_e changes by a median of 35 % (fewer pulses, the wrong-antenna ones removed). That is why those two, `fit-lpi-range-avg` and `fit-lpi-last-group` must be merged together.

## Scientific review (gate B, and memos with results)

Not done yet. Gate B waits for the supervisor's decision on the fix branches (TODO.tex, Q8), and is an L task (Henrik's go-ahead first).

## Full verification (when flagged, and always for calibration)

Required: this is a calibration change (REVIEW_PROCESS.md). Not done yet; part of gate B.

## Code review

- **No record, no guard test:** this file; `test_radar_timing.py` in `881bde0`.
- **The numbers disagreed** between the commit, the code comment and Memo 27 (18 vs 33 times; 29-37 % vs the benchmark's 26-30 %): the comment now quotes Memo 27's measure (the mean over the first 5 µs, 33 times on misa-l, 3 on zenith-l), and this record quotes the benchmark's ranges above. The commit message of `8e9184e` keeps its own (Memo 27's four other periods).
- **Residual spike inside the window:** the receiver filters smear the spike, so about 1 % of injection power remains on misa-l at 7830 in the range-Doppler path (none from about 7850; LPI below 2e-4). Noted in the code comment. **Open decision** (TODO.tex): 7830 or about 7850.
- **The non-default 2021 configuration** (`pass_band=0.018e6`, `filter_len=100`) would leave about 2.8 % at 7830: noted here; it is not in any config.
- **The catalogue's own mode table keeps 7800 on purpose** (its quiet window ends before the diode): stated in the comment so nobody synchronises the two.
- **Stored products** made with the old window must be marked out of date after the merge (REVIEW_PROCESS.md section 3 step 7; the Memo 8 redo).
- **Main has moved** (8ce21f2, `.claude` only): merge main before the pull request.

## Decision

- Gate A: passed on 2026-10-02 (code review answered, tests pass, every benchmark change explained or handed to gate B). Not merged: the change is product-changing.
- Gate B: waiting for the supervisor's decision on the fix branches (Q8), then the scientific review, then Henrik's approval (and Juha's, for calibration or physics).
