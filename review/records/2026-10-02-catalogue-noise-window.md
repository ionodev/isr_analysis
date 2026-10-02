# Review record: `catalogue-noise-window`

- **Gates:** A and B (product-changing)
- **Author:** Claude (Opus 5.5) sessions for Henrik: the fix on 2026-09-30, the gate A review fixes on 2026-10-01
- **Reviewer:** `isr-code-reviewer` (Sonnet 5.5, a fresh subagent), 2026-10-01; report in `documents/documents_logs/reviews/2026-10-01-gate-a-code-reviews/catalogue-noise-window.md`
- **Commits:** `8d7afb3`, `c423963` (with a merge of main 8206c84: `6777386`); based on `main` at `8206c84`. Main is now `8ce21f2`, which differs only in `.claude/agents/`; merge main before the pull request.

## What the change does

The single-pulse satellite catalogue took each pulse's noise from the 500 samples before the mode's `noise0`, where in the coded modes and mode 800 the noise diode is already on: noise 1.32-1.35 times (coded) and 1.17 times (mode 800) too high, SNRs 1.3 and 0.7 dB too low (Memo 28). The window (`quiet_window`) now ends at `last_echo`.

## Tests

`python3 -m pytest -q` (with the ion-line tables linked into `data/`): 25 passed. The noise-window test now checks `quiet_window` for every mode and decimation factor 1-16 (ends before the diode, whole decimated blocks only) and fails if the window goes back to `noise0`; a second test checks the mode table against `radar_timing.TMM`.

## Benchmark

`c423963`: exit 0, 79 of 79 files bit-identical; report `~/isr_project/regression/reports/c423963-c423963_vs_8206c84-8206c84.md`. This says nothing about the change: the benchmark does not run the catalogue. The change is product-changing (every detection's `noise_power`, `matched_filter_power`, `snr_db` and `cfar_ratio`; `doppler_hz` slightly; mode 300 too, its window moving from 7300-7800 to 7200-7700), so it needs gate B.

Together with the other six on branch `integration-all-fixes` (3fbc1ef, report `~/isr_project/regression/reports/3fbc1ef-3fbc1ef_vs_8206c84-8206c84.md`): exit 1, 0 of 81 files bit-identical, as expected. Every combined change is the sum of the single branches' changes, except one interaction: misa-l period 10875 (`pp-1712592750`), which contains an antenna switch. `antenna-switch-gap` alone loses that fit (108 of 108 gates); with `lpi-unconstrained-gates` it is recovered, and n_e changes by a median of 35 % (fewer pulses, the wrong-antenna ones removed). That is why those two, `fit-lpi-range-avg` and `fit-lpi-last-group` must be merged together.

## Scientific review (gate B, and memos with results)

Not done yet. Gate B waits for the supervisor's decision on the fix branches (TODO.tex, Q8), and is an L task (Henrik's go-ahead first).

## Full verification (when flagged, and always for calibration)

Not decided yet: the scientific review says whether it is needed.

## Code review

- **No record:** this file.
- **The test did not test the change:** fixed in `c423963`.
- **"Detections unchanged" holds in expectation, not exactly:** the DC offset now comes from a different window, so a cell at the CFAR threshold can flip. Gate B can compare detection counts on a sample, including mode 800.
- **Stored catalogue and scripts drift** (`validate_satellite_real.py`'s energy ratio) until the catalogue is rerun: decision open (rerun, or correct the stored SNRs), TODO.tex.
- **Pre-existing:** the CFAR search still ends at `noise0`, so the diode can sit in the training cells of the last valid rows; not a regression, out of scope.

## Decision

- Gate A: passed on 2026-10-02 (code review answered, tests pass, every benchmark change explained or handed to gate B). Not merged: the change is product-changing.
- Gate B: waiting for the supervisor's decision on the fix branches (Q8), then the scientific review, then Henrik's approval (and Juha's, for calibration or physics).
