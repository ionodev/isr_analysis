# Review record: `antenna-switch-gap`

- **Gates:** A and B (product-changing)
- **Author:** Claude (Opus 5.5) sessions for Henrik: the fix on 2026-09-30, the gate A review fixes on 2026-10-01
- **Reviewer:** `isr-code-reviewer` (Sonnet 5.5, a fresh subagent), 2026-10-01; report in `documents/documents_logs/reviews/2026-10-01-gate-a-code-reviews/antenna-switch-gap.md`
- **Commits:** `54eb8ea`, `8c15038`, and at gate B `bcf989f` (with merges of main: 8206c84 in `8fc50c1`, a7e3cb5 in `2dc023f`); based on `main` at `a7e3cb5`.

## What the change does

The antenna metadata change antenna at the event that opens the next record cycle, but the radar changes at the end of the previous one, a median of 8.6 s earlier (2.3-20.8 s; Memo 30), so a channel's inversion took in pulses sent on the other antenna or not at all. `get_antenna_select` now returns 0 (unknown) from `switch_guard_s` = 3.0 s before each closing event until the opening event, and every caller rejects 0. The guard was chosen on this recording: 2.5 s until gate B, 3.0 s since `bcf989f` (see the scientific review).

## Tests

`python3 -m pytest -q` (with the ion-line tables linked into `data/`): 29 passed, also after `bcf989f` (the tests now use the default guard). New: `test_millstone_radar_state.py` (a normal cycle, a recording that starts inside a cycle, the last event, the clamp on a short cycle, a change without a closing event, and that the values stay steps in {-1, 0, 1}); the start-of-recording test fails on `54eb8ea`.

## Benchmark

`bcf989f` (guard 3.0 s): not run yet; it needs a new benchmark before the merge. The numbers below are for 2.5 s.

`8c15038`: exit 1, 72 of 79 files bit-identical; report `~/isr_project/regression/reports/8c15038-8c15038_vs_8206c84-8206c84.md`. Before the fix (`8fc50c1`) it was 70 of 79: misa-l's first period and its fit were lost, because the first metadata event (+52.3 s) closes a cycle and its 0 was held back to the recording start. The remaining 7 files are the intended changes at antenna switches: one zenith-l 10-s period (`lpi-1712594080`) is no longer written, as none of its pulses passes the antenna test any more (Memo 30's periods made only of wrong pulses); misa-l `lpi-1712592900` loses 1183 values (pulses at the end of its cycle removed; P_tx changes by 12.6 %), and its fit `pp-1712592750` loses all 108 gates on this branch alone, which `lpi-unconstrained-gates` recovers (see the integration below).

Together with the other six on branch `integration-all-fixes` (3fbc1ef, report `~/isr_project/regression/reports/3fbc1ef-3fbc1ef_vs_8206c84-8206c84.md`): exit 1, 0 of 81 files bit-identical, as expected. Every combined change is the sum of the single branches' changes, except one interaction: misa-l period 10875 (`pp-1712592750`), which contains an antenna switch. `antenna-switch-gap` alone loses that fit (108 of 108 gates); with `lpi-unconstrained-gates` it is recovered, and n_e changes by a median of 35 % (fewer pulses, the wrong-antenna ones removed). That is why those two, `fit-lpi-range-avg` and `fit-lpi-last-group` must be merged together.

## Scientific review (gate B, and memos with results)

Reviewer: `isr-science-reviewer` (Opus 5.5, a fresh subagent), 2 October 2026. It judged the four linked branches, `lpi-unconstrained-gates`, `antenna-switch-gap`, `fit-lpi-range-avg` and `fit-lpi-last-group`, as one change. Report: `documents/documents_logs/reviews/2026-10-02-gate-b/four-linked-branches.md`; scripts in `scripts_four/` beside it.

- **Verdict: sound with changes.** Each of the four fixes does what it claims: the off-centre range window biased n_e high by about 4-5 %, the 2.5 s guard removed every pulse sent on the other antenna, and the recovered lags leave every value main already had unchanged. The required changes are in `fit-lpi-range-avg` (E1), `fit-lpi-last-group` (E2) and `antenna-switch-gap` (G1).
- **Spot checks** (report section 3):
  - The reviewer's own guard logic on the per-pulse `results/antenna_switch_lag.npz` (`guard.py`) reproduces Memo 30: at 2.5 s, 0 pulses sent on the other antenna and 288 (zenith-l) and 8 (misa-l) never sent; at 0 s, 43,887 and 16,027. **Agrees.**
  - Per-file calibrated power and P_tx (`pw.py`): the removed zenith-l file `lpi-1712594080` has power −0.077 at P_tx 448 kW, so it held only wrong pulses. **Agrees.**
- **Findings that concern this branch:**
  - **G1, the 2.5 s guard sat at the edge, and Memo 30's account of the leftover pulses is inaccurate.** Memo 30 says the never-sent pulses fall "elsewhere in the cycles". They do not: 273 come from one change, 8 April 05:51:33 UTC, where the transmitter stopped 3.85-6.28 s before the closing event, leaving one zenith-l period with about 28 % noise-only pulses; the other 23 sit at −2.57 to −2.50 s, right at the guard's start, in 5 changes. At 3.0 s: misa-l 0, zenith-l 273 (only the anomalous change), at a cost of about 0.5 s of good pulses per change. Tuning on this recording is acceptable: the guard is a measured property of the radar controller's timing from 592 changes, not fitted to the science result. *Done:* guard 3.0 s in `bcf989f` (Henrik), with the numbers in the docstring and comment. *Not done:* the reviewer's better alternative, a per-period check from `tx-h` (the fraction of accepted pulses whose transmitter power puts them in the wrong class; error-free in this recording by Memo 36), which would also catch the 05:51:33 change.
  - **G2, the "median 35 %" at period 10875** (`pp-1712592750`, the interaction in the Benchmark above). It is the median *absolute* relative change; the signed median ratio is +25 % over all gates and +28 % at 150-500 km height, 1.1σ per gate. Its parts: main took lag 0 only from `lpi-1712592900`, the file with wrong pulses, whose calibrated power rises by 1.31 and P_tx from 1109 to 1249 kW with the fix (its neighbours have 1237-1257 kW); and `mode300-injection-window` (+6.6 % in α) is in the integration too. The direction and size are plausible. *Open:* after cleaning, 592900 still has only 0.80 of 592750's power (residual contamination, MISA moving, or a real change).
- **The other findings of the joint review** (E1 range average, E2 incremental runs, G3 fit limits, G4 absolute scale) are in those branches' records; E1 and E2 are fixed there.
- **Claims the reviewer could not confirm:**
  - that every combined change is the sum of the single branches' changes, apart from period 10875 (the integration run was not decomposed);
  - Memo 30's per-file power recovery, and the noise-diode calibration inside the guarded stretch (never checked by anyone);
  - why `lpi-1712592900` (period 10890), after cleaning, still has only 0.80 of the echo power of `lpi-1712592750`;
  - a last group longer than one file (not exercised);
  - anything about the eclipse response: there are no eclipse-hour products from these branches.
- **Questions for the supervisor** (report section 4):
  1. Is the controller's timing known, from the closing event to the antenna change and the transmitter pause? Is the early transmitter stop on 8 April 05:51:33 UTC a known event? Either could replace the tuned guard.
  2. Should `fit_lpi` write values at debris-masked or empty gates at all? (For the range average, `fit-lpi-range-avg` now writes NaN there.)
  3. Fits at Te/Ti = 3: extend the tables, or flag? Can midday or dawn topside Te/Ti at Millstone exceed 3?
  4. Are incremental or near-real-time runs (`reanalyze` false) a real use case?
  5. Below 300 km the fit uses `acfs_g/2` with the variance of `acfs_e`: is that pairing intended?

**Henrik's decisions, 2 October 2026**, on the gate B reports:
- all seven fix branches are approved;
- four changes come first:
  - `fit-lpi-range-avg`: the range average counts only gates with data, uses the variance of the r²-weighted mean, and gives NaN at an empty or debris-masked centre gate;
  - `fit-lpi-last-group`: each fit stores its number of LPI files, and an incremental run refits a window that has grown;
  - `antenna-switch-gap`: a guard of 3.0 s;
  - `mode300-injection-window`: the injection window starts at 7850 µs;
- the satellite catalogue is rerun after the merge;
- a per-gate flag for fits at a fit limit (a parameter bound) comes later, in weeks 7-9.

## Full verification (when flagged, and always for calibration)

Called for by the reviewer, but scoped, and it need not block the code merge once the required changes are in. Two reasons:
- the merged set shifts the absolute densities (the range-average fix lowers zenith n_e near the F peak by about 3-5 %), and Memo 38's digisonde constant C was made from products built before the fixes (`lpi_30`, June 2024);
- nothing has been run on the eclipse hours or at night: the benchmark is one daytime hour plus a few single periods.

To do after the merge: redo Memo 38's C with the merged code (it should rise by about 3-5 %, and its 7 % spread should not grow), and fit the eclipse hours and a night hour. Not done yet.

## Code review

- **Real bug: the start of the recording counted as unknown** (60.6 s of MISA pulses dropped): fixed in `8c15038` (the step functions hold the first recorded value before the first event).
- **No test:** fixed in `8c15038`.
- **A missing `cycle_name` field** would have silently guarded on an arbitrary basis or crashed: `8c15038` falls back to the metadata as recorded, with a message.
- **Switches not preceded by a closing event** were left unguarded silently: `8c15038` reports how many.
- **The docstring's "8.2-20.8 s"** was not Memo 30's: corrected to 2.3-20.8 s (median 8.6 s).
- **A trailing closing event** would hold its antenna to the end: not present in this recording; left.
- **`tsys_harvest.py`'s comment** now documents 0.
- **The satellite catalogue reads the raw metadata itself** and keeps the 8-s lag: outside this branch, noted.
- **Interplay with `tx-delay-own-antenna`:** the delay is measured at the recording start; with the start fixed, both measure MISA's first cycle again.

## Decision

- Gate A: passed on 2026-10-02 (code review answered, tests pass, every benchmark change explained or handed to gate B). Not merged: the change is product-changing.
- Gate B: scientific review done on 2026-10-02 (above); approved by Henrik on 2026-10-02 with the guard at 3.0 s. `bcf989f` needs gate A (code review, benchmark) before the merge, together with `lpi-unconstrained-gates`, `fit-lpi-range-avg` and `fit-lpi-last-group`. The scoped full verification follows the merge.
