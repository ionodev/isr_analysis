# Review record: `catalogue-noise-window`

- **Gates:** A and B (product-changing)
- **Author:** Claude (Opus 5.5) sessions for Henrik: the fix on 2026-09-30, the gate A review fixes on 2026-10-01
- **Reviewer:** `isr-code-reviewer` (Sonnet 5.5, a fresh subagent), 2026-10-01; report in `documents/documents_logs/reviews/2026-10-01-gate-a-code-reviews/catalogue-noise-window.md`
- **Commits:** `8d7afb3`, `c423963` (with merges of main: 8206c84 in `6777386`, a7e3cb5 in `4b1eedf`); based on `main` at `a7e3cb5`. No code change at gate B.

## What the change does

The single-pulse satellite catalogue took each pulse's noise from the 500 samples before the mode's `noise0`, where in the coded modes and mode 800 the noise diode is already on: noise 1.32-1.35 times (coded) and 1.17 times (mode 800) too high, SNRs 1.3 and 0.7 dB too low (Memo 28). The window (`quiet_window`) now ends at `last_echo`.

## Tests

`python3 -m pytest -q` (with the ion-line tables linked into `data/`): 25 passed. The noise-window test now checks `quiet_window` for every mode and decimation factor 1-16 (ends before the diode, whole decimated blocks only) and fails if the window goes back to `noise0`; a second test checks the mode table against `radar_timing.TMM`.

## Benchmark

`c423963`: exit 0, 79 of 79 files bit-identical; report `~/isr_project/regression/reports/c423963-c423963_vs_8206c84-8206c84.md`. This says nothing about the change: the benchmark does not run the catalogue. The change is product-changing (every detection's `noise_power`, `matched_filter_power`, `snr_db` and `cfar_ratio`; `doppler_hz` slightly; mode 300 too, its window moving from 7300-7800 to 7200-7700), so it needs gate B.

Together with the other six on branch `integration-all-fixes` (3fbc1ef, report `~/isr_project/regression/reports/3fbc1ef-3fbc1ef_vs_8206c84-8206c84.md`): exit 1, 0 of 81 files bit-identical, as expected. Every combined change is the sum of the single branches' changes, except one interaction: misa-l period 10875 (`pp-1712592750`), which contains an antenna switch. `antenna-switch-gap` alone loses that fit (108 of 108 gates); with `lpi-unconstrained-gates` it is recovered, and n_e changes by a median of 35 % (fewer pulses, the wrong-antenna ones removed). That is why those two, `fit-lpi-range-avg` and `fit-lpi-last-group` must be merged together.

## Scientific review (gate B, and memos with results)

Reviewer: `isr-science-reviewer` (Opus 5.5, a fresh subagent), 2 October 2026, together with `tx-delay-own-antenna`. Report: `documents/documents_logs/reviews/2026-10-02-gate-b/txdelay-and-catalogue.md`, section 2; scripts and outputs in `documents/documents_logs/reviews/2026-10-02-gate-b-txdelay-catalogue/`. At the caller's request it also ran the detector on a sample of raw pulses (8 processes, under `systemd-run`).

- **Verdict: sound with changes.** The fix is correct and Memo 28's evidence supports it. "Detections unchanged" has to become "about 2 % of detections, all marginal, flip", and the stored catalogue should be rerun.
- **Full verification: not needed.** The catalogue's SNR enters no T_sys, density scale or pipeline product, and the reviewer's run on 360 pulses already recomputes the key claims independently.
- **Spot checks:** `quiet_window_profile.py`, `mode800_outlier_pulses.py` (where the diode switches on, and the noise factors); `cfar_window_flip_check.py`, `cfar_window_flip_summary.py` (detections, main against branch, 360 pulses).
- **Findings and what was done:**
  - **C1 (confirmed), where the diode switches on and whether the windows are clean** (60 own-antenna pulses per case). Switch-on in the mean power profile: 8280-8300 µs for the coded modes (6.5-8 times the background), about 7800 µs for mode 300, about 30100 µs for mode 800. The new windows are flat: 0.94-1.11 in 20 µs bins on 7 April at 18:31. The catalogue's estimator, old window over new, per pulse (median): coded 1.36 at 16:05 UTC and 1.32 at 18:31 (Memo 28: 1.324); mode 300 1.00 and 1.02 (Memo 28: 1.000); mode 800 1.14, mean 1.175 (Memo 28: 1.174). Caveat: mean-over-pulses profiles are pulled up by single interference pulses (in one mode-800 pulse of eight the background was 50 times normal).
  - **E1 (error), "which echoes the catalogue detects does not depend on this" (Memo 28 §2; Memo 35) is not exact.** `detect_pulse` from main (a7e3cb5) and from the branch (387d933) on 360 pulses (16:05 UTC, the bright pass of Memo 16, 18:13 UTC, and a mode-800 block); main reproduces the stored catalogue in 360 of 360. Detections go from 105 to 103: 2 lost, 0 gained, at 1.042 and 1.019 times the threshold (one coded, in the bright pass; one mode 800). Doppler changes reach 91 Hz (mode 300) and 42 Hz (coded), not "within 20 Hz" (one bin is about 460 Hz). Per-echo SNR change: coded median +1.2 dB (5-95 %: −0.3 to +2.5), mode 800 +0.9 dB (−0.5 to +2.3). *Done:* corrected here and in Memo 35 (Memo 28 still says "unchanged").
  - **G2 (gap), "no pipeline step uses snr_db" holds for the products, but analyses use it:** Memo 37's cuts at 15, 20 and 25 dB, `validate_satellite_real.py:66-67` (passes at 20 and 30 dB), and the energy ratios of Memos 16 and 28. These subsets are tied to the biased SNRs, and the coded and mode-300 detections are biased differently (about +1.2 dB against 0). *To do:* mark Memo 37's SNR cuts as made on the old SNRs.
  - **F3, rerun the stored catalogue.** The 2 % flips and the per-pulse scatter of the noise factor (16-84 %: 1.13-1.66; Memo 28: ±15 %) cannot be corrected after the fact. A correction would have to be per mode in linear units, SNR_new = (SNR_old + 1)·q − 1 with q = 1.324 (coded) or 1.174 (mode 800); a flat 1.3 dB errs by 0.1 dB at 10 dB and 0.01 dB at 20 dB. *Done:* Henrik decided to rerun the catalogue after the merge.
  - **S4 (suggestions, pre-existing):** the 62 decimated samples behind each noise estimate give a per-pulse scatter of about 1/(√62 · ln 2) ≈ 18 %, as observed. At the rerun, also consider the catalogue's fixed `--receiver-delay-samples 11`: the measured own-antenna delays are 10.6-10.8 µs (zenith-l) and 12.3-12.8 µs (misa-l), so `range_km` is off by −60 to +270 m. That would be a separate change.
- **Claims the reviewer could not confirm:** Memo 28's corrected energy ratios of 1.08-1.15 (not recomputed); its 3100-pulse medians (reproduced with 60-pulse samples only, and mode 800 from one 10-minute block only).
- **Questions for the supervisor:** is a rerun of the catalogue (hours) acceptable, and should it use per-channel receiver delays?

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

Not needed (scientific review): the catalogue's SNR enters no pipeline product, and the reviewer's detection run on 360 pulses already recomputed the key claims independently.

## Code review

- **No record:** this file.
- **The test did not test the change:** fixed in `c423963`.
- **"Detections unchanged" holds in expectation, not exactly:** the DC offset now comes from a different window, so a cell at the CFAR threshold can flip. Gate B can compare detection counts on a sample, including mode 800. Done at gate B: about 2 % of the detections flip (2 of 105 lost, both marginal; finding E1 below).
- **Stored catalogue and scripts drift** (`validate_satellite_real.py`'s energy ratio) until the catalogue is rerun: decision open (rerun, or correct the stored SNRs), TODO.tex. Decided at gate B: rerun after the merge (Henrik).
- **Pre-existing:** the CFAR search still ends at `noise0`, so the diode can sit in the training cells of the last valid rows; not a regression, out of scope.

## Decision

- Gate A: passed on 2026-10-02 (code review answered, tests pass, every benchmark change explained or handed to gate B). Not merged: the change is product-changing.
- Gate B: scientific review done on 2026-10-02 (above); approved by Henrik on 2026-10-02. No code change. After the merge: rerun the stored satellite catalogue, and mark Memo 37's SNR cuts as made on the old SNRs.
