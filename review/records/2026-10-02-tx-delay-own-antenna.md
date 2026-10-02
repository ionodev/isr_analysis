# Review record: `tx-delay-own-antenna`

- **Gates:** A and B (product-changing)
- **Author:** Claude (Opus 5.5) sessions for Henrik: the fix on 2026-09-30, the gate A review fixes on 2026-10-01
- **Reviewer:** `isr-code-reviewer` (Sonnet 5.5, a fresh subagent), 2026-10-01; report in `documents/documents_logs/reviews/2026-10-01-gate-a-code-reviews/tx-delay-own-antenna.md`
- **Commits:** `784f00b`, `c908db7` (with merges of main: 8206c84 in `e43b5fd`, a7e3cb5 in `94c1886`); based on `main` at `a7e3cb5`. No code change at gate B.

## What the change does

The channel-delay estimator looked 20 s ahead for coded pulses on the receiver's own antenna and otherwise used the other antenna's, whose leakage arrives 1.2 µs (zenith-l) or 2.3 µs (misa-l) away (Memo 29). The window now doubles, up to `max_search_s` = 2400 s, until it holds `n_pulses` own-antenna pulses; only then does it fall back.

## Tests

`python3 -m pytest -q` (with the ion-line tables linked into `data/`): 27 passed. New: `test_tx_delay.py` with fake readers (own pulses in the first window; a window that has to grow, which fails on main; the cap and the fallback; a recording shorter than the window).

## Benchmark

`c908db7`: exit 1, 64 of 79 files bit-identical; report `~/isr_project/regression/reports/c908db7-c908db7_vs_8206c84-8206c84.md`; the same differences as `e43b5fd` before the review fixes. zenith-l's delay moves from 11.88 to 10.62 µs (the recording starts in MISA's cycle, so main fell back to MISA's pulses). misa-l is bit-identical (its first 20 s hold own pulses).

*Corrected at gate B* (findings E1, E2, F3, F5 below). Gate A said here that the delay shifts every zenith-l range-Doppler spectrum by 1.3 samples (195 m) and moves the debris mask by 2 gates, and left two points unexplained. In fact:
- the effect on the spectra is a gain of about +0.04 % and a shift of about 14 m, not 195 m: the echo window is fixed, so the delay moves only the receive weighting inside it;
- no debris masking happens in the benchmark (`fit_lp` runs with `remove_space_objects` False); only the debris classifier's list changes;
- the single gates that change by 12-23 σ are fits at a parameter limit or ill-determined, not a range shift;
- the 0.1 µs between the pipeline's 10.62 µs and the independent 10.72 µs is a real drift of the delay over the day: the two take their pulses at different times.

Together with the other six on branch `integration-all-fixes` (3fbc1ef, report `~/isr_project/regression/reports/3fbc1ef-3fbc1ef_vs_8206c84-8206c84.md`): exit 1, 0 of 81 files bit-identical, as expected. Every combined change is the sum of the single branches' changes, except one interaction: misa-l period 10875 (`pp-1712592750`), which contains an antenna switch. `antenna-switch-gap` alone loses that fit (108 of 108 gates); with `lpi-unconstrained-gates` it is recovered, and n_e changes by a median of 35 % (fewer pulses, the wrong-antenna ones removed). That is why those two, `fit-lpi-range-avg` and `fit-lpi-last-group` must be merged together.

## Scientific review (gate B, and memos with results)

Reviewer: `isr-science-reviewer` (Opus 5.5, a fresh subagent), 2 October 2026, together with `catalogue-noise-window`. Report: `documents/documents_logs/reviews/2026-10-02-gate-b/txdelay-and-catalogue.md`, section 1; scripts and outputs in `documents/documents_logs/reviews/2026-10-02-gate-b-txdelay-catalogue/`.

- **Verdict: sound with changes.** The fix is right: the reference should be the own antenna's leakage, which is also what the lag-profile inversion uses (`outlier_lpi.py:433-436` takes `z_tx` from the echo channel itself). Its effect on the products is negligible. The changes needed are to the text: this record and Memo 35 described the effect wrongly.
- **Full verification: not needed.** No calibration is touched (T_sys bit-identical), and the key numbers check out against the existing outputs.
- **Spot checks:** `trace_fit_gates.py`, `classify_gates.py` (the >3σ gates, main against branch); `shift_test.py`, `shift_test2.py` (is the spectra change a range shift?); `window_model.py` (the delay's effect on the receive window); the existing `verify_combined_txdelay*.json` (the 0.1 µs gap and misa-l).
- **Findings and what was done:**
  - **E1 (error), "shifts every zenith-l RD spectrum by 1.3 samples (195 m)" is wrong for the range-Doppler path.** `range_dop_spec` (`avg_range_doppler_spec.py:98-106`) multiplies a fixed echo window by `hann(645)·conj(z_tx)`; the pulse in the sky is fixed by the transmission, so moving the reference moves only the receive weighting inside the window, not the range of the echo. Model (`window_model.py`, the pulse as a rectangle at 103-581 µs on tx-h): the weighting's centroid moves by −0.09 samples (−14 m), the gain by +0.04 %. Data (`shift_test2.py`, the 13 cached spectra files): best-fit shift −0.034 to +0.024 km, gain +0.016 to +0.023 %. At 280-450 km a real 0.19 km shift would change the profile by a median 0.19 %; the observed change is 0.04 %, and a pure gain accounts for 71 % of it. Memo 35's "off by 0.2-0.5 µs (30-75 m)" is overstated in the same way: for mode 300 that drift is about 2-5 m. *Done:* corrected in the Benchmark section and in Memo 35.
  - **E2 (error), "moves the debris mask by 2 gates" explains no fit change.** The benchmark runs `fit_lp` with `remove_space_objects=False` (`review/regression.py:103-104`), and `space_object_count` is zero everywhere in both runs: the list is recorded, nothing is blanked. What changes is the debris classifier's list (12 → 14 and 13 → 15 entries; new at 530.6, 584.6 and 589.1 km, 346.3 km lost), flipped by the ~0.02 % input change. That sensitivity matters in production, where blanking is on: one flip removes ±16 gates (±72 km) of a 10-s file. Not tested. *Done:* corrected here and in Memo 35; the classifier's instability is the same as in `mode300-injection-window`'s E2 (a fit deciding at its width bound).
  - **F3 (explained), the 12-23 σ gates are ill-conditioned fits, not the range shift.** Of the 73 gates that changed by more than 3σ, or changed with a NaN σ: 52 have a parameter at a `fit_spec` limit (Te/Ti 0.99 or 5, Ti 150 or 4000 K, v_i ±1500 m/s, ion fraction 0 or 1), 18 are ill-determined (dne/ne > 0.3, or a NaN error), and 3 are left: 220.3 km in both fits, where Te rises by 15-16 K against dTe = 3-4 K while the molecular fraction rises by 0.010-0.012 (the Te-composition trade-off; `dfrac` is NaN at every gate, so dTe leaves out the composition uncertainty), and 948.8 km, a low-SNR topside gate with Te = 856 K. In the well-determined gates at 170-500 km the rms change is 0.18-0.65 σ and does not correlate with gradient × shift (|corr| ≤ 0.16): convergence noise of the Nelder-Mead fit. *Done:* explained here; the fit-limit flag (weeks 7-9) will mark such gates.
  - **G4 (gap, pre-existing), the benchmark's 60-s `fit_lp` fits cannot judge the physics.** At 400-1050 km, 12-24 gates per 200 km band sit at a temperature, ratio or velocity limit, and 15-35 have the ion fraction at 0 or 1; at 600-800 km the median Te/Ti is 1.00-1.31, implausible in the sunlit topside at noon near solar maximum; the lowest gates (162-166 km) sit at Te ≈ Ti, near the 0.99 limit. Separately, the zenith-l `pp` files store `el` = 44.9°, `az` = −140°, MISA's pointing (`fit_lp.py:335`; metadata only). *Not changed here.*
  - **F5 (resolved), the 0.1 µs gap (10.62 against 10.72 µs) is a difference in sampling time, not in the estimator.** On the *same* first 100 own-antenna pulses (60.6-62.0 s into the recording) the independent estimator gives 10.609 µs against the pipeline's 10.620, less than the 1/64-sample grid apart. The 10.719 µs is the median of 120 pulses spread over 40 minutes; its 30-minute blocks are 10.67-10.77 µs, and 16-17 UTC gives 10.83 µs. So the delay really drifts by about 0.1 µs, significantly (the standard error of a 120-pulse median is about 0.005 µs). It does not matter for the products: about 1 m and a gain change of 3e-5 in mode 300.
  - **F6 (consistent), misa-l at 12.34 µs.** Its first 100 own pulses fall at 0-1.3 s; the first 2 h give 12.34-12.39 µs and 16-17 UTC 12.69-12.83 µs, a drift of 0.35-0.49 µs, matching Memo 35 (0.2-0.5 µs) and Memo 29 (12.66-12.85 µs). zenith-l drifts by 0.21 µs.
  - **S7 (suggestions, not done):** take the 100 pulses spread over the search window, not the first 100 in a row (1.4 s); pass the analysis start time to the estimator, as Memo 35 proposes. Until `antenna-switch-gap` is merged, a run that starts just after a real antenna switch still takes wrong-antenna pulses (the two exceptions in Memo 29).
  - **S8 (pre-existing, the reviewer's reading of the code, not checked on data):** `hann(645)` covers samples 0-644, but the pulse sits at about 114-592 in the echo channel's time base, so the weighted centroid is about 29 samples before the pulse centre, and each range-Doppler gate would be centred about 4.3 km below its `rgs_km` label. To settle: compare the range-Doppler and LPI ranges of the F peak or of a hard target. *Not done.*
- **Claims the reviewer could not confirm:** the estimator on real data at the 2400 s cap, at the fallback, in mode 800 and on other recordings (unit tests only); Memo 29's "118 of 120" (not rerun); whether the leakage is the true range zero.
- **Questions for the supervisor:**
  1. Is the leakage through the T/R switch the right range zero for both antennas? How long is the waveguide/feed path beyond it? Is a cross-antenna leakage offset of 1.2 / 2.3 µs what the hardware should give?
  2. Why does the leakage delay drift by 0.2-0.5 µs over a day (MISA also by 0.14 µs within an hour)? Should the range-Doppler path measure the delay per run start, or per period?
  3. Is centring `hann(645)` on samples 0-644 rather than on the pulse intended (S8)?
  4. `fit_lp` reports no composition error (`dfrac` is NaN), and many gates sit at limits in 60-s fits: should σ-based benchmark comparisons of range-Doppler changes be limited to well-determined gates?

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

Not needed (scientific review): the change does not touch calibration (T_sys is bit-identical in the benchmark), and the key numbers check out against the existing outputs.

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
- Gate B: scientific review done on 2026-10-02 (above); the record and Memo 35 corrected; approved by Henrik on 2026-10-02. No code change; it belongs with `antenna-switch-gap` in the merge order.
