# Review record: `fit-lpi-range-avg`

- **Gates:** A and B (product-changing)
- **Author:** Claude (Opus 5.5) sessions for Henrik: the fix on 2026-09-30, the gate A review fixes on 2026-10-01
- **Reviewer:** `isr-code-reviewer` (Sonnet 5.5, a fresh subagent), 2026-10-01; report in `documents/documents_logs/reviews/2026-10-01-gate-a-code-reviews/fit-lpi-range-avg.md`
- **Commits:** `446d239`, and at gate B `c05e758` (with merges of main: 8206c84 in `f52508e`, a7e3cb5 in `90957ac`); based on `main` at `a7e3cb5`.

## What the change does

In `fit_lpi`'s range average the variance slice lacked `max(0, ...)`, so the lowest `ra` gates got infinite variance, and both slices averaged 2 ra gates, half a gate off centre, instead of the stated 2 ra + 1. Both now use `[max(0, ri-ra), min(n, ri+ra+1))`.

At gate B (`c05e758`) the average itself was corrected, as `range_average()`: at each lag only the gates with finite data count in the r² weights (the weights of NaN gates pulled the average towards zero); the variance is Σw²v/(Σw)², the variance of the r²-weighted mean, instead of 1/Σ(1/v); and where the centre gate has no data at a lag (empty, or masked as space debris) the result is NaN. This is per lag: the lowest gates with data (8-14 here) hold only their short lags, and those lags are no longer filled from the gates above, which would put them off centre. Gate 8 (72 km), with 6 of 45 lags, is then not fitted, as without averaging.

## Tests

`python3 -m pytest -q` (with the ion-line tables linked into `data/`): 28 passed after `c05e758`. New in `c05e758`: `test_range_average.py` (the weights and the variance against a direct sum, the window edges, no bias next to a gate without data, NaN at an empty, a masked and a partly empty centre gate).

## Benchmark

`c05e758` (the gate B correction): not run yet; it needs a new benchmark before the merge. The numbers below are for `f52508e`.

`f52508e`: exit 1; report `~/isr_project/regression/reports/fit-lpi-range-avg-f52508e_vs_8206c84-8206c84.md`. The 8 `fit_lpi` files change and the other 71 are bit-identical. Each fit gains its lowest gate (finite where it was infinite variance); *corrected at gate B:* that gain was spurious, the 63 km gate has no data and its value was gate 8's scaled by about 0.43 (finding E1 below; `c05e758` makes it NaN). Parameters move with the recentred window (zenith-l n_e by a median relative 0.16-0.20 in these 60-s fits, error bars 1.3-7.7 % smaller). The empty fit of misa-l period 2 changes v_i from +100 to -100 m/s: explained, see below. The commit message's numbers (85 km gate, 0.2-0.7 σ) come from an earlier test against ba37540, not from the benchmark.

Together with the other six on branch `integration-all-fixes` (3fbc1ef, report `~/isr_project/regression/reports/3fbc1ef-3fbc1ef_vs_8206c84-8206c84.md`): exit 1, 0 of 81 files bit-identical, as expected. Every combined change is the sum of the single branches' changes, except one interaction: misa-l period 10875 (`pp-1712592750`), which contains an antenna switch. `antenna-switch-gap` alone loses that fit (108 of 108 gates); with `lpi-unconstrained-gates` it is recovered, and n_e changes by a median of 35 % (fewer pulses, the wrong-antenna ones removed). That is why those two, `fit-lpi-range-avg` and `fit-lpi-last-group` must be merged together.

## Scientific review (gate B, and memos with results)

Reviewer: `isr-science-reviewer` (Opus 5.5, a fresh subagent), 2 October 2026. It judged the four linked branches, `lpi-unconstrained-gates`, `antenna-switch-gap`, `fit-lpi-range-avg` and `fit-lpi-last-group`, as one change. Report: `documents/documents_logs/reviews/2026-10-02-gate-b/four-linked-branches.md`; scripts in `scripts_four/` beside it.

- **Verdict: sound with changes.** Each of the four fixes does what it claims: the off-centre range window biased n_e high by about 4-5 %, the 2.5 s guard removed every pulse sent on the other antenna, and the recovered lags leave every value main already had unchanged. The required changes are in `fit-lpi-range-avg` (E1), `fit-lpi-last-group` (E2) and `antenna-switch-gap` (G1).
- **Spot checks** (report section 3):
  - Window bias (`win.py`), on a smooth cubic-in-log model of the 300-s zero-lag profile: the old window +4.4 to +5.2 % at 300-650 km on zenith-l and +4.6 to +4.8 % on misa-l, the new one −0.3 to +1.4 %; new/old on the data, zenith-l: 0.98 (200-300 km), 0.94 (300-400 km), 0.95 (400-600 km). **Agrees** with Memo 35.
  - Error-bar shrink: with Memo 36's correlation of neighbouring gates, ρ₁ = −0.30, going from 2 to 3 gates truly shrinks σ by 0.72, against 0.82 reported, so the smaller error bars are not overstated.
  - The NaN-weight bias at the 63 km gate and at debris gates (`low.py`, `deb.py`): E1 below. The variance against the r²-weighted mean (`var.py`): S1 below.
  - The 300-s integration profiles are plausible for midday at solar maximum: hmF2 about 280 km; Te 1800-2500 K and Ti 1000-1200 K at 200-450 km; topside scale height about 200 km, against k(Te+Ti)/(m_O g) ≈ 210 km; v_i small.
- **Findings that concern this branch:**
  - **E1 (error), gates with no data entered the range average's weights.** The average is Σr²·acf/Σr², and the denominator included the r² of gates whose ACF is NaN, pulling the average towards zero. *Evidence:* the "recovered lowest gate" was spurious. The gate at 63 km (index 7) has no finite LPI value at any lag (`lpi-1712593940.h5`: gates 0-7 all NaN); its fitted value was gate 8's data reduced to about 0.43 (72²/(54²+63²+72²)), fitted as Te 9000/Ti 3000 or Te 450/Ti 150 K, the limits. The `max(0, ...)` itself changes nothing in these data, since gates 0-7 are always empty. And at debris-masked gates in single-file fits the fit wrote biased-low n_e: `zenith-l/pp-1712588710`, masked 549-576 km, 94, 136 and 41 against about 300 at the neighbouring gates; `pp-1712514700`, masked 306-333 km, 778, 868, 599 and 491 against about 1900 at 288 km. *Done in `c05e758`:* only finite gates count in the weights, and the average is NaN where the centre gate has no data or is masked (Henrik).
  - **S1 (pre-existing), the variance 1/Σ(1/v) does not match the r²-weighted mean** at the lowest gates: in the 300-s stretch at 60-150 km the error bar was up to 36 times too small (zenith-l) and 1.9 times (misa-l); at 150 km and above the two agree to 1 %. This explains misa-l's 0.35 error-bar ratio. *Done in `c05e758`:* the variance is Σw²v/(Σw)² (Henrik).
  - **G4 (gap), whether the new densities are closer to the truth.** There is evidence for the profile shape and the time series (the range-window shift is a removable bias, check above), not for the absolute n_e: Memo 35's "up to 8 %" includes `mode300-injection-window`, and the 2.2 times aperture-efficiency gap is untouched. The decisive comparison is Memo 38's C with the merged code (full verification below). Smaller point: Memo 35's account of the −4 % leaves out a term: the old window's r² normalisation adds about g/r (about 3 % near the peak, g = the 9 km gate spacing) on top of the gradient term g/2H, with H, the topside scale height, about 150-200 km here, not 100. *Done:* Memo 35 corrected.
- **The other findings of the joint review** (E2 incremental runs, G1 guard, G2 period 10875, G3 fit limits) are in those branches' records.
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

- **No record:** this file.
- **`fit_ionline.py` has the same bug** (a verbatim copy, lines 710-722): it looks unused; TODO.tex.
- **The ±100 m/s:** a pre-existing bug in `fit_acf`. With lag 0 lost the fit returns its start guess, and `fit_acf` flips the caller's guess in place, so a gate that raises leaves the next flat gate at -100 m/s. Not caused by this branch; TODO.tex (fix: copy the guess, NaN without lag 0). `lpi-unconstrained-gates` removes the case here.
- **The variance of the mean** (1/Σ(1/v)) does not match the r²-weighted ACF mean, and a NaN gate still counts in the ACF mean's weights: pre-existing; decision open (TODO.tex). Fixed at gate B in `c05e758` (Henrik's decision).
- **No unit test** of the slice: it would need the averaging pulled into its own function; left for now. Done at gate B: `range_average()` and `test_range_average.py` in `c05e758`.
- **Edge windows** hold fewer than 2 ra + 1 gates while `range_avg_window_km` stays nominal; the default `range_avg` does not average below 300 km.
- **Merge order:** `fit-lpi-correlated-gates` has the same slice edit; no conflict.

## Decision

- Gate A: passed on 2026-10-02 (code review answered, tests pass, every benchmark change explained or handed to gate B). Not merged: the change is product-changing.
- Gate B: scientific review done on 2026-10-02 (above); approved by Henrik on 2026-10-02 with the corrected average. `c05e758` needs gate A (code review, benchmark) before the merge, together with `lpi-unconstrained-gates`, `antenna-switch-gap` and `fit-lpi-last-group`. The scoped full verification follows the merge.
