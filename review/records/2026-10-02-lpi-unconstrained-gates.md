# Review record: `lpi-unconstrained-gates`

- **Gates:** A and B (product-changing)
- **Author:** Claude (Opus 5.5) sessions for Henrik: the fix on 2026-09-30, the gate A review fixes on 2026-10-01
- **Reviewer:** `isr-code-reviewer` (Sonnet 5.5, a fresh subagent), 2026-10-01; report in `documents/documents_logs/reviews/2026-10-01-gate-a-code-reviews/lpi-unconstrained-gates.md`
- **Commits:** `a999223`, `33022a1` (with merges of main: 8206c84 in `0877912`, a7e3cb5 in `5c59e97`); based on `main` at `a7e3cb5`. No code change at gate B.

## What the change does

With the default outlier rejection, the rejection can empty the clutter gates of a lag; that range gate is then an unknown no measurement constrains, the normal matrix is singular, and the whole lag was lost (36 of 78 MISA periods in a test hour, a median of 11 of 45 lags; Memo 29). `invert_normal` inverts for the constrained unknowns and leaves the unconstrained ones NaN; with every unknown constrained it is the plain inverse, as before.

## Tests

`python3 -m pytest -q` (with the ion-line tables linked into `data/`): 26 passed. New: `test_outlier_lpi.py` (the solve against least squares on the constrained columns; NaN exactly at the zero columns; bit-identical to `inv` when all are constrained; a NaN in the matrix is not taken for unconstrained).

## Benchmark

`33022a1`: exit 1, 74 of 79 files bit-identical; report `~/isr_project/regression/reports/33022a1-33022a1_vs_8206c84-8206c84.md`. The same differences as `0877912` before the review fixes (`lpi-unconstrained-gates-0877912_vs_main-8206c84.md`), so the refactor is output-neutral. All 5 changed files are misa-l: three LPI files get back every lost lag (8, 11 and 12 lags; every value main had is unchanged, only NaN became finite), and two fits change because of that. Period 2's fit (`pp-1712484020`) had no n_e on main; on the branch it fits at 107 of 108 gates, but many gates sit at a fit limit: finite is not the same as usable (for gate B). The counts given here at gate A (4 at the Te bound of 9000 K, 7 at Ti = 3000 K, 9 at |v_i| = 1500 m/s) missed the Te/Ti and Ti = 150 K hits; gate B counts 45 % of the gates at a limit, see below.

Together with the other six on branch `integration-all-fixes` (3fbc1ef, report `~/isr_project/regression/reports/3fbc1ef-3fbc1ef_vs_8206c84-8206c84.md`): exit 1, 0 of 81 files bit-identical, as expected. Every combined change is the sum of the single branches' changes, except one interaction: misa-l period 10875 (`pp-1712592750`), which contains an antenna switch. `antenna-switch-gap` alone loses that fit (108 of 108 gates); with `lpi-unconstrained-gates` it is recovered, and n_e changes by a median of 35 % (fewer pulses, the wrong-antenna ones removed). That is why those two, `fit-lpi-range-avg` and `fit-lpi-last-group` must be merged together.

## Scientific review (gate B, and memos with results)

Reviewer: `isr-science-reviewer` (Opus 5.5, a fresh subagent), 2 October 2026. It judged the four linked branches, `lpi-unconstrained-gates`, `antenna-switch-gap`, `fit-lpi-range-avg` and `fit-lpi-last-group`, as one change. Report: `documents/documents_logs/reviews/2026-10-02-gate-b/four-linked-branches.md`; scripts in `scripts_four/` beside it.

- **Verdict: sound with changes.** Each of the four fixes does what it claims: the off-centre range window biased n_e high by about 4-5 %, the 2.5 s guard removed every pulse sent on the other antenna, and the recovered lags leave every value main already had unchanged. The required changes are in `fit-lpi-range-avg` (E1), `fit-lpi-last-group` (E2) and `antenna-switch-gap` (G1).
- **Spot checks** (report section 3):
  - `invert_normal`, from the code: a zero diagonal of AᴴA means a zero column of the design matrix A, so leaving that unknown NaN is the exact least-squares solution for the rest. In the integration run the NaNs left are at 72-99 km only (`nan.py`). **Agrees.**
  - Fit-limit counts in every base and integration fit (`bounds.py`, `bh.py`): G3 below.
  - Per-file calibrated power and P_tx, and the period 10875 profiles (`pw.py`, `q2.py`, `q2b.py`): G2 below.
- **Findings that concern this branch:**
  - **G3, fits at a limit are common, not special to the recovered first period.** A limit is Te/Ti 0.99 or 3, Ti 150 or 3000 K, or |v_i| 1500 m/s. On main, 41-59 % of the gates of the 10-s fits sit on at least one; the recovered `pp-1712484020` has 45 %: Te/Ti = 3 at 32 gates, Ti = 150 K at 14, Ti = 3000 K at 7, |v_i| = 1500 m/s at 6. The hits are below 150 km and above 500 km height. At 150-500 km the period (10:00 UTC, before ground sunrise, with the F region already sunlit) is noisy but physical: Te/Ti 1.5-2.7, hmF2 about 300 km. *Done:* the fit is kept; the per-gate fit-limit flag (Memo 36's proposal) comes in weeks 7-9 (Henrik). Not blocking.
  - **G2, the "median 35 %" at period 10875** (`pp-1712592750`, the interaction in the Benchmark above). It is the median *absolute* relative change; the signed median ratio is +25 % over all gates and +28 % at 150-500 km height, 1.1σ per gate. Its parts: main took lag 0 only from `lpi-1712592900`, the file with wrong pulses, whose calibrated power rises by 1.31 and P_tx from 1109 to 1249 kW with the fix (its neighbours have 1237-1257 kW); and `mode300-injection-window` (+6.6 % in α) is in the integration too. The direction and size are plausible. *Open:* after cleaning, 592900 still has only 0.80 of 592750's power (residual contamination, MISA moving, or a real change).
  - **S1 (pre-existing, not changed here):** P_tx is a plain mean over files while the ACF is variance-weighted, so dropping a 448 kW file moved n_e by −2.7 % (the fix: normalise each file by its own P_tx). Each gate's start guess is the result of the gate below, so changes at the bottom gates spread upwards: the origin of the ±100 m/s case.
- **The other findings of the joint review** (E1 range average, E2 incremental runs, G1 guard, G4 absolute scale) are in those branches' records; E1, E2 and G1 are fixed there.
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

- **No test** of the solve: fixed in `33022a1` (helper `invert_normal`, `test_outlier_lpi.py`).
- **No review record:** this file.
- **A NaN diagonal counted as unconstrained** (`NaN > 0` is false), which would solve around NaN-contaminated rows silently: fixed in `33022a1`, only an exactly zero diagonal counts.
- **The benchmark summary's "no finite value at all" for period 2** was wrong (main held the fit's start values at one gate): corrected in the summary.
- **Log volume** (one line per affected lag and period): noted, left.

## Decision

- Gate A: passed on 2026-10-02 (code review answered, tests pass, every benchmark change explained or handed to gate B). Not merged: the change is product-changing.
- Gate B: scientific review done on 2026-10-02 (above); approved by Henrik on 2026-10-02. No code change here; to be merged together with `antenna-switch-gap`, `fit-lpi-range-avg` and `fit-lpi-last-group` once their gate B changes have passed gate A. The scoped full verification follows the merge.
