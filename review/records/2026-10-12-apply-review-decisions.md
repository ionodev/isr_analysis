# Review record: `apply-review-decisions`

- **Gates:** A and B (product-changing). Gate B is Henrik's own review (project
  policy of 3 October; pull request #7).
- **Author:** Claude session of 2026-10-12, for Henrik
- **Reviewer:** `isr-code-reviewer` subagent on Sonnet 5.5 (no part in the work)
- **Commits:** `da1bdac..0f7292f` with the merge of `main` (cd270de); built on
  pull request #4's branch (a39c642), `main` at 3659f19 after #4's merge

## What the change does

Applies Henrik's review of every change on `main` against Juha's repository
(`documents/overview/repo_comparison.tex`, categories 1-3, 12 October). What
he approved stays as it is. The branch applies the rest:
- the space-object cut of `fit_lpi.py` (and `fit_ionline.py`) is off by default
  (`space_object_cut=None`); `"absolute"` is Juha's test (width < 300 and its
  uncertainty < 100 m/s), `"significance"` the previous default (width + 3
  sigma < 300 m/s); `run_analysis.py` passes it from the `fit_lpi` section;
  the output records the cut used. Why: satellites are handled by the CFAR
  catalogue and the satellite columns, and the cut would otherwise remove
  satellite gates after every method under comparison;
- removed: the zenith beam constants (unused, superseded by Memo 37), the
  catalogue delay refinement (no effect on the plasma, Memo 29) and impulse
  blanking (to be replaced by impulses as unknowns in the inversion);
- category 4 (analysis and validation scripts) leaves `main`, kept on the
  branch `main-before-review-2026-10-03`; the unit tests stay.

## Tests

`python3 -m pytest -q`: 50 passed, none skipped (theory tables linked). New:
`test_space_object_cut.py` checks the default (off, in both fit modules),
both cuts at their boundaries, and that an unknown cut is refused.
`test_satellite_columns.py` is removed: it tested only the delay refinement.

## Benchmark

`python3 ~/isr_project/isr_analysis/review/regression.py apply-review-decisions
--base origin/fix-branches` (pull request #4's code, now `main`'s): exit status
2. Report: `~/isr_project/regression/reports/apply-review-decisions-da1bdac_vs_origin_fix-branches-a39c642.md`.

**Why 2, and why the run still answers the question.** The only reason is
`lpi_60/zenith-l/lpi-1712594080.h5` (period 11008, "made only of
wrong-antenna pulses"), missing on *both* sides: since pull request #4
(`antenna-switch-gap`) the pipeline rightly makes no output for it, and #4's
own benchmark (f8c96fe) shows it as an intended difference. The tool still
expects it, so every benchmark against `main` is now invalid for this reason
alone. 0 failed jobs; both sides produced the same 80 files. The tool's
expected list needs updating (a change to `review/`, which needs Henrik's
approval; not done here).

**What changes.** 70 of 80 files are bit-identical: every LPI file, the
range-Doppler spectra and the long-pulse fits. The 10 `fit_lpi` files differ:
- all 10 gain the dataset `space_object_cut` ("none");
- the 4 `misa-l` fits and one `zenith-l` fit had no flagged gate: nothing
  else changes;
- the other 5 `zenith-l` fits had 14 flagged gates in the base, 13 at
  495-567 km (the Starlink shells, Memo 29) and one at 324 km. Each flag
  blanked 4 gates in that file. Now nothing is blanked: 4 gates are recovered
  (NaN to finite) per fit where a whole group was blanked, and where only
  some files of a group were, those files now enter the average, so values
  at those gates change (relative changes up to 0.9 in Te, and large ones in
  vi and Sigma near zero values). No value is lost.

Every change is at a gate the cut flagged, as expected. **Consequence:** in a
run without the satellite columns, these gates now carry the satellite echo
into the fit instead of being blanked (the open TODO decision "Mask or flag
the zenith gates at 535-590 km?"). With the columns, the comparison methods
handle them.

## Scientific review (gate B, and memos with results)

Henrik's own review (policy of 3 October): the decisions are his; the
evidence for each is in `repo_comparison.tex` and the discussion of 12
October. No AI scientific review requested.

## Full verification (when flagged, and always for calibration)

Not needed: no calibration changes, and the effect is the removal of an
option's default, confined to the flagged gates.

## Code review

`isr-code-reviewer` on da1bdac: no blocking or major findings. Minor:
1. README still named removed files: read while being edited; the committed
   README (0f7292f) names none.
2. `.gitignore` comment named `plot_raw_voltage.py`: fixed in 0f7292f.
3. `fit_ionline` refusal of an unknown cut untested: test added in 0f7292f.
   No test runs `fit_lpifiles` end to end with the cut off: that needs the
   raw data; the benchmark covers it (no gate blanked).
4. `config/millstone_2023-09-05.json` has no `space_object_cut` key: the
   default applies; documented in README and `run_analysis.py`.
5. The bare `except:` around `fit_gaussian` in `fit_lpi.py` (existing, Juha's)
   hides a failed fit when the cut is on: unchanged, only matters when on.
The test commit 0f7292f changes `test_*.py` and documentation only, so the
benchmark of da1bdac stands.

## Changelog

Under "Waiting for gate B" in `CHANGELOG.md`; moved to the merge date before
the merge.

## Decision

- Gate A: tests and code review passed; benchmark valid but for the tool's
  expected-output list (above).
- Gate B: waiting for Henrik's approval of the pull request.
