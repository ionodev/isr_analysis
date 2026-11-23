# Review record: `clutfix-impulse-rows`

- **Gates:** A and B (product-changing for product-v3; output-neutral at the defaults)
- **Author:** Claude session (Opus 5.5), 22 November 2026, overnight, at Henrik's request
- **Reviewer:** `isr-code-reviewer` agent (Sonnet 5.5), a fresh subagent
- **Commits:** `7575c82` (the change), based on `weeks-10-12` at `88b9ef7`

## What the change does

With `impulse_rows=True`, `outlier_lpi` removes every lagged-product row
touched by a detected power-line impulse. It removed each impulse's rows
from the pulse and from its clutter partner (the other pulse of the same
code) for both kinds of lagged product, with one row set, one weight vector
and one solve. The partner's impulses are only in the clutter-subtracted
(differenced) products, `acfs_g`, which `fit_lpi` uses below
`clutter_cutoff_km` (150 km in product-v3); the plain products, `acfs_e`,
used above it, contain the pulse's own samples only. Rows hit only by the
partner's impulses were therefore removed from `acfs_e` without need: for
independent impulses about twice the rows needed. Now the plain products
lose only the rows of the pulse's own impulses
(`impulse_rows.touched_masks`); where the two masks differ at a lag, each
kind of product gets its own leave-one-out weights (the running mean over
its own kept rows, the row cap and the median floor as before) and its own
solve. The differenced products keep exactly their old rows, weights and
solve. `acfs_var` is then the plain products' variance, and the new
`acfs_var_g` the differenced one's, which `fit_lpi` takes for the gates it
takes from `acfs_g`. New key `impulse_rows_removed_fraction_e`; the existing
`impulse_rows_removed_fraction` keeps its meaning (the union, removed from
the differenced products).

## Tests

`python3 -m pytest -q`: 105 passed, none skipped (tables linked in). New tests
in `test_impulse_rows.py`: `test_touched_masks_without_partner_is_one_mask`;
`test_partner_impulses_touch_only_the_differenced_products` (an impulse in the
partner leaves every plain row of the pulse unchanged, and is present in the
differenced rows only it hits and absent from the rows the union mask keeps).

## Benchmark

`python3 ~/isr_project/isr_analysis/review/regression.py clutfix-impulse-rows --base weeks-10-12 --jobs 8`
(main's copy of the tool; the branch contains `main`): exit status 0, "Every
output is bit-identical", 80 output files compared, run 22 November 23:40 UTC.
Report: `~/isr_project/regression/reports/clutfix-impulse-rows-7575c82_vs_weeks-10-12-88b9ef7.md`.

Why no change, although the change alters products: the benchmark runs
`outlier_lpi` at the pipeline defaults, where `impulse_rows=False`, so the new
masks, weights and solves are never reached; with the option off the code
takes the one joint solve, unchanged. The change takes effect only in runs
with `impulse_rows=True` (product-v3 and the rfi_rows variants), which the
benchmark does not cover (REVIEW_PROCESS.md section 1: options off by
default). It was shown in two other ways: (1) a check run of 8 April
periods 10471-10473 (the benchmark's power-line period and the next two)
with the product-v3 configuration, `weeks-10-12` against this branch
(`/mnt/data/henrik/clutfix_2026-11/check/{old,new}`): `acfs_g`, `noise_g`,
`sat_amp_g` bit-identical, the new `acfs_var_g` equal to the old `acfs_var`,
rows removed from the plain products 29, 27, 23 % -> 16, 15, 13 %, the plain
ACF's error bar at lag 10 us 0.94-0.96 of the old (median 150-600 km), mean
change of the plain ACF +0.016, -0.026, -0.013 sigma (lags 10-60 us); (2) the
full rerun of both days and of the injection test (memo 55), where the same
identities are checked on every file.

## Summary for Henrik's review (gate B)

What was wrong: above 150 km product-v3 fits the plain ACF, from which the
impulse removal also took the rows hit only by the clutter partner's
impulses, about twice the rows needed; the cost of the removal (error bars,
Table 1's share of rows) was overstated, with no bias. What changes: with
`impulse_rows=True`, `acfs_e`, `acfs_var`, `noise_e`, `sat_amp_e`,
`sat_amp_var` and `sat_resid` at lags where the masks differ; `acfs_g` and
everything below 150 km do not. The size, on both days, the injection test
and the plasma parameters: memo 55. Not checked: other readers of `acfs_var`
(below).

## Scientific review (only if Henrik asked for one)

Not requested.

## Code review

`isr-code-reviewer` (Sonnet 5.5), 22 November: passes, no must-fix; bit
identity with the option off and in the non-split case confirmed by reading;
the split case's rows, weights and variance consistent; the real-data check
verified on the h5 files. Caveats, answered here, none changed in the code:
- With `outlier_rejection=True` and a split lag, `acfs_g` no longer loses the
  rows where only `mm_e` failed the ratio or power cut (the old joint row set
  required both finite), so it would differ from the old output there. No
  product runs rejection with the impulse rows (all configs have
  `outlier_rejection: false`); each product now keeping its own rows is the
  consistent rule.
- In split mode the rejection statistics (`outl_*`, `sat_rows_rejected`) use
  the union weights, and `retained_measurement_fraction` counts rows kept in
  both products, i.e. the union mask: conservative for the plain products,
  whose share removed is in `impulse_rows_removed_fraction_e`.
- The tests check the masks and the physical claim, not the split solve;
  that is covered by the real-data check and the rerun.
- `fit_ionline.py`, `ion_line_fits.py` and `review/regression.py` read
  `acfs_var` for `acfs_g` gates; in split lags that is now the plain
  products' variance, which understates the differenced products' error.
  None is in the product path (`fit_lpi` is updated); the benchmark has no
  split lags.

## Changelog

Entry under "Waiting for gate B" in `CHANGELOG.md`.

## Decision

- Gate A: tests, benchmark, code review and record done on 22-23 November 2026.
- Gate B: waiting for Henrik (whether to adopt the fix in the paper and the product).
