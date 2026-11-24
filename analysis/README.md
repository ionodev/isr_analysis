# Analysis scripts

Scripts that test, evaluate or plot the pipeline's results for the memos.
They are not part of the pipeline: nothing in the repository's top folder
imports them, and the products do not depend on them. They import the
pipeline modules from the folder above.

Run them from the repository's top folder, e.g.
`python3 analysis/compare_weights.py ...`; the usage is in each script's
docstring (written before the move, so it says `python3 compare_weights.py`).
To rerun a memo exactly, check out the commit or product tag its Reproduce
section names (`product-v1` to `product-v3.1` have the scripts in the top
folder).

| Script | What it does | Memo |
|---|---|---|
| `inject_satellites.py` | synthetic satellite pass in a quiet interval, through each LPI variant | 42 |
| `compare_injection.py` | bias of the LPI variants against the injected truth | 42 |
| `compare_variants.py` | comparison measures of the five LPI variants | 42 |
| `compare_fits.py` | the LPI variants compared at the level of the fitted parameters | 44 |
| `interval_measures.py` | precision, honesty and echo residual on the test intervals | 42, 45, 48 |
| `change_table.py` | how the fitted parameters of one run differ from another's, by altitude band | 45, 48 |
| `column_checks.py` | per-column checks: is there a coherent point target in each satellite column | 46, 48 |
| `random_control.py` | the pure-plasma control: satellite columns where there is no echo | 46, 48 |
| `rejection_measures.py` | what the outlier rejection removes next to the satellite columns | 48 |
| `confirmation_count.py` | false alarms after the tracking confirmation | 48 |
| `end_to_end_injection.py` | synthetic passes through detection, tracking and the inversion | 45, 48 |
| `e2e_threshold_figure.py` | figure of the end-to-end test at three detection thresholds | 48 |
| `compare_impulse_rows.py` | the impulse-row removal on minutes with power-line sparks | 47 |
| `inject_impulses.py` | synthetic power-line impulses, with and without the row removal | 51 |
| `weight_bias_mc.py` | Monte Carlo of the bias of data-estimated variance weights | 51, 54 |
| `compare_weights.py` | variance weights against leave-one-out weights in the product | 53 |
| `eclipse_timing.py` | the eclipse's start, maximum, end and obscuration at Millstone Hill | 52 |
| `eclipse_analysis.py` | 8 against 7 April: NmF2, hmF2, Te, Ti | 51, 52 |
| `eclipse_plots.py` | altitude-time plots with uncertainties and quality flags | 44, 45 |
| `eclipse_difference.py` | altitude-time plot of one run minus another | 45 |
| `ion_line_fits.py` | measured and fitted ACF and spectrum at a few altitudes | 45, 55 |
| `range_lag_doppler.py` | range-lag and range-Doppler images of one fit window | 45, 46 |

The plasma-line calibration (`plasma_line_picker.py`,
`plasma_line_density_scale.py`) stays in the top folder: it produces the
density constants in `fit_quality.py`, which every product uses.
