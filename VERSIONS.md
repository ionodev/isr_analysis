# Product versions

The pipeline is a test bed that improves in steps; each step that made a
product used in the memos gets a version number (Henrik, 1 November 2026).
The numbers follow the memos' "satellite echo removal version 1 / 2".
Memo 42's LPI variants 1-5 are method variants for the comparison, not
product versions.

| Version | What it adds | Memos | Commit | Runs (`/mnt/data/henrik/`) | Viewer exports |
|---|---|---|---|---|---|
| **v0** | LPI, variance weights and background column; no satellite or impulse handling (variant 3) | 42, 43 | weeks-7-9 of 12 October (Memo 43 lists c35327e) | `density_scale_2026-10/v3_variance_background` | `isr-nosat` |
| (variant 5) | satellite columns on the detected echoes only (P_FA 1e-10 catalogue) | 42 | | `variants_2026-10/v5_variance_background_satellites` | `isr-detected` |
| **v1** | + pass tracking: columns on every pulse of a detected pass (satellite echo removal version 1) | 45 | 31025ca (tag `product-v1`) | `tracked_2026-10/eclipse_v5_tracked` | `isr-tracked` |
| **v2** | + P_FA 1e-6 catalogue with tracking confirmation (version 2), + impulse-hit rows removed; the final product of 27 October | 47, 48, 51 | c0ff2d6 (tag `product-v2`) | `sat_v2_rfi_2026-10/eclipse_040{7,8}` | `isr-v2rfi` |
| v2-190s | v2 refitted with one window per zenith cycle (about 190 s) | 51, 53 | as v2 | `sat_v2_rfi_2026-10/eclipse_040{7,8}_190s` | `isr-v2rfi-190s` |
| v2-loo | v2 with leave-one-out weights (`weights_exclude_self`) | 51, 53 | 0ee02eb (weeks-7-9) | `sat_v2_rfi_2026-10/eclipse_040{7,8}_loo` | `isr-v2rfi-loo` |
| **v3** | v2's inversions refitted: one fit per zenith cycle, the plain ACF above 150 km, 63-km range windows from 230 km, correlated-gate errors from 300 km, a finer table to Te/Ti 5 (candidate; Henrik decides) | 54 | 407b40b (weeks-10-12) | `v3_2026-11/runs/v3_040{7,8}` | `isr-v3-v1` |
| v3-60s | v3 with 60-s fits, for the eclipse | 54 | as v3 | `v3_2026-11/runs/v3_60s_040{7,8}` | `isr-v3-60s-v1` |
| v3-loo | v3's refit of v2-loo's inversions (leave-one-out weights); 60-s variant `v3-loo-60s` | 54 | as v3 | `v3_2026-11/runs/v3_loo{,_60s}_040{7,8}` | `isr-v3-loo-v1`, `isr-v3-loo-60s-v1` |
| v3-loo100 | leave-one-out weights over 100 pulses (`weights_width` 100), LPI rerun, v3's fit; 60-s variant `v3-loo100-60s` | 54 | b671886 (weeks-10-12) | `v3_2026-11/runs/v3_loo100{,_60s}_040{7,8}` | `isr-v3-loo100-v1`, `isr-v3-loo100-60s-v1` |

The `-v1` / `-v2` at the end of an export's name (`isr-v2rfi-v2`) is the
export's own revision (v2 exports carry the quality flags), not the
pipeline version.

All versions use `fit_quality.DENSITY_CONSTANT` = 7.15e8 (the digisonde,
Memo 43) unless their memo says otherwise; Memo 53 measured the plasma-line
constant of v2 (7.63e8) and v2-loo (6.52e8).

## v3: how it is optimised

Every setting changed for v3 is judged by a target fixed in advance: the
measures of Memo 53 against Madrigal, the baseline (retention, scatter,
residual satellite contamination at the Starlink shells, density against
the digisonde and plasma lines), and the injection tests of Memos 48 and 51
for bias. One change at a time, each with a stated physical or statistical
reason; it is kept if it improves the target and makes nothing worse than
v2 or Madrigal. Settings are chosen on 8 April and confirmed on 7 April.

Done 1-2 November (Memo 54). v3 changes only the fit (`fit_lpi`), so its
runs link to v2's 10-second inversions. Its configuration (the `fit_lpi`
step; top-level `table_dir`):

    "max_dt": 300 (60 for v3-60s), "range_avg": [1, 3, 3, 5],
    "range_limits": [0, 230, 300, 700, 1500], "clutter_cutoff_km": 150,
    "range_avg_rho": [-0.30, -0.08, -0.06], "range_avg_rho_min_km": 300,
    "table_dir": "/mnt/data/henrik/v3_2026-11/tables/wide"

The table (Te/Ti 1-5 in 0.2, Ti 100-4000 K in 100 K, ion fraction in 0.05,
densities 1e11 and 1e12 only, the one the fit reads) is made by
`documents/analysis/v3_2026-11/make_table.py`. Every tried change, kept or
not, is in Memo 54's Table 1. Tried and not kept: `fit_select_best`,
`fit_restarts`, `max_lag_us` (switches left in, off by default).

Open (Henrik): the weights. With the plain ACF above 150 km, v3 with
variance weights carries their SNR-dependent bias at the F peak (its density
scale drifts, impulse removal biases n_e by -4.5 %); v3-loo is unbiased but
loses cells at low SNR. Density constants from the plasma lines: v3
7.95e8, v3-loo 6.18e8, v3-loo100 6.20e8 (the exports carry 7.15e8).
Leave-one-out over 100 pulses instead of 9 (2 November) stays unbiased but
recovers neither precision nor low-SNR retention: the gap is the variance
weights' shrinkage, not the leave-one-out weights' noise (Memo 54).
