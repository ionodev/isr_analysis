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
| **v3** | optimised settings (in progress, weeks 10-12) | | weeks-10-12 | | |

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
