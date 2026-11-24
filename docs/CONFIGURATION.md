# Configuration reference

`run_analysis.py` reads one JSON file. This page lists every setting, its
default and what it does. The [user guide](USER_GUIDE.md) explains how to
run the pipeline; `config/example_basic.json` and `config/example_v3.json`
are complete examples. The column "v3" gives the value of pipeline version 3
where it differs from the default.

```json
{
 "experiment": "my_run",
 "data_dir": "/path/to/recording",
 "output_dir": "/path/to/output",
 "steps": {
  "lpi":        { "enabled": true, ... },
  "fit_lpi":    { "enabled": true, ... },
  "long_pulse": { "enabled": false, ... },
  "fit_lp":     { "enabled": false, ... }
 }
}
```

## Top level

| Key | Default | Meaning |
|---|---|---|
| `data_dir` | required | the recording: the folder holding `rf_data/` and `metadata/` |
| `output_dir` | `data_dir` | where the results are written |
| `experiment` | — | a name for the run (not used by the pipeline) |
| `max_time_s` | all | analyse only the first this many seconds of the recording (`lpi`, `long_pulse`) |
| `radar_freq_hz` | 440.2e6 | the radar frequency, for the theory tables and the fits |
| `table_dir` | `data/` | the folder of the theory tables (v3: the wider table, `data/v3/wide`) |
| `steps` | — | the steps below; a step runs when its `enabled` is `true` |

Keys starting with `_` are ignored, so they can hold notes.

## `steps.lpi`: the lag-profile inversion

Writes `<output_dir>/lpi_<range_gate_us>/<channel>/lpi-<t>.h5`.

| Key | Default | v3 | Meaning |
|---|---|---|---|
| `channel` | required | | `zenith-l` or `misa-l` |
| `range_gate_us` | required | 60 | the range gate in µs (60 µs = 9 km) |
| `pass_band_hz` | required | 100000 | the low-pass filter's half bandwidth (Hz) |
| `filter_len` | 20 | | the low-pass filter's length (samples) |
| `avg_dur` | 10 | | the integration period: seconds of data inverted together |
| `min_tx_frac` | 0.5 | | how much of the transmitted pulse may be missing at the lowest gates; sets the lowest range gate |
| `max_range_delay_us` | 7000 | | the longest delay analysed (µs); sets the highest range |
| `lag_avg` | 1 | | average this many neighbouring lags |
| `save_acf_images` | true | false | write a PNG of the ACFs of every period |
| `reanalyze` | false | | recompute periods whose output exists |
| `periods` | all | the zenith periods | the integration periods to analyse (10-s periods from the start of the recording) |
| `weighting` | `"variance"` | | `"variance"`: each lagged product weighted by its local noise level; `"equal"`: ordinary least squares |
| `background_column` | true | | an extra unknown per lag for a background that is the same at every range |
| `outlier_rejection` | true | false | reject lagged products by a ratio test and a power cut; switch it off with the satellite unknowns, which model what it would reject |
| `impulse_rows` | false | true | detect power-line impulses and leave out every lagged product they touch |
| `weights_exclude_self` | false | true | leave-one-out weights: each product's own power is left out of its weight, which keeps the weights unbiased |
| `weights_width` | 10 | 100 | the number of pulses in the running mean that gives the weights |
| `weights_row_cap` | null | 20 | with leave-one-out weights: a product's variance is at least its own power divided by this, so a single product with strong interference is still weighted down |
| `satellites` | off | see below | satellite echoes as unknowns |

### `steps.lpi.satellites`

| Key | Default | v3 | Meaning |
|---|---|---|---|
| `catalogue` | required | | the satellite catalogue (`pipeline/single_pulse_satellite_cfar.py`, user guide section 5) |
| `channel_delay_us` | measured | | the delay from the transmit sample to the receiver channel (µs); left out or `null`, it is measured from the data |
| `track` | false | true | link the detections into passes and give an unknown to every pulse of a pass, detected or not |
| `track_options` | {} | `{"min_weak_pulses": 3}` | options of the tracking; `min_weak_pulses`: a pass of weak detections only is kept when it holds at least this many pulses |
| `strong_pfa` | — | 1e-10 | with a catalogue made at a higher false-alarm probability: the probability whose threshold marks a detection strong |
| `min_snr_db` | — | | use only echoes at least this strong |
| `enabled` | true | | false switches the section off |

## `steps.fit_lpi`: the fit of the ACFs

Reads the `lpi` step's files (it needs `steps.lpi` in the configuration,
for `range_gate_us`) and writes `pp-<t>.h5` next to them.

| Key | Default | v3 | Meaning |
|---|---|---|---|
| `channel` | required | | as in `lpi` |
| `max_dt` | 300 | | consecutive integration periods starting within this many seconds are fitted together |
| `range_avg` | [1, 3, 5] | [1, 3, 3, 5] | per height band, the ACF is averaged over 2·n + 1 gates (n = 1: 3 gates) |
| `range_limits` | [0, 300, 700, 1500] | [0, 230, 300, 700, 1500] | the heights (km) between which each `range_avg` applies; one more than `range_avg` |
| `range_avg_rho` | null | [-0.30, -0.08, -0.06] | the correlation between the errors of gates 1, 2, 3 apart, for the error of the range average; null: independent |
| `range_avg_rho_min_km` | 0 | 300 | apply `range_avg_rho` only in the bands from this height up |
| `clutter_cutoff_km` | 300 | 150 | below this range the ground-clutter-subtracted ACF is fitted, above it the plain ACF |
| `first_lag` | 0 | | the first lag fitted (index) |
| `max_lag_us` | null | | fit only the lags up to this (µs); null: all |
| `space_object_cut` | null | | leave out gates whose ACF looks like a hard target: `"absolute"` (Doppler width below 300 m/s and its error below 100 m/s) or `"significance"` (width plus three errors below 300 m/s) |
| `fit_select_best` | false | | keep the best of the fit's starting points |
| `fit_restarts` | 0 | | restart the fit this many times from the best start |
| `reanalyze` | false | | refit windows whose output exists |

## `steps.ruisdap_export`: the RUISDAP viewer

After `fit_lpi`, writes the fits one folder per UTC day in the format of the
[RUISDAP](https://github.com/jvierine/ruisdap) Millstone viewer
(`pipeline/ruisdap_export.py`; user guide section 6). The viewer shows each
day next to its plasma lines, which come from RUISDAP's own export of that
day, so only days RUISDAP has exported can be exported.

| Key | Default | Meaning |
|---|---|---|
| `output_dir` | required | the folder the day folders go into |
| `name` | required | the day folder's name, with `{date}` (e.g. `"millstone-{date}-zenith-isr"`); an existing folder is skipped |
| `label` | required | what the run is, shown in the viewer |
| `plasma` | required | per day, `{"YYYY-MM-DD": {"manifest": ..., "url": ...}}`: the day's RUISDAP export's `manifest.json` and the URL of its plasma spectrum file; the days listed are the days exported |
| `magic_constant` | fit_quality's of the day | the density constant |
| `antenna` | `"ZENITH"` | the antenna, as in the plasma manifest |
| `quality_flags` | true | write fit_quality's flags per fit and gate |
| `echo_catalogue` | — | a satellite catalogue for the SATELLITE flag (echoes of at least 20 dB) |
| `source_commit` | — | the isr_analysis commit, for the record |

## `steps.long_pulse`: the long-pulse spectra

Writes `<output_dir>/range_doppler_<mode><postfix>/<channel>/il_<t>.h5`.

| Key | Default | Meaning |
|---|---|---|
| `channel` | required | as in `lpi` |
| `mode` | 300 | the pulse type (`sweepid`) to use |
| `avg_dur` | 10 | seconds averaged into one spectrum |
| `step` | 10 | seconds between spectra |
| `avg_type` | `"outlier_removal"` | how the pulses are averaged: `"outlier_removal"`, `"median"` or `"mean"` |
| `postfix` | `"_outlier"` | appended to the output folder's name |
| `min_tx_pulses` | 100 | the fewest transmitted pulses for a spectrum |
| `tx_delay_us` | measured | the transmit-to-receiver delay (µs); null measures it from the data |
| `reanalyze` | false | recompute spectra whose output exists |

## `steps.fit_lp`: the fit of the long-pulse spectra

Reads the `long_pulse` step's files (it needs `steps.long_pulse` in the
configuration) and writes `pp-<t>.h5` next to them.

| Key | Default | Meaning |
|---|---|---|
| `channel` | required | as in `lpi` |
| `avg_dur` | 600 | seconds of spectra averaged into one fit |
| `ridx` | [35, 230] | the first and last range gate fitted (indices) |
| `fit_bandwidth_hz` | 50000 | fit the spectrum within ± this frequency |
| `notch_bands_hz` | 20–28 kHz | frequency bands left out of the fit, as `[[low, high], ...]` in Hz; `[]` fits the whole band. Interference belongs to a recording: review this per campaign |
| `pulse_length_us` | from the pulse type | the transmitted pulse length (µs), for another experiment |
| `remove_space_objects` | false | leave out gates that look like hard targets |
| `reanalyze` | false | refit windows whose output exists |
