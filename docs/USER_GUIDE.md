# User guide

This guide explains how to install the pipeline, run it on a Millstone Hill
recording, and read what it writes. The [README](../README.md) gives the
overview and the folder layout.

- [1. What the pipeline does](#1-what-the-pipeline-does)
- [2. Installation](#2-installation)
- [3. The input data](#3-the-input-data)
- [4. Running the pipeline](#4-running-the-pipeline)
- [5. Satellite and impulse removal (version 3)](#5-satellite-and-impulse-removal-version-3)
- [6. The output files](#6-the-output-files)
- [7. Reading the results in Python](#7-reading-the-results-in-python)
- [8. The density scale](#8-the-density-scale)
- [9. Known limitations](#9-known-limitations)

## 1. What the pipeline does

The pipeline reads the raw receiver voltage of an incoherent scatter radar
(ISR) and estimates, for every range gate and time window, the electron
density $n_e$, the electron and ion temperatures $T_e$ and $T_i$, and the
line-of-sight ion velocity $v_i$.

`run_analysis.py` runs up to four steps, each switched on in the
configuration:

| Step | Module | What it does | Writes |
|---|---|---|---|
| `lpi` | `pipeline/outlier_lpi.py` | lag-profile inversion of the coded pulses: the autocorrelation function (ACF) at every range gate, every integration period (10 s) | `lpi-*.h5` |
| `fit_lpi` | `pipeline/fit_lpi.py` | averages the ACFs in time and range and fits a theoretical ion-line ACF to them | `pp-*.h5` |
| `long_pulse` | `pipeline/avg_range_doppler_spec.py` | range–Doppler power spectra of the uncoded long pulse | `il_*.h5` |
| `fit_lp` | `pipeline/fit_lp.py` | fits a theoretical spectrum to the long-pulse spectra | `pp-*.h5` |
| `ruisdap_export` | `pipeline/ruisdap_export.py` | the fits in the format of the RUISDAP web viewer, one folder per day | `parameters.f32`, `manifest.json`, ... |

The two branches are independent: the coded pulses (`lpi`, `fit_lpi`) give
the finer range resolution and are the main product; the long pulse
(`long_pulse`, `fit_lp`) is optional. The satellite and impulse removal
(section 5) work inside the `lpi` step.

## 2. Installation

The pipeline needs Python 3, the FFTW library and an MPI implementation:

```bash
sudo apt-get install libfftw3-dev libopenmpi-dev openmpi-bin   # Debian/Ubuntu
git clone https://github.com/ionodev/isr_analysis.git
cd isr_analysis
bash tools/setup_env.sh          # a virtual environment in ~/venv/isr_analysis
source ~/venv/isr_analysis/bin/activate
```

`setup_env.sh` installs the packages of `requirements.txt` into the virtual
environment and takes numpy, scipy, matplotlib and h5py from the system where
it has them. Without the script: `pip install -r requirements.txt`. The tests
and a few scripts outside the pipeline need `requirements-optional.txt` as
well.

**Theory tables.** The fit compares the measured ACFs with precomputed
theoretical spectra. The first run builds these tables in `data/` (10–30
minutes, in parallel under MPI); later runs reuse them. Version 3 uses a
wider table that is built separately (section 5).

Check the installation with the tests:

```bash
pip install -r requirements-optional.txt
python3 -m pytest tests
```

Seven tests need the theory tables and are skipped until they exist.

## 3. The input data

The input is a recording in the [Digital RF](https://github.com/MITHaystack/digital_rf)
format as Millstone Hill writes it: one folder, `data_dir`, that holds

| Path | Contents |
|---|---|
| `rf_data/zenith-l/` | complex voltage of the zenith antenna's receiver, 1 MHz |
| `rf_data/misa-l/` | the same for the steerable MISA antenna |
| `rf_data/tx-h/` | a sample of every transmitted pulse |
| `metadata/id_metadata/` | the pulse type of every pulse (`sweepid`: 1–32 alternating code, 300 long pulse, 800 other) |
| `metadata/powermeter/` | the transmitter power |
| `metadata/antenna_control_metadata/` | which antenna transmits, and its pointing |

The receiver is calibrated from the noise diode switched on at the end of
every pulse, inside the recording; no separate calibration file is needed
for the inversion and the fit.

## 4. Running the pipeline

### The configuration

A run is described by one JSON file. Start from an example:

| Example | What it runs |
|---|---|
| `config/example_basic.json` | the plain pipeline: inversion and fit of the coded pulses, the long pulse switched off |
| `config/example_v3.json` | version 3: the same with satellite and impulse removal (section 5) |

```bash
cp config/example_basic.json config/my_run.json
```

and set at least:

- `data_dir`: the recording (section 3);
- `output_dir`: where the results go (default: inside `data_dir`);
- `channel` in each step: `zenith-l` or `misa-l`.

Other top-level settings: `max_time_s` (analyse only the first this many
seconds, for a quick test), `radar_freq_hz` (default 440.2 MHz), `table_dir`
(the theory tables; default `data/`). Keys starting with `_` (such as
`_notes`) are ignored. Each step's settings are listed in
[the configuration reference](CONFIGURATION.md).

### Running

```bash
mpirun -np 24 python3 run_analysis.py config/my_run.json
```

The steps run in the order `lpi`, `fit_lpi`, `ruisdap_export`, `long_pulse`, `fit_lp`; the
work of each step is shared among the MPI processes. Without `mpirun` the
pipeline runs in one process. Use as many processes as the machine has
cores, and keep each process single-threaded:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
```

**Resuming.** A step skips the time windows whose output already exists,
so an interrupted run continues where it stopped. Set `"reanalyze": true` in
a step to compute everything again.

**Part of a recording.** `steps.lpi.periods` takes a list of integration
period numbers (10-s periods counted from the start of the recording) and
analyses only those, e.g. the times when the zenith antenna transmits.

## 5. Satellite and impulse removal (version 3)

Version 3 detects satellite echoes and power-line impulses in every pulse
and removes them inside the inversion:

- **satellite echoes** are found by a separate detector run, which writes a
  catalogue; the `lpi` step links the detections into passes and gives every
  echo its own unknown in the inversion;
- **impulses** are found in the voltage by the `lpi` step itself, and the
  measurements (lagged products) they touch are left out;
- **leave-one-out weights** keep the inversion's weights unbiased when
  measurements are left out.

`config/example_v3.json` has these settings. Two things are needed before
it runs:

**1. The satellite catalogue.** The detector searches every pulse for echoes
with a matched filter and a constant false-alarm-rate test, and writes a
Digital Metadata channel:

```bash
mpirun -np 24 python3 pipeline/single_pulse_satellite_cfar.py \
    --data /path/to/data_dir --output /path/to/satellite_catalogue --pfa 1e-6
```

`--pfa` is the false-alarm probability per searched cell (version 3:
`1e-6`; the default, `1e-10`, keeps only strong echoes). The run can be
restarted; `--help` lists the other options. Put the catalogue's path in
`steps.lpi.satellites.catalogue`.

**2. The wider theory table.** Version 3 fits $T_e/T_i$ up to 5 on a finer
grid:

```bash
mpirun -np 24 python3 tools/make_table.py data/v3       # writes data/v3/wide
```

and the example's `table_dir` is `data/v3/wide`.

`steps.lpi.satellites.channel_delay_us` is the delay between the transmit
sample and the receiver channel; `null` measures it from the data.

## 6. The output files

Each step writes into a folder under `output_dir`:

```
output_dir/
  lpi_60/zenith-l/              lpi step (60 = range_gate_us) and fit_lpi step
    lpi-<t>.h5, lpi-<t>.png     one per integration period (10 s); <t> its start, Unix time
    pp-<t>.h5, pp-<t>.png       one per fit window (max_dt)
    pp_fit_<t>.png              the fitted ACFs of that window
  range_doppler_300_outlier/zenith-l/   long_pulse and fit_lp steps
    il_<t>.h5, il_<t>.png       long-pulse spectra
    pp-<t>.h5, pp_lp_<t>.png    long-pulse fits
```

### The fits: `pp-<t>.h5`

One file per fit window, one value per range gate:

| Key | Meaning | Unit |
|---|---|---|
| `rgs` | range of each gate | km |
| `t0`, `t1` | start and end of the window | Unix s |
| `ne`, `dne` | electron density and its standard error, **relative** (section 8) | arbitrary |
| `Te`, `dTe` | electron temperature | K |
| `Ti`, `dTi` | ion temperature | K |
| `vi`, `dvi` | line-of-sight ion velocity (same sign convention as Madrigal's, positive away from the radar; not independently verified) | m/s |
| `dTe_Ti` | standard error of the fitted $T_e/T_i$ | |
| `O_frac` | O⁺ fraction of the O⁺/H⁺ fit above 700 km (NaN below) | |
| `Sigma`, `Sigma_params` | full covariance of the fitted parameters ($T_e/T_i$, $T_i$, $v_i$, zero-lag power) per gate | |
| `ne_const` | `ne` = zero-lag power × `ne_const`, for propagating `Sigma` to `ne` | |
| `te_ti_max` | the upper $T_e/T_i$ of the theory table used | |

Gates where the fit failed are NaN.

### The ACFs: `lpi-<t>.h5`

One file per 10-s integration period, the input of the fit:

| Key | Meaning |
|---|---|
| `rgs_km`, `lags` | range of each gate (km), lag of each column (s) |
| `acfs_e`, `acfs_var` | the ACF and its variance |
| `acfs_g`, `acfs_var_g` | the ACF after ground-clutter subtraction, and its variance (the fit uses it below `clutter_cutoff_km`) |
| `T_sys`, `alpha`, `P_tx` | system temperature (K); the noise diode's scale from receiver power to temperature; transmit power (W) |
| `sat_*` | the satellite unknowns, when enabled: each echo's pulse (`sat_keys`), delay, Doppler shift and fitted amplitude per lag |
| `impulse_*` | the impulse removal, when enabled: thresholds and the share of rows removed per lag |
| `t0`, `t1` | start and end of the period (Unix s) |

### For the RUISDAP viewer

The step `ruisdap_export` writes the fits in the format of the
[RUISDAP](https://github.com/jvierine/ruisdap) web viewer, one folder per
day; `tools/export_ruisdap.py` does the same for a finished run. The viewer
shows each day next to its plasma lines, from RUISDAP's own export of that
day, so the step needs that export (see the configuration reference). A new
folder appears in the viewer once it is added to the viewer's list of
analyses on its server.

## 7. Reading the results in Python

`pipeline/fit_quality.py` reads the fits of a run into arrays and marks the
gates that should not be used:

```python
import sys
sys.path.insert(0, "pipeline")
import fit_quality as fq

F = fq.load_fits("/path/to/output_dir")     # reads lpi_60/zenith-l/pp-*.h5
# other settings: fq.load_fits(d, subdir="lpi_30/misa-l"); another recording's
# density constant: fq.load_fits(d, constant=6.2e8)
# F["t0"], F["t1"]: window times; F["r"]: ranges (km);
# F["ne"], F["Te"], F["Ti"], F["vi"] and their errors "dne", ...: [window, gate]
flags = fq.flags(F)                         # one bit field per gate
good = fq.good(flags)                       # True where the fit can be used
```

`load_fits` scales `ne` to m⁻³ by the density constant (section 8). A gate
is **good** when it is not:

- **FAILED**: a parameter or its uncertainty is not finite;
- **BOUND**: the fit ended at a limit of its search ($T_e/T_i$ at its
  lowest or highest value, $T_i$ at 150 or 3000 K, $|v_i|$ at 1500 m/s);
- **NOISY**: $dn_e/n_e$ or $dT_e/T_e$ above 0.5.

`flags` can also mark, for information, the gates near catalogued satellite
echoes (`echoes=`) and near satellite unknowns that left residual power
(`columns=`, from `fq.residual_columns`).

## 8. The density scale

The fit gives $n_e$ in relative units: the echo power scaled for range and
transmit power. A **density constant** $C$ turns it into m⁻³:
$n_e = C \times$ `ne`. `fit_quality.load_fits` applies it.

`fit_quality.DENSITY_CONSTANT_BY_DAY` holds the constants measured for the
recording of 7–8 April 2024 from the plasma lines, and `DENSITY_CONSTANT`
the value used for every other day. **For another recording, or other
processing settings, measure the constant again**, and pass it to
`load_fits(..., constant=C)`; it depends on the radar's state and on the
processing. The tools are in `calibration/`:

1. `plasma_line_picker.py` finds the plasma-line ridge in the recording's
   plasma-line spectra automatically (`plasma_line_clicker.py` lets you
   click it by hand instead);
2. `plasma_line_density_scale.py` turns the picks into densities and
   compares them with the fitted peak density, rejecting outliers, to give
   $C$ with its uncertainty; it can also compare with a digisonde.

Each script's usage is in its docstring.

## 9. Known limitations

- The density constants in `fit_quality.py` are for 7–8 April 2024
  (section 8); `fit_quality.residual_columns` numbers the integration
  periods from the start of that recording.
- The scripts in `analysis/` reproduce the project's evaluations on the
  April 2024 recording; they find its data through environment variables
  (`analysis/README.md`).
- Satellite echoes too faint for the detector are not removed.
