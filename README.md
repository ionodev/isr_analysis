# ISR analysis

This repository turns raw incoherent scatter radar (ISR) recordings from the
Millstone Hill observatory (Massachusetts, USA) into profiles of the
ionosphere: electron density, electron and ion temperature, and ion drift
velocity against altitude and time. The main dataset is a 49-hour recording,
7–9 April 2024, which spans the total solar eclipse of 8 April 2024.

It is a fork of Juha Vierinen's
[isr_analysis](https://github.com/jvierine/isr_analysis), extended in a
student project.

![Lag-profile inversion example](figs/lpi_example_2023-09-05.png)

## What the radar measures

The radar transmits pulses at 440.2 MHz. A very small fraction of each pulse
is scattered back by thermal fluctuations of the electron density in the
ionosphere. The strength of this echo gives the electron density. Its
spectrum, or equivalently its autocorrelation function (ACF), gives the
temperatures and the drift velocity.

The recording holds the complex baseband voltages of three channels, in
[Digital RF](https://github.com/MITHaystack/digital_rf) format:

| Channel    | What it is |
|------------|------------|
| `zenith-l` | receiver of the fixed, vertically pointing 68 m antenna |
| `misa-l`   | receiver of the steerable 46 m MISA antenna |
| `tx-h`     | a sample of the transmitted pulse |

The experiment interleaves several pulse types ("modes"):

- **coded pulses** (sweeps 1–32): fine range resolution low in the ionosphere;
- **an uncoded long pulse** (mode 300): more signal high up, at coarser resolution;
- **mode 800**: MISA's low-elevation scans.

At the end of every pulse a noise diode of known temperature is switched on,
which calibrates the receiver.

## How it works

```
raw voltages (Digital RF)
 ├─ coded pulses ─► outlier_lpi.py ─► ACF per range gate ─► fit_lpi.py ─┐
 │                  (lag-profile inversion)                              ├─► ne, Te, Ti, vi profiles
 └─ long pulse ───► avg_range_doppler_spec.py ─► spectra ─► fit_lp.py ──┘
```

1. **Calibration.** The noise diode's window measures the receiver gain, which
   converts the echo power to a system temperature and then to density.
2. **Lag-profile inversion** (`outlier_lpi.py`). The coded pulses are
   decoded into an ACF for each 60 µs range gate, by solving a linear
   least-squares problem over every pulse of a 10 s period.
3. **Range–Doppler spectra** (`avg_range_doppler_spec.py`). The same for the
   long pulse, as averaged power spectra.
4. **Fitting** (`fit_lpi.py`, `fit_lp.py`). A theoretical ISR spectrum
   (`isr_spec.py`, precomputed into tables by `il_interp.py`) is fitted to each
   ACF or spectrum.

Unwanted signals are handled on the way:

- **Satellites and debris.** `single_pulse_satellite_cfar.py` searches every
  pulse for bright point targets and writes a catalogue. `satellite_columns.py`
  lets the inversion fit catalogued echoes as extra unknowns, instead of
  discarding the affected data.
- **Outliers** in the lagged products are rejected statistically. Ground
  clutter is removed pulse to pulse.
- **Hard targets in the fit.** `fit_lpi.py` can also leave out range gates
  whose ACF is too narrow for plasma (`space_object_cut` in the `fit_lpi`
  section of the configuration). It is off by default, since the catalogue and
  the satellite columns handle satellites earlier. `"absolute"` is the
  original test (Doppler width below 300 m/s and its uncertainty below
  100 m/s); `"significance"` asks for the width to be below 300 m/s at three
  standard deviations.

`run_analysis.py` runs the steps enabled in a JSON configuration file; the
modules it uses are in `pipeline/`.

## Quick start

Requirements: Python 3, the FFTW library (`libfftw3-dev`) and an MPI
implementation (e.g. OpenMPI); the Python packages are in `requirements.txt`.

```bash
git clone https://github.com/ionodev/isr_analysis.git
cd isr_analysis
bash tools/setup_env.sh                 # creates ~/venv/isr_analysis from requirements.txt
cp config/example_basic.json config/my_run.json
# edit data_dir, output_dir, channel and the steps to run
mpirun -np 24 ~/venv/isr_analysis/bin/python3 run_analysis.py config/my_run.json
```

The first run builds the theory lookup tables in `./data/` (10–30 minutes);
later runs reuse them. `docs/appendix_analysis_guide.tex` explains how to get an
account on the analysis server and run the pipeline there.

### With satellite and impulse removal

`config/example_v3.json` has the settings of pipeline version 3: satellite
echoes as unknowns with pass tracking, the rows hit by power-line impulses
removed, and leave-one-out weights. It needs two things first:

```bash
# 1. the satellite catalogue: every pulse searched for echoes
python3 pipeline/single_pulse_satellite_cfar.py --data /path/to/data --output /path/to/satellite_catalogue --pfa 1e-6
# 2. the wider theory table the version 3 fit uses (written to data/v3/wide)
mpirun -np 24 python3 tools/make_table.py data/v3
```

Then set `data_dir`, `output_dir` and `steps.lpi.satellites.catalogue` in a
copy of the example and run it as above.

### Tests

```bash
pip install -r requirements-optional.txt
python3 -m pytest tests
```

The tests also run on GitHub for every push and pull request.

## Repository layout

`run_analysis.py`, the script that runs everything, is the only script in
the top folder.

| Folder | Purpose |
|--------|---------|
| `pipeline/` | the modules `run_analysis.py` uses (below) |
| `calibration/` | the density scale: plasma-line picking and the calibration constant |
| `analysis/` | evaluation, test and figure scripts for the results (not part of the pipeline; see its README) |
| `plotting/` | diagnostic plots of raw data, ACFs and profiles |
| `tests/` | unit tests (`python3 -m pytest` in the top folder) |
| `tools/` | environment setup and the regression benchmark |
| `config/` | run configurations |
| `data/` | the theory lookup tables (built on the first run) |
| `docs/`, `figs/` | the analysis-server guide; figures |

The modules in `pipeline/`:

| Files | Purpose |
|-------|---------|
| `outlier_lpi.py`, `fit_lpi.py` | coded pulses: inversion and fit |
| `avg_range_doppler_spec.py`, `fit_lp.py` | long pulse: spectra and fit |
| `isr_spec.py`, `il_interp.py`, `fit_ionline.py` | ISR theory, its lookup tables, ion-line model helpers |
| `radar_timing.py`, `millstone_radar_state.py`, `tx_delay.py`, `tx_power.py` | pulse timing, antenna pointing, transmitter state |
| `raw_reader.py`, `stuffr.py` | fast reader for the raw data; helpers |
| `single_pulse_satellite_cfar.py`, `satellite_columns.py` | satellite detection, tracking and fitting |
| `impulse_rows.py` | power-line impulses: the rows they hit, removed from the inversion |
| `fit_quality.py` | density constants and quality flags of the fitted products |

Scripts outside `pipeline/` find its modules by adding the folder to their
import path, so they run from anywhere, e.g. `python3 analysis/eclipse_plots.py`.

## Documentation

Methods, tests and results are written up in a series of processing memos,
with a running task list. These are kept outside this repository and are
available to project contributors.

## Reusable range–Doppler processing

The Rust [range–Doppler estimator](crates/range-doppler) handles configurable
within-pulse spectra, overlapping fast-time gates, UTC integration windows and
signed-spectrum normalization. It has no instrument file-format dependencies.
The [SYISR uncoded-long-pulse driver](syisr/lower_hybrid) handles SYISR binary
records, experiment checks, beam geometry, transient screening and HDF5 products.
Animation export uses Rustplotlib's reusable streaming `AnimationWriter` API.

The [SYISR 2 ms ion-line profile analysis](syisr/ulp2ms/) provides raw pulse
screening, signed lag products, full range/lag covariance, and joint
5000-km profile inversion with either knot composition or a global
sigmoid-H/Gaussian-He model. Scientific products remain external HDF5 files.
