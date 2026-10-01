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
- **Power-line spark impulses** can be blanked (`impulse_blanking.py`, off by
  default).

`run_analysis.py` runs the steps enabled in a JSON configuration file.

## Quick start

```bash
git clone https://github.com/ionodev/isr_analysis.git
cd isr_analysis
bash setup_env.sh                       # creates ~/venv/isr_analysis
cp config/millstone_eclipse2024.json config/my_run.json
# edit data_dir, output_dir, channel and the steps to run
mpirun -np 24 ~/venv/isr_analysis/bin/python3 run_analysis.py config/my_run.json
```

The first run builds the theory lookup tables in `./data/` (10–30 minutes);
later runs reuse them. `appendix_analysis_guide.tex` explains how to get an
account on the analysis server and run the pipeline there.

## Repository layout

| Files | Purpose |
|-------|---------|
| `run_analysis.py`, `config/` | pipeline driver and run configurations |
| `outlier_lpi.py`, `fit_lpi.py` | coded pulses: inversion and fit |
| `avg_range_doppler_spec.py`, `fit_lp.py` | long pulse: spectra and fit |
| `isr_spec.py`, `il_interp.py` | ISR theory and its lookup tables |
| `radar_timing.py`, `millstone_radar_state.py`, `tx_delay.py`, `tx_power.py` | pulse timing, antenna pointing, transmitter state |
| `raw_reader.py` | fast reader for the raw data |
| `single_pulse_satellite_cfar.py`, `satellite_columns.py` | satellite detection and fitting |
| `impulse_blanking.py`, `sky_noise_model.py` | interference blanking, sky noise |
| `plot_*.py` | diagnostic plots of raw data, ACFs and profiles |
| `validate_*.py` | tests of the methods on real and injected data |
| `test_*.py` | unit tests (`python3 -m pytest`) |

## Documentation

Methods, tests and results are written up in a series of processing memos,
with a running task list. These are kept outside this repository and are
available to project contributors.
