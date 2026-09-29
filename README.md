# ISR Analysis Pipeline

Tools for collective Thomson scatter radar ionospheric plasma-parameter analysis. The pipeline handles space-object contamination, radio-frequency interference, and coded/uncoded long-pulse modes.

![Lag-profile inversion example](figs/lpi_example_2023-09-05.png)

## Quick start

```bash
# 1. Clone the repo
git clone https://github.com/jvierine/isr_analysis.git ~/src/isr_analysis
cd ~/src/isr_analysis

# 2. Create the Python virtual environment (once per machine)
bash setup_env.sh

# 3. Edit or create a config file under config/
cp config/millstone_2023-09-05.json config/my_experiment.json
# ... edit data_dir, output_dir, radar_freq_hz, etc.

# 4. Run (replace 24 with your core count)
mpirun -np 24 ~/venv/isr_analysis/bin/python3 \
    run_analysis.py config/my_experiment.json
```

## Setup details

### Python environment

`setup_env.sh` creates `~/venv/isr_analysis` with `--system-site-packages` so
it inherits system-installed numpy/scipy/matplotlib/h5py, then pip-installs:

| Package | Purpose |
|---------|---------|
| `pyfftw` | Fast FFT via FFTW3 (needs `libfftw3-dev`) |
| `digital_rf` | Read raw DigitalRF HDF5 voltage data |
| `mpi4py` | MPI parallelisation |
| `jcoord` | Geographic coordinate conversion |

Activate the environment (optional, only needed for interactive use):
```bash
source ~/venv/isr_analysis/bin/activate
```

### ISR theory lookup tables

Plasma-parameter fitting requires a precomputed spectral interpolation table.
The table is cached in `./data/` and regenerated automatically if it does not
exist or if you change `radar_freq_hz` in the config.  The `fit_lpi`, `fit_lp`,
and legacy `fit_ionline` paths all use this same lazy loader.  When launched
under `mpirun`, table generation is collective and the spectra are divided
among the MPI ranks; rank 0 writes the completed cache after an all-rank
reduction. Generation takes ~10–30 minutes depending on the machine and rank
count; subsequent runs load the cached file.

Override the cache directory with `"table_dir": "/path/to/tables"` in your
JSON config.

## Configuration

Config files live in `config/`. Key fields:

```json
{
  "experiment":    "my_experiment",
  "data_dir":      "/path/to/raw/digitalrf/experiment",
  "output_dir":    "/path/to/results",         // optional
  "max_time_s":    600,                          // omit for full dataset
  "radar_freq_hz": 440200000,                    // default 440.2 MHz
  "table_dir":     "./data",                     // optional, default ./data

  "steps": {
    "lpi":      { "enabled": true,  "channel": "zenith-l", "range_gate_us": 60, ... },
    "fit_lpi":  { "enabled": true,  "channel": "zenith-l", "max_dt": 300, ... },
    "long_pulse":{ "enabled": false, ... },
    "fit_lp":   { "enabled": false, ... }
  }
}
```

See `config/millstone_2023-09-05.json` for a complete annotated example.

## Analysis modes

### Coded long pulse (LPI)

High-range-resolution bottom-side analysis. Run in order:

1. `outlier_lpi.py` — lag-profile inversion; estimates ACFs per range gate
2. `fit_lpi.py` — fits ACFs to Te, Ti, vi, ne profiles

Driven by the `lpi` and `fit_lpi` steps in the config.

### Uncoded long pulse

Low range resolution, optimised for the topside where SNR is low. Run in order:

1. `avg_range_doppler_spec.py` — range-Doppler spectral averaging
2. `fit_lp.py` — fits Doppler spectra to Te, Ti, vi, ne

Driven by the `long_pulse` and `fit_lp` steps in the config.

### Single-pulse satellite CFAR

`single_pulse_satellite_cfar.py` searches every eligible pulse with the
simultaneously measured `tx-h` waveform.  It dispatches sweeps 1--32, 300, and
800 to their correct gates, selects `zenith-l` or `misa-l` from the recorded
transmit/receive antenna state, and coherently sums non-overlapping groups of
eight complex64 voltage samples before decimation.  The range search therefore
has eight times fewer gates (1.199 km spacing), while the pulse-time keys and
the 11-sample receiver-delay correction remain on the original 1 MHz clock.
The detector uses complex64 FFTW transforms with at least 4x waveform zero
padding and applies two-dimensional CA-CFAR.  Separated echoes in one pulse are
retained.

The output is a compact Digital Metadata channel.  Its keys are the original
absolute pulse samples (Unix microseconds at 1 MHz), and each record contains
arrays of detected range, Doppler, SNR, CFAR, and provenance values.  MPI ranks
first write independent compressed HDF5 chunks, making an interrupted run
restartable; rank 0 creates the metadata channel only when all chunks exist.

```bash
mpirun -np 50 ~/venv/isr_analysis/bin/python3 \
    single_pulse_satellite_cfar.py \
    --data /path/to/usrp-rx0-r_20240407T100000_20240409T110000 \
    --output /path/to/usrp-rx0-r_20240407T100000_20240409T110000/metadata/satellite_detections \
    --receiver-delay-samples 11 --decimation-factor 8 --fft-padding 4 \
    --max-doppler-hz 60000
```

The fixed receiver delay is subtracted from the receive-window start before
range conversion; the transmit array is not circularly rolled.  Positive
Doppler denotes increasing complex baseband phase.  Catalogue association is
performed as a later stage so detections remain independent of TLE coverage.

## Output files

Results go to `output_dir` (or alongside the raw data if unset):

```
lpi_<rg_us>/zenith-l/
  lpi-<unix_t>.h5    ACF per range gate and lag
  lpi-<unix_t>.png   diagnostic image

lpi_<rg_us>/zenith-l/
  pp-<unix_t>.h5     Te, Ti, vi, ne profiles
  pp-<unix_t>.png    diagnostic plot
```

## Plotting raw voltage

`plot_raw_voltage.py` plots the undecoded complex voltage samples one figure per
interpulse period, using the transmit pulse metadata (`metadata/id_metadata`) to
find the IPP boundaries. Useful for checking timing, the transmit pulse and
ground clutter extent, the noise injection window, and RFI before running any
analysis.

The dataset comes from `config/millstone_eclipse2024.json`, read automatically
(the `CONFIG` constant at the top of the script picks the file): `data_dir` says
where the raw data is and `output_dir` is where the plots go. There is nothing
to pass on the command line.

```bash
# first 10 IPPs of every channel
python3 plot_raw_voltage.py

# 20 IPPs of one channel, starting 100 IPPs into the recording
python3 plot_raw_voltage.py --channel zenith-l --n-ipp 20 --start-ipp 100
```

Two PNGs are written per IPP per channel, with the transmit pulse (red) and
noise injection (green) windows shaded and the sweepid in the title:

```
<output_dir>/IPP/<channel>/raw voltage/ipp-<index>-<unix_t>.png   |z|, real, imaginary (linear)
<output_dir>/IPP/<channel>/power/ipp-<index>-<unix_t>.png         10*log10(|z|^2), dB
```

When more than one channel is plotted, the power of all of them is additionally
drawn in one figure per IPP, for comparing the transmit samples against the echo
channel:

```
<output_dir>/IPP/combined/ipp-<index>-<unix_t>.png                all channels, dB
<output_dir>/IPP/combined/close_up/ipp-<index>-<unix_t>.png       the same, transmit pulse only
<output_dir>/IPP/combined/tx_delay-<unix_t>.png                   matched filter delay vs IPP
```

The close ups isolate the transmit pulse (the `tx0`–`tx1` window of the sweepid,
plus a small margin). Those isolated pulses are also matched filtered against
each other:

```
xc(tau) = sum_t z_tx(t) conj(z_echo(t - tau))
```

that is, the transmit channel convolved with the time reversed complex conjugate
of the echo channel. The cross correlation spectrum is zero padded before the
inverse transform, which band limited interpolates `xc` onto a 1/`--oversample`
sample grid, so the lag of the peak resolves the transmit pulse delay between
the two channels to a fraction of a sample. `--max-lag-us` bounds the search.
The delay is printed per IPP and plotted against IPP number. It is per
experiment: 10.50 us for the 2023-09-05 data, and for eclipse2024 12.39 us
(misa-l) and 11.98 us (zenith-l). Interpulse periods transmitting the uncoded
long pulse (sweepid 300) give a slightly different value than the coded ones,
because the autocorrelation of an unmodulated pulse has no sharp peak to
locate.

`--start-ipp` and `--t0` select where to start, `--outdir` overrides the output
location from the config, and `--remove-dc` subtracts the known USRP DC offset.

`--fix-tx-delay` shifts each echo channel onto the `tx-h` time base by its
measured channel delay (see *Transmit channel delay* below) and writes the same
set of plots to `IPP_offset_fix/` instead of `IPP/`, so the two can be compared
directly. The delay plot then shows the residual, which sits at zero for the
coded pulses; the uncoded long pulse interpulse periods keep a ~0.2 us residual
because the estimator cannot locate their broad correlation peak, not because
they are corrected any differently.

```bash
python3 plot_raw_voltage.py --fix-tx-delay
```

## Transmit channel delay

`tx-h` and the echo channels (`zenith-l`, `misa-l`) are digitised by separate
receiver chains, each with its own filter and decimation group delay, so a
transmit waveform read from `tx-h` is offset by 10-12 us from the same pulse as
it appears in an echo channel. That is far too long to be cable length: it is
dominated by the DDC filter chain, and it differs per channel and per
experiment.

`tx_delay.py` measures it, by matched filtering `tx-h` against the transmit
pulse leaking into the echo channel:

```python
import tx_delay as txd
delay_us, spread_us, n_pulses = txd.estimate_channel_delay(datadir, "misa-l")
```

Only the phase coded pulses are used. An uncoded long pulse autocorrelates to a
broad triangle with no sharp peak to locate and gives a delay biased by ~0.2 us;
the delay is a property of the receiver chain, not of the waveform, so the coded
pulse estimate is the one to apply to every mode. Measured values:

| Experiment | Channel | Delay |
|---|---|---|
| 2023-09-05 | zenith-l | 10.49 us |
| eclipse2024 | misa-l | 12.34 us |
| eclipse2024 | zenith-l | 11.88 us |

`avg_range_doppler_spec.py` uses this to put the transmit waveform onto the echo
channel time base before estimating the range-Doppler ambiguity function. It
used to apply a hardcoded `n.roll(z_tx, 11)`, eyeballed from the leakthrough and
quantised to whole samples (1 us, 150 m in range); the correction is now
measured per run and per channel using the first 100 coded pulses, and applied as a band limited interpolation so
it is not quantised. The value used is stored as `tx_delay_us` in each output
file. Set `"tx_delay_us"` in the `long_pulse` config step to override the
measurement with a fixed value.

`outlier_lpi.py` needs no such correction: it takes its transmit envelope from
the pulse leaking into the echo channel itself, so the delay cancels.

This is a software workaround. The real fix requires a hardware modification: Interleave the transmit sample
into the echo channel with an analog switch, so that both pass through one
receiver chain and no correction is needed at all.

## Reading the raw data

The raw voltage is one HDF5 file per second per channel, on spinning disks.
Reading it pulse by pulse through `DigitalRFReader` costs several small reads
per pulse, and is limited by disk seeks. `raw_reader.py` reads whole seconds
instead, with one read each, and returns identical samples:

```python
from raw_reader import RawReader, pulse_index
r = RawReader("%s/rf_data/" % dirname)               # r.read_vector_1d(key, 10000, "zenith-l")
keys = pulse_index(dirname)["key"]                    # every pulse, cached after the first call
```

`pulse_index` builds a list of all pulses once (about 10 minutes for the
eclipse recording) in `~/.cache/isr_analysis/`, and later loads it in under a
second. For scripts that read a lot of data, sort the work by time, give each
process a contiguous stretch, and let 8 to 16 processes read at once; more
processes read slower. `outlier_lpi.py` uses the reader by default
(`fast_read=True`), which changes little there because the inversion is
limited by computation. `validate_raw_reader.py` checks that the output is bit
identical.

## Plotting plasma parameters

After the fitting step has produced `pp-*.h5` files, use `plot_pp.py` to
produce a 4-panel time-vs-range summary (Te, Ti, vi, Ne):

```bash
python3 plot_pp.py <pp_directory> [channel]
```

For example, to plot the LPI fit results from the 2023-09-05 run:

```bash
python3 plot_pp.py \
    /mnt/data/juha/millstone_hill/results/2023-09-05/lpi_60/zenith-l \
    zenith-l
```

The script saves a PNG (`ppar-<channel>-<start_time>.png`) and an HDF5 file
(`ppar-<channel>-<start_time>.h5`) in the current directory with merged
time-series arrays (range, Te, Ti, vi, ne, uncertainties, space-object masks).

**Electron density calibration** — the raw `ne` values in the `pp-*.h5` files
are uncalibrated. A `magic_constant` (default `75800958.63 × 10`) converts
them to absolute density in m⁻³. If you have run the calibration workflow
(using `plasma_line_clicker.py` and `estimate_magic_constant.py`), place
`magic_const.h5` in the same directory as the `pp-*.h5` files and it will be
picked up automatically.

## System temperature against the 440 MHz sky

Three tools compare the measured system noise temperature with a global sky
model, which tests the noise injection calibration end to end.

```bash
python3 tsys_harvest.py config/millstone_eclipse2024.json -o tsys.h5
python3 compare_tsys_sky.py tsys.h5 -o tsys_sky/
```

`tsys_harvest.py` collects the `T_sys` that every integration period already
carries, from all analysis output sets found under `data_dir`, and attaches the
antenna pointing, which the analysis does not record because the plasma
parameter fits do not need it.

`sky_noise_model.py` evaluates PyGDSM at the radar frequency. The default is
GSM2008 on its `haslam` basemap, which is locked to the Haslam 408 MHz survey
and carries its 1 degree resolution, so 440 MHz is an 8 per cent extrapolation
from the anchor map. The CMB is included: PyGDSM leaves it out by default
because its models describe galactic and extragalactic emission, but an antenna
does not know the difference. Run the module directly to see how much sky signal
there is to find.

`compare_tsys_sky.py` masks periods when a bright source is in the beam, rejects
interference by its shape in time rather than the size of its residual, and
fits `T_sys = T_0 + eta * T_sky`.

### Result, eclipse2024, zenith antenna

| | value |
|---|---|
| `T_0` | 156.2 +- 1.6 K |
| `eta` | 0.661 +- 0.055 |
| variance explained | 0.76 |

`T_0` is receiver plus spillover plus atmosphere. `eta` is **not** beam
efficiency alone: `T_sys` is derived as `noise/alpha` with `alpha` proportional
to `1/T_INJECTION`, so the whole series scales with the assumed 1172 K, and this
fit cannot separate a beam efficiency of 0.66 from a `T_INJECTION` that should
be 1773 K. Uncertainties are from a 30 minute moving block bootstrap; the least
squares covariance is optimistic by more than an order of magnitude here,
because consecutive ten second periods see the same sky.

Two independent estimators agree on both numbers to 1 per cent:
`range_doppler_300_outlier` gives the values above and `lpi_30`, once its bias
is undone (see `lpi_bias_factor`), gives `T_0 = 155.3 +- 1.7 K` and
`eta = 0.654 +- 0.060`.

### Two things this exposed

**The zenith antenna is not at the zenith.** Its boresight is at elevation
88.16, azimuth 172.9, which puts it at declination 40.79 against Cygnus A's
40.734. Cyg A therefore transits the beam centre once per sidereal day and takes
`T_sys` from 170 K to 1600 K. This reproduces, from a radio source, the 1.84
degree offset the hard target work measured from satellite echoes.

**MISA cannot be fitted from the recorded metadata.** Its azimuth alternates
between integration periods and the antenna control metadata does not locate the
switches well enough to say which azimuth a ten second period belongs to. It is
reported but flagged as not a measurement.

> Code is still under active development.
