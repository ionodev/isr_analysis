# isr_analysis

Incoherent scatter radar (ISR) analysis for the Millstone Hill radar, from the
raw receiver voltage to profiles of electron density, electron and ion
temperature and ion velocity. Satellite echoes and power-line interference are
removed inside the lag-profile inversion.

Python, parallel with MPI; reads Digital RF recordings and writes HDF5. A fork
of Juha Vierinen's [isr_analysis](https://github.com/jvierine/isr_analysis),
extended in a student project at UiT The Arctic University of Norway.

[User guide](docs/USER_GUIDE.md) ·
[Configuration reference](docs/CONFIGURATION.md) ·
[Changelog](CHANGELOG.md)

![Lag-profile inversion example](figs/lpi_example_2023-09-05.png)

The ACF of one 10-second integration period against range and lag, from the
lag-profile inversion (Millstone Hill, 5 September 2023).

## Features

- **Lag-profile inversion** of the alternating-code pulses: the
  autocorrelation function (ACF) at every 9 km gate by weighted least squares
  ([Virtanen et al. 2008](https://doi.org/10.5194/angeo-26-571-2008)), with
  ground-clutter subtraction and receiver calibration from the noise diode.
- **Satellite echoes as unknowns:** a single-pulse detector (matched filter,
  CFAR) finds the echoes, pass tracking links them, and each echo gets its own
  unknown in the inversion instead of being masked.
- **Power-line impulses:** detected in the voltage; the lagged products they
  touch are left out of the inversion, with leave-one-out weights that stay
  unbiased.
- **Fits** of $n_e$, $T_e$, $T_i$ and $v_i$ with full covariance, and quality
  flags for every gate.
- **Long-pulse** range–Doppler spectra and their fit.

## Quick start

```bash
git clone https://github.com/ionodev/isr_analysis.git
cd isr_analysis
bash tools/setup_env.sh                 # needs libfftw3-dev and an MPI library
cp config/example_basic.json config/my_run.json
# set data_dir, output_dir and the channel
mpirun -np 24 ~/venv/isr_analysis/bin/python3 run_analysis.py config/my_run.json
```

`config/example_v3.json` adds the satellite and impulse removal; it needs a
satellite catalogue and a wider theory table first
([user guide, section 5](docs/USER_GUIDE.md#5-satellite-and-impulse-removal-version-3)).
Tests: `python3 -m pytest tests`.

## Repository layout

`run_analysis.py` runs everything and is the only script in the top folder.

| Folder | Contents |
|---|---|
| `pipeline/` | the modules `run_analysis.py` uses |
| `calibration/` | the density scale from the plasma lines |
| `analysis/` | evaluation and figure scripts (not part of the pipeline; [README](analysis/README.md)) |
| `plotting/` | diagnostic plots |
| `tests/` | unit tests |
| `tools/` | environment setup, the wider theory table, the regression benchmark |
| `config/` | example configurations |
| `docs/` | the user guide, the configuration reference, the analysis-server guide |
| `syisr/`, `crates/` | Juha Vierinen's SYISR analyses and Rust range–Doppler estimator (below) |

## Documentation

- [User guide](docs/USER_GUIDE.md): installation, input data, running,
  output files, reading the results, the density scale
- [Configuration reference](docs/CONFIGURATION.md): every setting of every step
- [Analysis scripts](analysis/README.md)
- [Running on the analysis server](docs/appendix_analysis_guide.tex)

## SYISR and Rust range–Doppler processing

From [jvierine/isr_analysis](https://github.com/jvierine/isr_analysis):

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
