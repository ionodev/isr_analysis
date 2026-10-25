#!/usr/bin/env python3
"""
Driver for Millstone Hill ISR analysis pipeline.

Usage:
    python run_analysis.py <config.json>
    mpirun -np 8 python run_analysis.py <config.json>

Steps in the pipeline:
    lpi         - Lag-profile inversion (coded long pulse, ACF estimation),
                  optionally with the catalogued satellite echoes as unknowns
    fit_lpi     - Fit ACFs to plasma parameters
    long_pulse  - Range-Doppler spectra (uncoded long pulse)
    fit_lp      - Fit Doppler spectra to plasma parameters
"""

import sys
import os
import json
import numpy as n

# Ensure local modules (stuffr, il_interp, millstone_radar_state, …) are found
# when MPI workers are spawned in a different working directory.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Fall back to Agg only when no display is available (headless / MPI on server).
# If DISPLAY or WAYLAND_DISPLAY is set, use whatever backend matplotlib picks normally.
if not os.environ.get("DISPLAY") and not os.environ.get("WAYLAND_DISPLAY") and not os.environ.get("MPLBACKEND"):
    os.environ["MPLBACKEND"] = "Agg"

# Many ranks read the same interpolation tables and lpi files at once. HDF5
# locks even read-only opens, which makes those collide. Must be set before
# h5py loads the HDF5 library, so keep this above any import that pulls in h5py.
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")


def barrier():
    """
    Wait until every MPI rank has finished the current step. Each rank does
    its own share of a step and then goes on to the next, so without this a
    rank starts fit_lpi (or fit_lp) while others are still writing the files
    it groups: the fit windows then depend on timing, and some are fitted
    from incomplete groups (seen on 12 October in the variant runs). Without
    mpi4py, or with one process, it does nothing.
    """
    try:
        from mpi4py import MPI
    except ImportError:
        return
    MPI.COMM_WORLD.Barrier()


def load_config(path):
    with open(path) as f:
        return json.load(f)


def satellite_detections(datadir, cfg):
    """
    The satellite stage of the lpi step: the catalogued echoes of the
    channel, for outlier_lpi to model as extra unknowns, one column per echo.
    None when the configuration has no enabled "satellites" section.

        "satellites": {
            "catalogue": path of the single pulse CFAR catalogue (Digital
                         Metadata, single_pulse_satellite_cfar.py),
            "channel_delay_us": tx-h to echo channel delay; null measures it
                         from the data (tx_delay.py),
            "min_snr_db": only echoes at least this strong (optional),
            "track": link the detections into passes and give a column to
                     every pulse of a pass, detected or not
                     (satellite_columns.TrackedDetections; default false),
            "track_options": keyword arguments of TrackedDetections
                     (optional; {"min_weak_pulses": 3} confirms the weak
                     detections of a lower-threshold catalogue by tracking),
            "strong_pfa": with a catalogue made at a lower false-alarm
                     probability, the probability whose threshold marks a
                     detection strong (e.g. 1e-10, the default catalogue's);
                     the others are weak (optional),
            "random_columns": keyword arguments of
                     satellite_columns.RandomControlColumns, e.g.
                     {"per_period": 200, "seed": 1}: extra columns at random
                     places without an echo, the pure-plasma control
                     (optional; for tests only),
            "enabled": true (optional)
        }
    """
    sat = cfg.get("satellites")
    if not sat or not sat.get("enabled", True):
        return None
    import satellite_columns as satcol
    delay = sat.get("channel_delay_us")
    if delay is None:
        import tx_delay
        delay = tx_delay.estimate_channel_delay(datadir, cfg["channel"], verbose=False)[0]
        if delay is None:
            raise RuntimeError("the channel delay of %s could not be measured; give satellites.channel_delay_us"
                               % cfg["channel"])
    print("satellite stage: catalogue %s, channel delay %.2f us" % (sat["catalogue"], delay))
    det = satcol.CatalogueDetections(sat["catalogue"], cfg["channel"], delay,
                                     min_snr_db=sat.get("min_snr_db"), strong_pfa=sat.get("strong_pfa"))
    if sat.get("track", False):
        det = satcol.TrackedDetections(det, satcol.zenith_pulse_keys(datadir), **sat.get("track_options", {}))
    if sat.get("random_columns"):
        # the pure-plasma control: columns where no echo is
        det = satcol.RandomControlColumns(det, satcol.zenith_pulse_keys(datadir), **sat["random_columns"])
    return det


def run_lpi(datadir, output_base, cfg, max_time_s=None):
    import outlier_lpi as olpi
    olpi.lpi_files(
        dirname=datadir,
        channel=cfg["channel"],
        rg=cfg["range_gate_us"],
        avg_dur=cfg.get("avg_dur", 10),
        min_tx_frac=cfg.get("min_tx_frac", 0.5),
        pass_band=cfg["pass_band_hz"],
        filter_len=cfg.get("filter_len", 20),
        maximum_range_delay=cfg.get("max_range_delay_us", 7000),
        save_acf_images=cfg.get("save_acf_images", True),
        lag_avg=cfg.get("lag_avg", 1),
        reanalyze=cfg.get("reanalyze", False),
        output_base=output_base,
        max_time_s=max_time_s,
        # the LPI variant (defaults: variance weights, background column,
        # outlier rejection)
        weighting=cfg.get("weighting", "variance"),
        background_column=cfg.get("background_column", True),
        outlier_rejection=cfg.get("outlier_rejection", True),
        # the satellite stage, and the periods to analyse (None for all)
        satellite_detections=satellite_detections(datadir, cfg),
        periods=cfg.get("periods", None),
    )


def run_fit_lpi(datadir, output_base, cfg, lpi_cfg, radar_freq_hz=440.2e6, table_dir=None):
    import fit_lpi as flpi
    postfix = "_%d" % lpi_cfg["range_gate_us"]
    flpi.fit_lpifiles(
        dirn=datadir,
        channel=cfg["channel"],
        postfix=postfix,
        max_dt=cfg.get("max_dt", 300),
        plot=False,
        first_lag=cfg.get("first_lag", 0),
        reanalyze=cfg.get("reanalyze", False),
        range_avg=n.array(cfg.get("range_avg", [1, 3, 5])),
        output_base=output_base,
        radar_freq_hz=radar_freq_hz,
        table_dir=table_dir,
        # null (off), "absolute" or "significance": fit_lpi.space_object_flagged
        space_object_cut=cfg.get("space_object_cut", None),
    )


def run_long_pulse(datadir, output_base, cfg, max_time_s=None):
    import avg_range_doppler_spec as ards
    ards.avg_range_doppler_spectra(
        dirname=datadir,
        channel=cfg["channel"],
        mode=cfg.get("mode", 300),
        avg_dur=cfg.get("avg_dur", 10),
        step=cfg.get("step", 10),
        avg_type=cfg.get("avg_type", "outlier_removal"),
        postfix=cfg.get("postfix", "_outlier"),
        min_tx_pulses=cfg.get("min_tx_pulses", 100),
        reanalyze=cfg.get("reanalyze", False),
        output_base=output_base,
        max_time_s=max_time_s,
        tx_delay_us=cfg.get("tx_delay_us", None),
    )


def run_fit_lp(datadir, output_base, cfg, lp_cfg, radar_freq_hz=440.2e6, table_dir=None):
    import fit_lp as flp
    postfix = "_%d%s" % (lp_cfg.get("mode", 300), lp_cfg.get("postfix", "_outlier"))
    flp.fit_spectra(
        dirname=datadir,
        channel=cfg["channel"],
        postfix=postfix,
        avg_dur=cfg.get("avg_dur", 600),
        ridx=cfg.get("ridx", [35, 230]),
        remove_space_objects=cfg.get("remove_space_objects", False),
        reanalyze=cfg.get("reanalyze", False),
        output_base=output_base,
        radar_freq_hz=radar_freq_hz,
        table_dir=table_dir,
        fit_bandwidth_hz=cfg.get("fit_bandwidth_hz", 50e3),
        notch_bands_hz=cfg.get("notch_bands_hz", None),
        pulse_length_us=cfg.get("pulse_length_us", None),
    )


def main():
    if len(sys.argv) < 2:
        print("Usage: python run_analysis.py <config.json>")
        sys.exit(1)

    config = load_config(sys.argv[1])
    datadir = config["data_dir"]
    output_base = config.get("output_dir", None)
    max_time_s = config.get("max_time_s", None)
    radar_freq_hz = config.get("radar_freq_hz", 440.2e6)
    table_dir = config.get("table_dir", os.path.join(os.path.dirname(os.path.abspath(__file__)), "data"))
    steps = config.get("steps", {})

    lpi_cfg = steps.get("lpi", {})
    lp_cfg = steps.get("long_pulse", {})

    if lpi_cfg.get("enabled", False):
        print("=== Step: LPI (lag-profile inversion) ===")
        run_lpi(datadir, output_base, lpi_cfg, max_time_s=max_time_s)
        barrier()

    if steps.get("fit_lpi", {}).get("enabled", False):
        if not lpi_cfg:
            print("ERROR: fit_lpi requires lpi config for postfix")
            sys.exit(1)
        print("=== Step: fit_lpi (ACF fitting) ===")
        run_fit_lpi(datadir, output_base, steps["fit_lpi"], lpi_cfg,
                    radar_freq_hz=radar_freq_hz, table_dir=table_dir)

    if lp_cfg.get("enabled", False):
        print("=== Step: long_pulse (range-Doppler spectra) ===")
        run_long_pulse(datadir, output_base, lp_cfg, max_time_s=max_time_s)
        barrier()

    if steps.get("fit_lp", {}).get("enabled", False):
        if not lp_cfg:
            print("ERROR: fit_lp requires long_pulse config for postfix")
            sys.exit(1)
        print("=== Step: fit_lp (Doppler spectrum fitting) ===")
        run_fit_lp(datadir, output_base, steps["fit_lp"], lp_cfg,
                   radar_freq_hz=radar_freq_hz, table_dir=table_dir)


if __name__ == "__main__":
    main()
