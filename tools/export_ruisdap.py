#!/usr/bin/env python3
"""
Export one day of a finished run into the RUISDAP Millstone viewer's format
(pipeline/ruisdap_export.py; run_analysis.py's step "ruisdap_export" does the
same at the end of a run).

Usage:
    python3 tools/export_ruisdap.py RUN OUTPUT --date 2024-04-08 \\
        --plasma-manifest BASE/manifest.json --plasma-url /guisdap/data/BASE/plasma.f32 \\
        --label "..." [--magic-constant C] [--postfix lpi_60] [--channel zenith-l] \\
        [--echo-catalogue CATALOGUE | --echo-rows ROWS.npy] [--no-quality-flags]
"""
import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pipeline"))
import ruisdap_export as rex  # noqa: E402

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('run', type=Path, help="the run's output_dir (contains <postfix>/<channel>)")
p.add_argument('output', type=Path, help='the new export folder')
p.add_argument('--date', required=True, help='the UTC day, YYYY-MM-DD')
p.add_argument('--plasma-manifest', type=Path, required=True, help="manifest.json of the day's RUISDAP export")
p.add_argument('--plasma-url', required=True, help="URL of that export's plasma spectrum file")
p.add_argument('--label', required=True, help='what the run is, shown in the viewer')
p.add_argument('--magic-constant', type=float, default=None, help="density constant (default: fit_quality's of the day)")
p.add_argument('--calibration-method', default=None, help='how the constant was found, for the manifest')
p.add_argument('--postfix', default='lpi_60')
p.add_argument('--channel', default='zenith-l')
p.add_argument('--antenna', default='ZENITH')
p.add_argument('--source-commit', default='', help='isr_analysis commit of the run')
p.add_argument('--echo-catalogue', default=None, help='satellite catalogue for the SATELLITE flag')
p.add_argument('--echo-rows', type=Path, default=None, help='instead: a .npy of its rows (period, range km, SNR dB)')
p.add_argument('--no-quality-flags', action='store_true')
a = p.parse_args()
rex.export_day(a.run, a.output, a.date, a.plasma_manifest, a.plasma_url, a.label, postfix=a.postfix,
               channel=a.channel, antenna=a.antenna, magic_constant=a.magic_constant,
               calibration_method=a.calibration_method, source_commit=a.source_commit,
               quality=not a.no_quality_flags, echo_catalogue=a.echo_catalogue, echo_rows=a.echo_rows)
