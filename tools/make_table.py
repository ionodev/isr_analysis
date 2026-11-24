#!/usr/bin/env python3
"""
The wider ion-line table of pipeline version 3 (config/example_v3.json).

The default table (isr_spec.il_table, built on the first run) has Te/Ti 1-3
in steps of 0.4, Ti 100-4000 K in steps of 205 K and the ion fraction in
steps of 0.11. The fit interpolates linearly between these nodes, so fitted
values collect at them, and many topside fits end at the edge Te/Ti = 3.
This writes, into <out>/wide, a table with
    Te/Ti 1-5 in steps of 0.2, Ti 100-4000 K in steps of 100 K, ion fraction
    0-1 in steps of 0.05.
The densities are 1e11 and 1e12 m^-3 only: fit_lpi evaluates the table at
1e12 (Debye-length effects are ignored in the fit), so no other density node
is ever read. Each file also stores the ACFs (il_interp computes them on
first load otherwise, and many ranks would write them at once).

Usage: mpirun -np N python3 tools/make_table.py <out dir>
then set "table_dir": "<out dir>/wide" in the configuration.
"""
import os
import sys

import numpy as np
from mpi4py import MPI

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pipeline"))
import isr_spec  # noqa: E402
import il_interp  # noqa: E402

out = os.path.join(sys.argv[1], "wide")
comm = MPI.COMM_WORLD
os.makedirs(out, exist_ok=True)
grid = dict(nes=np.array([1e11, 1e12]), te_ti_ratios=np.linspace(1, 5, 21),
            tis=np.linspace(100, 4000, 40), frs=np.linspace(0, 1, 21))
for m0, m1 in ((32, 16), (16, 1)):
    isr_spec.il_table(mass0=m0, mass1=m1, radar_freq=440.2e6, outdir=out, **grid)
    comm.Barrier()
    if comm.Get_rank() == 0:
        # load once: computes and stores the ACFs ("A")
        il_interp.ilint(radar_freq=440.2e6, ion_mass1=m0, ion_mass2=m1, table_dir=out, verbose=False)
        print("wrote ion_line_interpolate_%d_%d_440.2.h5" % (m0, m1), flush=True)
    comm.Barrier()
