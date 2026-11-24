#!/bin/bash
# Create a virtual environment for the isr_analysis ISR analysis pipeline.
# Run once on each machine. Activate before running: source ~/venv/isr_analysis/bin/activate
set -e

VENV_DIR="${HOME}/venv/isr_analysis"

python3 -m venv --system-site-packages "${VENV_DIR}"

"${VENV_DIR}/bin/pip" install --upgrade pip

# The pipeline's packages (requirements.txt); numpy, scipy, matplotlib and
# h5py come from the system when it has them. pyfftw needs libfftw3-dev, mpi4py
# an MPI implementation (e.g. OpenMPI).
"${VENV_DIR}/bin/pip" install -r "$(dirname "$0")/../requirements.txt"

echo ""
echo "Environment ready at ${VENV_DIR}"
echo "Activate with:  source ${VENV_DIR}/bin/activate"
echo "Run analysis:   mpirun -np 24 ${VENV_DIR}/bin/python3 run_analysis.py config/my_run.json"
