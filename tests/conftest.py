import os
import sys

# The tests import the pipeline modules (pipeline/), the calibration
# (calibration/) and run_analysis.py (the top folder) by their plain names.
TOP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for d in ("calibration", "pipeline"):
    sys.path.insert(0, os.path.join(TOP, d))
sys.path.insert(0, TOP)
