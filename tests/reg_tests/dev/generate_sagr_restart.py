"""Generate the restart state of the SA-noft2-Gamma-Retheta regression tests.

Converges the tutorial wing from freestream with the options and flow conditions of
reg_transition.py and writes the double-precision volume solution (grid included) that the tests
read as restart file:

    input_files/mdo_tutorial_sagr_dp.cgns

This is a full flow solve on the real build; run it in parallel, for example

    mpirun -np 16 python generate_sagr_restart.py

and retrain the reference files afterwards (see ../README_transition.md).
"""

from restart_utils import generateRestart

if __name__ == "__main__":
    generateRestart("sagr", __doc__)
