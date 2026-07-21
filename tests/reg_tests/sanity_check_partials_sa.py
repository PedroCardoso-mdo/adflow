#!/usr/bin/env python
"""Sanity check (stages 2/3 of the plain-SA partial-derivative isolation
ladder, see docs/current-task.md): "my script" (reg_sagr.py's block-splitting
helpers, generalized to the plain-SA nw=6 layout) against either the
tutorial-wing SA mesh (stage 2, --mesh tutorial) or the AR5 mesh (stage 3,
--mesh ar5), both plain SA.

Re-verifies, at a single fixed restart state (getResidual() once, no
reconvergence), the two checks that already exist and are already
documented for SA-GR:
  - dot-product (fwd/rev transpose) consistency -- reg_test_utils'
    assert_dot_products_allclose (generic, nw-agnostic, unmodified)
  - reverse-fast (_fast_b) vs full-reverse (_b) consistency, block-seeded --
    reg_sagr.assert_bwdfast_blocks_allclose (now generalized to accept
    nw=6/plain-SA in getStateBlocks)

This is a quick sanity check for these two rungs -- it isn't the thing
under investigation (that's the 3-way FD/CS/AD forward check, built
separately). --mesh ar5 previously found a real w->R dot-product mismatch
(rel err ~1.4e-5 vs tol 2e-10) -- --mesh tutorial isolates whether that's
mesh-specific or a script bug.

Usage: mpirun -n 2 /home/mdo/packages_v2/mach/bin/python sanity_check_partials_sa.py --mesh {tutorial,ar5}
"""
import argparse
import os
import copy

from mpi4py import MPI
from adflow import ADFLOW

import reg_test_utils as utils
import reg_sagr
from reg_default_options import adflowDefOpts

parser = argparse.ArgumentParser()
parser.add_argument("--mesh", choices=["tutorial", "ar5"], default="ar5")
args = parser.parse_args()

baseDir = os.path.dirname(os.path.abspath(__file__))
rank = MPI.COMM_WORLD.rank

options = copy.copy(adflowDefOpts)
options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])

if args.mesh == "ar5":
    from test_adjoint_ar5_sa import ar5SAOptions, ap_ar5_sa

    options.update(copy.deepcopy(ar5SAOptions))
    ap = copy.deepcopy(ap_ar5_sa)
else:
    # same rans_tut_wing fixture as test_functionals.py's stage-1 case --
    # known-good SA tutorial-wing mesh+restart. Options copied verbatim from
    # test_functionals.py's rans_tut_wing test_params entry (not a trimmed
    # subset) to rule out any option difference as the cause of a mismatch.
    from reg_aeroproblems import ap_tutorial_wing
    from reg_default_options import defaultAeroDVs

    options.update(
        {
            "gridfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
            "restartfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
            "mgcycle": "sg",
            "equationtype": "RANS",
            "smoother": "DADI",
            "cfl": 1.5,
            "cflcoarse": 1.25,
            "resaveraging": "never",
            "nsubiter": 3,
            "nsubiterturb": 3,
            "ncyclescoarse": 100,
            "ncycles": 1000,
            "monitorvariables": ["cpu", "resrho", "resturb", "cl", "cd", "cmz", "yplus", "totalr"],
            "usenksolver": True,
            "l2convergence": 1e-14,
            "l2convergencecoarse": 1e-4,
            "nkswitchtol": 1e-3,
            "adjointl2convergence": 1e-14,
            "frozenturbulence": False,
        }
    )
    ap = copy.deepcopy(ap_tutorial_wing)
    for dv in defaultAeroDVs:
        ap.addDV(dv)

if rank == 0:
    print("\n=== mesh: %s ===\n" % args.mesh, flush=True)

CFDSolver = ADFLOW(options=options, debug=True)

# propagate restart state through the code once -- no reconvergence
CFDSolver.getResidual(ap)


class _PrintHandler:
    """Minimal stand-in for BaseRegTest that actually asserts (not just
    prints) and does the same MPI reduction BaseRegTest.par_add_sum does
    (baseclasses/testing/pyRegTest.py:222, comm.reduce(np.sum(values)) on
    rank 0) -- these are PARTIAL per-rank sums (state arrays are
    domain-decomposed), so comparing them without reducing across ranks
    first is meaningless and was the cause of an earlier false-alarm
    "mismatch" here."""

    def __init__(self, comm):
        self._stored = {}
        self.comm = comm

    def root_print(self, msg):
        if rank == 0:
            print(msg, flush=True)

    def par_add_norm(self, key, vec, **kwargs):
        import numpy

        reducedSum = self.comm.reduce(numpy.sum(vec**2))
        if rank == 0:
            n = numpy.sqrt(reducedSum)
            print("  norm(%s) = %.6e" % (key, n), flush=True)

    def par_add_sum(self, key, val, compare=False, rtol=1e-10, atol=1e-10, **kwargs):
        import numpy

        reducedSum = self.comm.reduce(numpy.sum(val))
        if rank != 0:
            return
        if not compare:
            self._stored[key] = reducedSum
            print("  %s (fwd) = %.6e" % (key, reducedSum), flush=True)
        else:
            ref = self._stored[key]
            print("  %s (rev) = %.6e" % (key, reducedSum), flush=True)
            numpy.testing.assert_allclose(reducedSum, ref, rtol=rtol, atol=atol, err_msg=key)


handler = _PrintHandler(MPI.COMM_WORLD)

if rank == 0:
    print("\n=== dot-product consistency (generic, plain SA) ===\n", flush=True)
utils.assert_dot_products_allclose(handler, CFDSolver, tol=2e-10)

if rank == 0:
    print("\n=== reverse-fast (_fast_b) vs full-reverse (_b), block-seeded ===\n", flush=True)
reg_sagr.assert_bwdfast_blocks_allclose(CFDSolver, seed=314, atol=1e-16)

if rank == 0:
    print("\n=== SANITY CHECK PASSED (mesh=%s) ===\n" % args.mesh, flush=True)
