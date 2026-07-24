#!/usr/bin/env python
"""Dev diagnostic: SA-BCM FULL partial-derivative ladder -- Tapenade forward (_d) vs FD vs
complex-step, Tapenade reverse (_b) vs reverse-fast (_fast_b), and fwd/rev dot-product
(transpose) consistency -- all block-by-block (meanflow / nuTilde), for either SABCM_Exp
variant. Not the registered testflo suite (no ref files, no BaseRegTest) -- raw printed tables,
run directly. Mirrors the sibling SA-GR repo's 3-stage ladder (dot-product, _b-vs-_fast_b,
AD-vs-FD/CS) in dev-script form, generalized to SA-BCM's 2-block (not 4-block) state layout.

Five modes, matched to the installed build (ad/fd/bwd/dot need the REAL build; cs needs the
COMPLEX build -- they cannot run in the same process):

  --mode ad   REAL build. For each column block, seed a state perturbation masked to that
              block, compute dR/dw * wDot via Tapenade's forward-mode (_d) routines, print the
              row-block-masked norms (4 combos). Writes the table to --out (default: a JSON next
              to this script) so fd/cs mode can diff against it.

  --mode fd   REAL build, same sweep via finite difference (mode="FD", h=...). Prints its own
              table AND, if --ref resolves an ad-mode JSON, an abs/rel-error comparison.

  --mode cs   COMPLEX build, same sweep via complex step (h=1e-40). Same ref-comparison
              behavior as fd mode -- this is the decisive, step-free ground truth.

  --mode bwd  REAL build. For each ROW block, seed resBar masked to that block, compute wBar via
              the full reverse (_b) routines and wBarFast via the state-only reverse-fast
              (_fast_b) routines, print ||wBar - wBarFast|| per block (should be ~0 -- these are
              two independently Tapenade-generated code paths that must agree exactly).

  --mode dot  REAL build. For each (column, row) block pair, seed wDot on the column block and
              dwBar on the row block, and check the transpose identity
              dwBar . (dR/dw wDot) == wDot . (dR/dw^T dwBar) -- i.e. forward (_d) and reverse
              (_b) are consistent with each other. Proves consistency, not correctness (fwd/rev
              can both be wrong identically and still pass) -- cs mode is the correctness check.

State: linearizes about the CGNS written by dev/run_bcm_case.py's un-restarted run (not a formal
generate_bcm_restart.py restart -- this is a dev check, any state works for a Jacobian-vector-
product comparison, converged or not). Override with --restartfile if you have a better one.

examples:
  # real build: forward AD/FD, reverse consistency, dot products -- smooth variant
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_partials_bcm.py --variant smooth --mode ad
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_partials_bcm.py --variant smooth --mode fd
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_partials_bcm.py --variant smooth --mode bwd
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_partials_bcm.py --variant smooth --mode dot

  # complex build, decisive CS check against the ad-mode table
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_partials_bcm.py --variant smooth --mode cs

  # same five, hard variant
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_partials_bcm.py --variant hard --mode ad
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_partials_bcm.py --variant hard --mode fd
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_partials_bcm.py --variant hard --mode bwd
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_partials_bcm.py --variant hard --mode dot
  mpirun -np 2 --bind-to core /home/mdo/packages_v2/mach/bin/python dev/diag_partials_bcm.py --variant hard --mode cs
"""
import os
import sys
import copy
import json
import argparse

import numpy
from mpi4py import MPI

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from reg_default_options import adflowDefOpts
from reg_bcm import (
    ap_bcm_tut_wing,
    bcmBaseOptionsSmooth,
    bcmBaseOptionsHard,
    bcmAeroDVs,
    getStateBlocks,
    maskStateVector,
)

baseDir = os.path.dirname(os.path.abspath(__file__))
COMM = MPI.COMM_WORLD
RANK = COMM.rank


def rprint(*a):
    if RANK == 0:
        print(*a, flush=True)


def _default_restart(variant):
    return os.path.join(baseDir, "..", "output_files", "mdo_tutorial_bcm_%s_000_vol.cgns" % variant)


def _default_ref(variant):
    return os.path.join(baseDir, "partials_ad_bcm_%s.json" % variant)


def buildSolver(variant, complex_build, restartfile):
    if complex_build:
        from adflow import ADFLOW_C as ADF
    else:
        from adflow import ADFLOW as ADF

    variantOptions = bcmBaseOptionsSmooth if variant == "smooth" else bcmBaseOptionsHard

    options = copy.copy(adflowDefOpts)
    options["outputdirectory"] = os.path.join(baseDir, "../output_files")
    options.update(copy.deepcopy(variantOptions))
    options["restartfile"] = restartfile
    # only want the residual/Jacobian operator at this state, not a solve
    options["ncycles"] = 1
    options["usenksolver"] = False
    options["useanksolver"] = False

    solver = ADF(options=options, debug=True)

    ap = copy.deepcopy(ap_bcm_tut_wing)
    for dv in bcmAeroDVs:
        ap.addDV(dv)

    solver.getResidual(ap)  # propagate the restart state
    return solver, ap


def sweepBlocks(solver, seed, mode, h):
    """Column-by-column forward products dR/dw * e[block], row-block-masked. Returns
    {rowName: {colName: normValue}}."""
    nw, blocks = getStateBlocks(solver)
    wDotFull = solver.getStatePerturbation(seed)

    extraArgs = {} if mode is None else {"mode": mode, "h": h}

    table = {row: {} for row in blocks}
    for colName, colOffsets in blocks.items():
        wDot = maskStateVector(wDotFull, nw, colOffsets)
        if mode == "CS":
            # re-seat the real state between CS columns -- avoids imaginary-buffer
            # contamination from the prior column (see reg_bcm.py's identical note)
            solver.setStates(numpy.real(solver.getStates()))
        resDot = solver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True, **extraArgs)
        for rowName, rowOffsets in blocks.items():
            rowDot = maskStateVector(resDot, nw, rowOffsets)
            norm = numpy.sqrt(COMM.allreduce(numpy.sum(rowDot.astype(numpy.complex128).real ** 2) +
                                              numpy.sum(numpy.imag(rowDot) ** 2), op=MPI.SUM))
            table[rowName][colName] = float(norm)
    return table


def bwdConsistency(solver, seed):
    """Row-block-seeded ||wBar - wBarFast|| -- full reverse (_b) vs reverse-fast (_fast_b)."""
    nw, blocks = getStateBlocks(solver)
    dwBarFull = solver.getStatePerturbation(seed)

    results = {}
    for rowName, rowOffsets in blocks.items():
        dwBar = maskStateVector(dwBarFull, nw, rowOffsets)
        wBar = solver.computeJacobianVectorProductBwd(resBar=dwBar, wDeriv=True)
        wBarFast = solver.computeJacobianVectorProductBwdFast(resBar=dwBar)
        diff = wBar - wBarFast
        diffNorm = numpy.sqrt(COMM.allreduce(numpy.sum(diff ** 2), op=MPI.SUM))
        wBarNorm = numpy.sqrt(COMM.allreduce(numpy.sum(wBar ** 2), op=MPI.SUM))
        results[rowName] = (float(diffNorm), float(wBarNorm))
    return results


def dotProductCheck(solver, seed):
    """(column, row) block-paired transpose identity: dwBar.(dR/dw wDot) == wDot.(dR/dw^T dwBar)."""
    nw, blocks = getStateBlocks(solver)
    wDotFull = solver.getStatePerturbation(seed)
    dwBarFull = solver.getStatePerturbation(seed + 1)

    results = {}
    for colName, colOffsets in blocks.items():
        wDot = maskStateVector(wDotFull, nw, colOffsets)
        dwDot = solver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True)
        for rowName, rowOffsets in blocks.items():
            dwBar = maskStateVector(dwBarFull, nw, rowOffsets)
            wBar = solver.computeJacobianVectorProductBwd(resBar=dwBar, wDeriv=True)

            dotLocal1 = numpy.sum(dwDot * dwBar)
            dotLocal2 = numpy.sum(wDot * wBar)
            dot1 = COMM.allreduce(dotLocal1, op=MPI.SUM)
            dot2 = COMM.allreduce(dotLocal2, op=MPI.SUM)
            results[(colName, rowName)] = (float(dot1), float(dot2))
    return results


def print_bwd_table(title, results):
    rprint("\n==== %s: ||wBar - wBarFast|| by row block, _b vs _fast_b ====" % title)
    hdr = "%-12s %18s %18s %11s" % ("row block", "||wBar||", "||diff||", "rel diff")
    rprint(hdr)
    for row, (diffNorm, wBarNorm) in results.items():
        rel = diffNorm / max(wBarNorm, 1e-30)
        rprint("%-12s %18.6e %18.6e %11.2e" % (row, wBarNorm, diffNorm, rel))


def print_dot_table(title, results):
    rprint("\n==== %s: fwd/rev transpose consistency by (col, row) block pair ====" % title)
    hdr = "%-28s %18s %18s %11s" % ("w[col] -> R[row]", "dwBar.(dR/dw w)", "wDot.(dR/dw^T b)", "rel diff")
    rprint(hdr)
    for (col, row), (dot1, dot2) in results.items():
        rel = abs(dot1 - dot2) / max(abs(dot1), abs(dot2), 1e-30)
        rprint("%-28s %18.6e %18.6e %11.2e" % ("w[%s] -> R[%s]" % (col, row), dot1, dot2, rel))


def print_table(title, table):
    rprint("\n==== %s ====" % title)
    cols = list(next(iter(table.values())).keys())
    hdr = "%-12s" % "row \\ col" + "".join("%18s" % c for c in cols)
    rprint(hdr)
    for row, vals in table.items():
        rprint("%-12s" % row + "".join("%18.6e" % vals[c] for c in cols))


def print_comparison(title, ref, got):
    rprint("\n==== %s vs ref ====" % title)
    hdr = "%-28s %16s %16s %11s %11s" % ("dR[row]/dw[col]", "ref (AD)", "this run", "|abs err|", "rel err")
    rprint(hdr)
    rprint("-" * len(hdr))
    worst = (-1.0, "")
    for row, vals in got.items():
        for col, v in vals.items():
            r = ref.get(row, {}).get(col)
            label = "dR[%s]/dw[%s]" % (row, col)
            if r is None:
                rprint("%-28s %16s %16.6e %11s %11s" % (label, "(no ref)", v, "-", "-"))
                continue
            ad = abs(v - r)
            rd = ad / max(abs(r), 1e-30)
            rprint("%-28s %16.6e %16.6e %11.2e %11.2e" % (label, r, v, ad, rd))
            if rd > worst[0]:
                worst = (rd, label)
    if worst[0] >= 0:
        rprint("-" * len(hdr))
        rprint("worst relative error: %.3e at %s" % worst)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--variant", choices=["smooth", "hard"], required=True)
    p.add_argument("--mode", choices=["ad", "fd", "cs", "bwd", "dot"], required=True)
    p.add_argument("--restartfile", default=None)
    p.add_argument("--seed", type=int, default=314)
    p.add_argument("--fd-h", type=float, default=1e-7)
    p.add_argument("--out", default=None, help="ad mode: where to write the reference table (json)")
    p.add_argument("--ref", default=None, help="fd/cs mode: ad-mode table to compare against")
    args = p.parse_args()

    restartfile = args.restartfile or _default_restart(args.variant)
    if not os.path.isfile(restartfile):
        rprint("ERROR: restart file not found: %s" % restartfile)
        rprint("Run dev/run_bcm_case.py --variant %s first, or pass --restartfile." % args.variant)
        sys.exit(1)

    out = args.out or _default_ref(args.variant)
    ref_path = args.ref or _default_ref(args.variant)

    if args.mode == "ad":
        solver, ap = buildSolver(args.variant, complex_build=False, restartfile=restartfile)
        table = sweepBlocks(solver, args.seed, mode=None, h=None)
        print_table("AD (Tapenade forward _d), variant=%s" % args.variant, table)
        if RANK == 0:
            json.dump(table, open(out, "w"), indent=2)
            rprint("\nwrote AD reference table -> %s" % out)

    elif args.mode == "fd":
        solver, ap = buildSolver(args.variant, complex_build=False, restartfile=restartfile)
        table = sweepBlocks(solver, args.seed, mode="FD", h=args.fd_h)
        print_table("FD (h=%.0e), variant=%s" % (args.fd_h, args.variant), table)
        if os.path.isfile(ref_path):
            ref = json.load(open(ref_path))
            print_comparison("FD", ref, table)
        else:
            rprint("\n(no AD reference at %s -- run --mode ad first for a comparison)" % ref_path)

    elif args.mode == "cs":
        solver, ap = buildSolver(args.variant, complex_build=True, restartfile=restartfile)
        table = sweepBlocks(solver, args.seed, mode="CS", h=1e-40)
        print_table("CS (complex-step), variant=%s" % args.variant, table)
        if os.path.isfile(ref_path):
            ref = json.load(open(ref_path))
            print_comparison("CS", ref, table)
        else:
            rprint("\n(no AD reference at %s -- run --mode ad first for a comparison)" % ref_path)

    elif args.mode == "bwd":
        solver, ap = buildSolver(args.variant, complex_build=False, restartfile=restartfile)
        results = bwdConsistency(solver, args.seed)
        print_bwd_table("BWD vs BWDFast, variant=%s" % args.variant, results)

    else:  # dot
        solver, ap = buildSolver(args.variant, complex_build=False, restartfile=restartfile)
        results = dotProductCheck(solver, args.seed)
        print_dot_table("Dot products, variant=%s" % args.variant, results)


if __name__ == "__main__":
    main()
