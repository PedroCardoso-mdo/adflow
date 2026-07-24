#!/usr/bin/env bash
#
# run_sagr_tests.sh -- one entry point for the SA-noft2-Gamma-Retheta (SA-GR)
# derivative regression suite.
#
# The SA-GR partials are validated on three levels (see
# docs/VERIFICATION/three-stage-verification.md):
#   Stage 1  dot-product consistency   forward _d  <->  reverse _b
#   Stage 2  fast-reverse consistency  reverse _b  <->  reverse-fast _fast_b
#   Stage 3  ground truth              forward AD  vs  complex-step (CS) / FD
#
# Stages 1-2 and the AD/FD half of Stage 3 run on the REAL build; the decisive
# CS half of Stage 3 runs on the COMPLEX build (ADFLOW_C). This script drives
# both and prints a compact per-stage summary.
#
# On top of the partials (dR/dw, dR/dXv Jacobian-vector products) the suite
# also validates the FULL TOTAL derivatives df/dx -- the complete adjoint,
# exactly like the SA test_adjoint.py: TestAdjointSAGR solves the SA-GR
# adjoint (reverse _b, 8-state) for df/d{alpha,mach,twist,span,shape} and
# TestCmplxStepSAGR re-converges the complex build and checks it by CS. This
# is the "adjoint" stage below (test_adjoint_sagr.py).
#
# On top of the derivatives, the suite also checks the PRIMAL residual operator
# itself: test_blockette_sagr.py asserts the cache-blocked "blockette" residual
# (blocketteResCore) equals the reference "block" residual (saGammaReTheta_block)
# for the same state w, across all 8 variables -- guarding the inlined SA-GR
# kernels in blockette.F90 against drift (this is the "blockette" stage below).
#
# Usage:
#   ./run_sagr_tests.sh              run the whole suite (real + complex)
#   ./run_sagr_tests.sh real        real-build stages only (1, 2, AD/FD)
#   ./run_sagr_tests.sh cs          complex-build CS ground truth only
#   ./run_sagr_tests.sh adjoint     full total-derivative adjoint vs CS
#   ./run_sagr_tests.sh blockette   blockette residual == block residual
#   ./run_sagr_tests.sh train       regenerate the JSON reference files
#   ./run_sagr_tests.sh genw        regenerate the converged restart state (w)
#
# Everything the case depends on -- mesh, restart (w), AeroProblem, options,
# crossflow on/off -- lives in reg_sagr.py. To move to a different mesh: point
# sagrGridFile/sagrRestartFile there at the new CGNS, run `genw` if you need a
# fresh converged state, then `train` to rebuild the JSON, then run the suite.
#
set -uo pipefail

cd "$(dirname "$(readlink -f "$0")")"

# --- knobs (override from the environment) --------------------------------
PY=${PY:-/home/mdo/packages_v2/mach/bin/python}
NP=${NP:-2}          # MPI ranks for the real-build stages
NP_CS=${NP_CS:-1}    # MPI ranks for the complex CS stage
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}   # never oversubscribe (see CLAUDE memory)

FWD=test_jacVecProdFWD_sagr.py
BWD=test_jacVecProdBWDFast_sagr.py
ADJ=test_adjoint_sagr.py
BLK=test_blockette_sagr.py

hr()  { printf '%.0s-' {1..72}; echo; }
head() { hr; echo ">>> $*"; hr; }

run_real() {
    head "Real build -- Stage 1 (dot products), Stage 2 (_b vs _fast_b), Stage 3 AD/FD"
    "$PY" -m testflo -n "$NP" "$FWD" "$BWD" -v
}

run_cs() {
    head "Complex build -- Stage 3 CS ground truth (16 coupling blocks, wDot/xVDot/xDvDot)"
    "$PY" -m testflo -n "$NP_CS" "$FWD" -m "cmplx_test_*" -v
}

run_adjoint() {
    head "Full total derivatives df/dx -- adjoint (real) then CS check (complex)"
    "$PY" -m testflo -n "$NP" "$ADJ" -v
    "$PY" -m testflo -n "$NP" "$ADJ" -m "cmplx_test_*" -v
}

run_blockette() {
    head "Blockette residual == block residual (SA-GR, same w, all 8 vars)"
    "$PY" -m testflo -n "$NP" "$BLK" -v
}

do_train() {
    head "Retraining JSON reference files (crossflow-converged state)"
    "$PY" -m testflo -n "$NP" "$FWD" "$BWD" "$ADJ" -m "train*" -v
    echo "refs written: refs/jacvecfwd_sagr_tut_wing.json  refs/jacvecbwd_sagr_tut_wing.json  refs/adjoint_sagr_tut_wing.json"
}

do_genw() {
    head "Regenerating converged restart state (w) via dev/generate_sagr_restart.py"
    echo "NOTE: non-standard step (no download server) -- see dev/README.md"
    mpirun -np "$NP" --bind-to core "$PY" dev/generate_sagr_restart.py "${@:2}"
}

case "${1:-all}" in
    real)      run_real ;;
    cs)        run_cs ;;
    adjoint)   run_adjoint ;;
    blockette) run_blockette ;;
    train)     do_train ;;
    genw)      do_genw "$@" ;;
    all)
        run_real;     rc_real=$?
        run_cs;       rc_cs=$?
        run_adjoint;  rc_adj=$?
        run_blockette; rc_blk=$?
        hr
        echo "SUMMARY"
        echo "  real build (Stage 1/2/3-AD-FD): $([ $rc_real -eq 0 ] && echo PASS || echo FAIL)"
        echo "  complex build (Stage 3 CS)    : $([ $rc_cs   -eq 0 ] && echo PASS || echo FAIL)"
        echo "  full adjoint df/dx (real + CS): $([ $rc_adj  -eq 0 ] && echo PASS || echo FAIL)"
        echo "  blockette == block residual   : $([ $rc_blk  -eq 0 ] && echo PASS || echo FAIL)"
        echo "  FD residual tests are @expectedFailure (metric noise on the"
        echo "  13-order residual); CS is the enforced ground truth."
        hr
        [ $rc_real -eq 0 ] && [ $rc_cs -eq 0 ] && [ $rc_adj -eq 0 ] && [ $rc_blk -eq 0 ]
        ;;
    *)
        echo "usage: $0 [all|real|cs|adjoint|blockette|train|genw]" >&2; exit 2 ;;
esac
