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
# Usage:
#   ./run_sagr_tests.sh              run the whole suite (real + complex)
#   ./run_sagr_tests.sh real        real-build stages only (1, 2, AD/FD)
#   ./run_sagr_tests.sh cs          complex-build CS ground truth only
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

do_train() {
    head "Retraining JSON reference files (crossflow-converged state)"
    "$PY" -m testflo -n "$NP" "$FWD" "$BWD" -m "train*" -v
    echo "refs written: refs/jacvecfwd_sagr_flatplate.json  refs/jacvecbwd_sagr_flatplate.json"
}

do_genw() {
    head "Regenerating converged restart state (w) via dev/generate_sagr_restart.py"
    echo "NOTE: non-standard step (no download server) -- see dev/README.md"
    mpirun -np "$NP" --bind-to core "$PY" dev/generate_sagr_restart.py "${@:2}"
}

case "${1:-all}" in
    real)  run_real ;;
    cs)    run_cs ;;
    train) do_train ;;
    genw)  do_genw "$@" ;;
    all)
        run_real; rc_real=$?
        run_cs;   rc_cs=$?
        hr
        echo "SUMMARY"
        echo "  real build (Stage 1/2/3-AD-FD): $([ $rc_real -eq 0 ] && echo PASS || echo FAIL)"
        echo "  complex build (Stage 3 CS)    : $([ $rc_cs   -eq 0 ] && echo PASS || echo FAIL)"
        echo "  FD residual tests are @expectedFailure (metric noise on the"
        echo "  13-order residual); CS is the enforced ground truth."
        hr
        [ $rc_real -eq 0 ] && [ $rc_cs -eq 0 ]
        ;;
    *)
        echo "usage: $0 [all|real|cs|train|genw]" >&2; exit 2 ;;
esac
