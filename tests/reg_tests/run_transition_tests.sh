#!/usr/bin/env bash
#
# run_transition_tests.sh -- convenience driver for the transition-model regression tests.
# It only wraps the testflo commands listed in README_transition.md.
#
# Usage:
#   ./run_transition_tests.sh <model> <stage>
#
#   model: sagr | sgamma | sabcm
#   stage:
#     real       fixed-state tests on the real build: forward and reverse products vs the refs,
#                dot products, reverse vs reverse-fast, AD vs FD (and blockette vs block for sagr)
#     cs         fixed-state tests on the complex build: forward products vs complex step
#     train      retrain the refs of the fixed-state tests (real build)
#     adjoint    total derivatives: adjoint vs ref (real build). Solves the adjoint: use a cluster.
#     adjoint-cs total derivatives: complex step vs adjoint ref (complex build). Converges the flow
#                once per design variable: use a cluster.
#     train-adjoint  retrain the adjoint ref (real build). Solves the adjoint: use a cluster.
#     genw       regenerate the converged restart state with dev/generate_<model>_restart.py.
#                Converges the flow: use a cluster. Extra arguments are passed to the script.
#
# Environment:
#   PY      python interpreter (default: python)
#   NP      MPI ranks for genw (default: 2). The tests set their own rank count (N_PROCS).
#   NTEST   number of tests testflo runs concurrently (default: 1)
#
# sabcm has no adjoint test and no restart state (see README_transition.md).
#
set -uo pipefail

cd "$(dirname "$(readlink -f "$0")")"

PY=${PY:-python}
NP=${NP:-2}
NTEST=${NTEST:-1}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

usage() {
    echo "usage: $0 <sagr|sgamma|sabcm> <real|cs|train|adjoint|adjoint-cs|train-adjoint|genw> [genw args]" >&2
    exit 2
}

[ $# -ge 2 ] || usage
MODEL=$1
STAGE=$2

case "$MODEL" in
    sagr | sgamma | sabcm) ;;
    *) usage ;;
esac

FWD=test_jacVecProdFWD_$MODEL.py
BWD=test_jacVecProdBWDFast_$MODEL.py
ADJ=test_adjoint_$MODEL.py
FIXED="$FWD $BWD"
[ "$MODEL" = sagr ] && FIXED="$FIXED test_blockette_sagr.py"

flo() { "$PY" -m testflo -n "$NTEST" -v "$@"; }

needs_solve() {
    if [ "$MODEL" = sabcm ]; then
        echo "stage '$STAGE' is not available for sabcm" >&2
        exit 2
    fi
}

case "$STAGE" in
    real) flo $FIXED ;;
    cs) flo "$FWD" -m "cmplx_test_*" ;;
    train) flo "$FWD" "$BWD" -m "train*" ;;
    adjoint)
        needs_solve
        flo "$ADJ"
        ;;
    adjoint-cs)
        needs_solve
        flo "$ADJ" -m "cmplx_test_*"
        ;;
    train-adjoint)
        needs_solve
        flo "$ADJ" -m "train*"
        ;;
    genw)
        needs_solve
        mpirun -np "$NP" "$PY" "dev/generate_${MODEL}_restart.py" "${@:3}"
        ;;
    *) usage ;;
esac
