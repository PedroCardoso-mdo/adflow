#!/usr/bin/env bash
#
# run_sgamma_tests.sh -- one entry point for the SA-noft2-Gamma (SA-sgamma)
# derivative regression suite. Same ladder as run_sagr_tests.sh, see
# docs/VERIFICATION/VERIF_00_three_stage_verification.md:
#   Stage 1  dot-product consistency   forward _d  <->  reverse _b
#   Stage 2  fast-reverse consistency  reverse _b  <->  reverse-fast _fast_b
#   Stage 3  ground truth              forward AD  vs  complex-step (CS) / FD
# plus the full total derivatives df/dx (adjoint vs CS, test_adjoint_sgamma.py).
# Stages 1-2 and AD/FD run on the REAL build, CS on the COMPLEX build.
# There is no SA-sgamma blockette test (blockettes are forced off for this model).
#
# Usage:
#   ./run_sgamma_tests.sh           run the whole suite (real + complex)
#   ./run_sgamma_tests.sh real      real-build stages only (1, 2, AD/FD)
#   ./run_sgamma_tests.sh cs        complex-build CS ground truth only
#   ./run_sgamma_tests.sh adjoint   full total-derivative adjoint vs CS
#   ./run_sgamma_tests.sh train     regenerate the JSON reference files
#   ./run_sgamma_tests.sh genw      regenerate the converged restart state (w)
#
# Everything the case depends on -- mesh, restart (w), AeroProblem, options --
# lives in reg_sgamma.py.
#
set -uo pipefail

cd "$(dirname "$(readlink -f "$0")")"

# --- knobs (override from the environment) --------------------------------
PY=${PY:-/home/mdo/packages_v2/mach/bin/python}
NP=${NP:-2}          # MPI ranks for the real-build stages
NP_CS=${NP_CS:-1}    # MPI ranks for the complex CS stage
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}   # never oversubscribe (see CLAUDE memory)

FWD=test_jacVecProdFWD_sgamma.py
BWD=test_jacVecProdBWDFast_sgamma.py
ADJ=test_adjoint_sgamma.py

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

do_train() {
    head "Retraining JSON reference files"
    "$PY" -m testflo -n "$NP" "$FWD" "$BWD" "$ADJ" -m "train*" -v
    echo "refs written: refs/jacvecfwd_sgamma_tut_wing.json  refs/jacvecbwd_sgamma_tut_wing.json  refs/adjoint_sgamma_tut_wing.json"
}

do_genw() {
    head "Regenerating converged restart state (w) via dev/generate_sgamma_restart.py"
    echo "NOTE: non-standard step (no download server) -- see dev/README.md"
    mpirun -np "$NP" --bind-to core "$PY" dev/generate_sgamma_restart.py "${@:2}"
}

case "${1:-all}" in
    real)      run_real ;;
    cs)        run_cs ;;
    adjoint)   run_adjoint ;;
    train)     do_train ;;
    genw)      do_genw "$@" ;;
    all)
        run_real;     rc_real=$?
        run_cs;       rc_cs=$?
        run_adjoint;  rc_adj=$?
        hr
        echo "SUMMARY"
        echo "  real build (Stage 1/2/3-AD-FD): $([ $rc_real -eq 0 ] && echo PASS || echo FAIL)"
        echo "  complex build (Stage 3 CS)    : $([ $rc_cs   -eq 0 ] && echo PASS || echo FAIL)"
        echo "  full adjoint df/dx (real + CS): $([ $rc_adj  -eq 0 ] && echo PASS || echo FAIL)"
        echo "  FD residual tests are @expectedFailure (metric noise on the"
        echo "  13-order residual); CS is the enforced ground truth."
        hr
        [ $rc_real -eq 0 ] && [ $rc_cs -eq 0 ] && [ $rc_adj -eq 0 ]
        ;;
    *)
        echo "usage: $0 [all|real|cs|adjoint|train|genw]" >&2; exit 2 ;;
esac
