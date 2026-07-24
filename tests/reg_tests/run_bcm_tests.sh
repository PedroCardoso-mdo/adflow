#!/usr/bin/env bash
#
# run_bcm_tests.sh -- one entry point for the SA-BCM derivative regression suite.
# Mirrors the sibling SA-GR repo's run_sagr_tests.sh 1:1.
#
# The SA-BCM partials are validated on the same three levels as SA-GR:
#   Stage 1  dot-product consistency   forward _d  <->  reverse _b
#   Stage 2  fast-reverse consistency  reverse _b  <->  reverse-fast _fast_b
#   Stage 3  ground truth              forward AD  vs  complex-step (CS) / FD
#
# Stages 1-2 and the AD/FD half of Stage 3 run on the REAL build; the decisive
# CS half of Stage 3 runs on the COMPLEX build (ADFLOW_C).
#
# On top of the partials, the suite also validates the FULL TOTAL derivatives
# df/dx via test_adjoint_bcm.py (the "adjoint" stage), and the PRIMAL residual
# operator itself via test_blockette_bcm.py (the "blockette" stage) -- see that
# file's docstring for why this one is safety-critical for SA-BCM specifically
# (useBlockettes is NOT force-disabled here, unlike SA-GR).
#
# Both SA-BCM variants (SABCM_Exp False/True, "smooth"/"hard") are parameterized
# into every test file via parameterized_class, so a single testflo invocation
# per file already covers both -- no separate flag needed.
#
# Usage:
#   ./run_bcm_tests.sh              run the whole suite (real + complex)
#   ./run_bcm_tests.sh real        real-build stages only (1, 2, AD/FD)
#   ./run_bcm_tests.sh cs          complex-build CS ground truth only
#   ./run_bcm_tests.sh adjoint     full total-derivative adjoint vs CS
#   ./run_bcm_tests.sh blockette   blockette residual == block residual (both variants)
#   ./run_bcm_tests.sh train       regenerate the JSON reference files
#   ./run_bcm_tests.sh genw        regenerate the converged restart states (both variants)
#
# Everything the case depends on -- mesh, restarts (w), AeroProblem, options --
# lives in reg_bcm.py. To move to a different mesh: point bcmGridFile/
# bcmRestartFileSmooth/Hard there at the new CGNS, run `genw` for a fresh
# converged state, then `train` to rebuild the JSON, then run the suite.
#
set -uo pipefail

cd "$(dirname "$(readlink -f "$0")")"

# --- knobs (override from the environment) --------------------------------
PY=${PY:-/home/mdo/packages_v2/mach/bin/python}
NP=${NP:-2}          # MPI ranks for the real-build stages
NP_CS=${NP_CS:-1}    # MPI ranks for the complex CS stage
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}   # never oversubscribe

FWD=test_jacVecProdFWD_bcm.py
BWD=test_jacVecProdBWDFast_bcm.py
ADJ=test_adjoint_bcm.py
BLK=test_blockette_bcm.py

hr()  { printf '%.0s-' {1..72}; echo; }
head() { hr; echo ">>> $*"; hr; }

run_real() {
    head "Real build -- Stage 1 (dot products), Stage 2 (_b vs _fast_b), Stage 3 AD/FD (smooth+hard)"
    "$PY" -m testflo -n "$NP" "$FWD" "$BWD" -v
}

run_cs() {
    head "Complex build -- Stage 3 CS ground truth (smooth+hard)"
    "$PY" -m testflo -n "$NP_CS" "$FWD" -m "cmplx_test_*" -v
}

run_adjoint() {
    head "Full total derivatives df/dx -- adjoint (real) then CS check (complex), smooth+hard"
    "$PY" -m testflo -n "$NP" "$ADJ" -v
    "$PY" -m testflo -n "$NP" "$ADJ" -m "cmplx_test_*" -v
}

run_blockette() {
    head "Blockette residual == block residual (SA-BCM, same w, all 6 vars, smooth+hard)"
    "$PY" -m testflo -n "$NP" "$BLK" -v
}

do_train() {
    head "Retraining JSON reference files (smooth+hard)"
    "$PY" -m testflo -n "$NP" "$FWD" "$BWD" "$ADJ" -m "train*" -v
    echo "refs written: refs/jacvecfwd_bcm_{smooth,hard}_tut_wing.json  refs/jacvecbwd_bcm_{smooth,hard}_tut_wing.json  refs/adjoint_bcm_{smooth,hard}_tut_wing.json"
}

do_genw() {
    head "Regenerating converged restart states (w) via generate_bcm_restart.py, both variants"
    echo "NOTE: expect this NOT to converge cleanly first try -- use dev/run_bcm_case.py to iterate first"
    mpirun -np "$NP" --bind-to core "$PY" generate_bcm_restart.py --variant smooth "${@:2}"
    mpirun -np "$NP" --bind-to core "$PY" generate_bcm_restart.py --variant hard "${@:2}"
}

case "${1:-all}" in
    real)      run_real ;;
    cs)        run_cs ;;
    adjoint)   run_adjoint ;;
    blockette) run_blockette ;;
    train)     do_train ;;
    genw)      do_genw "$@" ;;
    all)
        run_blockette; rc_blk=$?
        run_real;     rc_real=$?
        run_cs;       rc_cs=$?
        run_adjoint;  rc_adj=$?
        hr
        echo "SUMMARY"
        echo "  blockette == block residual   : $([ $rc_blk  -eq 0 ] && echo PASS || echo FAIL)  (run FIRST -- see test_blockette_bcm.py docstring)"
        echo "  real build (Stage 1/2/3-AD-FD): $([ $rc_real -eq 0 ] && echo PASS || echo FAIL)"
        echo "  complex build (Stage 3 CS)    : $([ $rc_cs   -eq 0 ] && echo PASS || echo FAIL)"
        echo "  full adjoint df/dx (real + CS): $([ $rc_adj  -eq 0 ] && echo PASS || echo FAIL)"
        echo "  FD residual tests may be noisy near the tanh-blend/KS-max kinks; CS is the"
        echo "  enforced ground truth."
        hr
        [ $rc_blk -eq 0 ] && [ $rc_real -eq 0 ] && [ $rc_cs -eq 0 ] && [ $rc_adj -eq 0 ]
        ;;
    *)
        echo "usage: $0 [all|real|cs|adjoint|blockette|train|genw]" >&2; exit 2 ;;
esac
