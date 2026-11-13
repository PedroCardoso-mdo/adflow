#!/bin/bash
# Two-iteration smoke test of every transition example on a workstation.
#
# It proves that script + mesh + options + build work together; it does NOT converge anything.
# Each case runs alone (never in parallel), from free stream, and stops after two solver
# iterations (--nCycles caps "Iter Tot": 4 is enough for the first two ANK steps, 5 for the
# SA-BCM wing which starts at CFL 1). The optimisation scripts run with --baselineOnly:
# the problem is built and the start design evaluated; no optimiser iteration, no adjoint.
#
# Output: <SMOKE_DIR>/<case>/{run.log,result.json,...} and <SMOKE_DIR>/SMOKE_REPORT.md, which
# compares iteration 0-2 of every case with the cluster run of the same case
# (smoke_reference.json holds the first rows of the cluster logs).
#
# Usage:   ./smoke_all.sh                 # every case
#          ./smoke_all.sh sickle ar5      # only the cases whose name contains one of the words
# Environment (defaults in brackets):
#          NP [6]          ranks per case
#          SMOKE_DIR       [<repo>/showcase/smoke_runs]
#          MPIRUN          ["mpirun -np $NP --bind-to core"]
#          SKIP_LARGE=1    skip the meshes above 3 M cells (sickle wing 4.4 M: ~60 GB; AR5 L1 3.7 M)
# The c376 optimisation mesh (13.7 M cells) never runs here; its set-up is exercised on the L4
# mesh of the same family.
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
NP=${NP:-6}
SMOKE_DIR=${SMOKE_DIR:-$REPO/showcase/smoke_runs}
MPIRUN=${MPIRUN:-"mpirun -np $NP --bind-to core"}
MESH=${TRANSITION_MESH_ROOT:-$REPO/showcase/meshes}
export OMP_NUM_THREADS=1
cd "$HERE"
mkdir -p "$SMOKE_DIR"

python -c "import adflow; from adflow import gitInfo as g; print('adflow:', adflow.__file__, g.GIT_HASH, g.GIT_BRANCH, 'dirty' if g.GIT_DIRTY else 'clean')" || exit 1

# run <case> <script> <args...>
run() {
  local name=$1; shift
  if [ $# -gt 0 ] && [ -n "${FILTER:-}" ]; then
    local keep=0; for w in $FILTER; do [[ "$name" == *"$w"* ]] && keep=1; done
    [ $keep = 1 ] || return 0
  fi
  local out=$SMOKE_DIR/$name
  rm -rf "$out"; mkdir -p "$out"
  echo -n "$(date +%H:%M:%S)  $name ... "
  local t0=$SECONDS
  $MPIRUN python "$@" --outputDir "$out" > "$out/run.log" 2>&1
  local rc=$?
  echo "$rc" > "$out/exit_code"
  echo "rc=$rc  $((SECONDS - t0)) s  $(grep -c -E '^ +1 +[0-9]+ +[0-9]+ +' "$out/run.log") iteration rows"
}
FILTER="$*"
N=4       # two ANK iterations
A="--nCycles $N"

# ---- run_airfoil.py: NACA 0012 (three models, mesh-study point; article point A)
run airfoil_naca0012_gr_scalar          run_airfoil.py --airfoil naca0012 --model gr     --alpha 2.7 $A
run airfoil_naca0012_sgamma_scalar      run_airfoil.py --airfoil naca0012 --model sgamma --alpha 2.7 $A
run airfoil_naca0012_bcm_scalar         run_airfoil.py --airfoil naca0012 --model bcm    --alpha 2.7 $A
run airfoil_naca0012_gr_scalar_pointA   run_airfoil.py --airfoil naca0012 --model gr     --tu 0.3 --alpha 4.52 $A
run airfoil_naca0012_bcmo_scalar_pointA run_airfoil.py --airfoil naca0012 --model bcmo   --tu 0.3 --alpha 4.52 $A
# ---- run_airfoil.py: NLF0416 and S809 validation points, scalar and matrix
run airfoil_nlf0416_gr_scalar           run_airfoil.py --airfoil nlf0416 --model gr     --alpha 4 $A
run airfoil_nlf0416_gr_matrix           run_airfoil.py --airfoil nlf0416 --model gr     --dissipation matrix --alpha 4 $A
run airfoil_nlf0416_sgamma_scalar       run_airfoil.py --airfoil nlf0416 --model sgamma --nuinf 2.79e-7 --alpha 4 $A
run airfoil_nlf0416_sgamma_matrix       run_airfoil.py --airfoil nlf0416 --model sgamma --nuinf 2.79e-7 --dissipation matrix --alpha 4 $A
run airfoil_nlf0416_bcm_scalar          run_airfoil.py --airfoil nlf0416 --model bcm    --alpha 4 $A
run airfoil_nlf0416_bcmo_scalar         run_airfoil.py --airfoil nlf0416 --model bcmo   --alpha 4 $A
run airfoil_nlf0416_bcmo_matrix         run_airfoil.py --airfoil nlf0416 --model bcmo   --dissipation matrix --alpha 4 $A
run airfoil_nlf0416_sa_scalar           run_airfoil.py --airfoil nlf0416 --model sa     --alpha 4 $A
run airfoil_s809_gr_scalar              run_airfoil.py --airfoil s809    --model gr     --alpha 4 2 $A      # two alphas: the polar loop
run airfoil_s809_gr_matrix              run_airfoil.py --airfoil s809    --model gr     --dissipation matrix --alpha 4 $A
# ---- run_swept_nlf0415.py: crossflow, matrix (validated) and scalar (counter-example, as it ran: Re on sqrt(2))
run swept_nlf0415_matrix                run_swept_nlf0415.py --reynolds 2.0e6 $A
run swept_nlf0415_scalar                run_swept_nlf0415.py --reynolds 2.0e6 --dissipation scalar --reynoldsLength 1.41421356 $A
# ---- run_ar5_wing.py
run ar5_gr_L4                           run_ar5_wing.py --model gr --level L4 $A
run ar5_bcm                             run_ar5_wing.py --model bcm --nCycles 5
run ar5_sa                              run_ar5_wing.py --model sa --alpha 1.5 --nCycles 5
# ---- opt_airfoil.py: build the problem and evaluate the NACA 0012 at point A (Tu 0.3 %, cl 0.5)
B="--points 0.3:0.5 --baselineOnly $A"
run opt_airfoil_gr_scalar               opt_airfoil.py --model gr     $B
run opt_airfoil_gr_matrix               opt_airfoil.py --model gr     --dissipation matrix $B
run opt_airfoil_sgamma_matrix           opt_airfoil.py --model sgamma --dissipation matrix --nuinf 2.79e-7 --opt sgammaArtVisc=5.0 $B
run opt_airfoil_bcmo_matrix             opt_airfoil.py --model bcmo   --dissipation matrix $B
run opt_airfoil_sa_scalar               opt_airfoil.py --model sa     $B
run opt_airfoil_gr_matrix_multipoint    opt_airfoil.py --model gr     --dissipation matrix --points 0.3:0.5,0.6:0.5,0.3:0.7 --baselineOnly $A
# ---- opt_ar5_wing.py
run opt_ar5_bcm_article                 opt_ar5_wing.py --case bcm_article --baselineOnly --nCycles 5
run opt_ar5_gr_c376_on_L4               opt_ar5_wing.py --case gr_c376 --baselineOnly $A --mesh "$MESH/ar5_wing/corrected_foil_ar5_vol_L4.cgns"
# ---- large meshes last
if [ "${SKIP_LARGE:-0}" != 1 ]; then
  run ar5_gr_L1                         run_ar5_wing.py --model gr --level L1 $A
  run sickle_wing                       run_sickle_wing.py $A
fi

python smoke_compare.py "$SMOKE_DIR"
