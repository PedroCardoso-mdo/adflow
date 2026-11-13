#!/bin/bash
#SBATCH --job-name=transition_examples
#SBATCH --account=CHANGE_ME
#SBATCH --partition=normal-x86
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --signal=USR2@300
#SBATCH --no-requeue
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Generic Slurm template for the transition examples: ONE node of 128 cores, several cases side by
# side, one srun per case. Written for Deucalion (normal-x86: 128 cores, 242 GB, billed per node
# whatever is requested), so the node is filled: 4 x 32 ranks below.
#
#   mkdir -p logs && sbatch job_template.sh
#
# Rules built into the template
#   * every case is its own subshell with its own timeout: a hung case cannot block the others;
#   * the timeout sends SIGUSR2 first: ADflow writes the solution and leaves the solve, the script
#     still writes result.json; SIGKILL follows 300 s later;
#   * to stop a job by hand:  scancel --signal=USR2 <jobid>   (a plain scancel loses the flow state);
#   * OMP_NUM_THREADS=1 and contiguous cores (-m block:block --cpu-bind=cores);
#   * ADflow is memory-bound: for 3D cases give more ranks to one case rather than packing more
#     cases (sickle wing: 64 ranks and ~100 GB; AR5 L1: 64 ranks and ~105 GB; AR5 L0 and the c376
#     optimisation need 2 and 4 nodes).
#   * every log starts with the ADflow banner (version | git hash | install time): check it.
set -uo pipefail

# ---- environment: edit these three lines -------------------------------------------------------
source "$HOME/env_adflow.sh"                       # modules + python environment with this ADflow build
EX=$SLURM_SUBMIT_DIR                               # examples/transition of the repository
export TRANSITION_MESH_ROOT=$EX/../../showcase/meshes
# -------------------------------------------------------------------------------------------------
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
J=${SLURM_JOB_ID:-nojob}
cd "$EX"; mkdir -p logs runs
python -c "import adflow; from adflow import gitInfo as g; print('adflow:', adflow.__file__, g.GIT_HASH, g.GIT_BRANCH, 'DIRTY' if g.GIT_DIRTY else 'clean')"

# lane <tag> <ranks> <time limit> <script> <arguments...>
lane() {
  local tag=$1 n=$2 limit=$3; shift 3
  (
    echo "START $tag n=$n limit=$limit $(date)"
    timeout -s USR2 -k 300 "$limit" \
      srun --exact -n "$n" -c 1 -m block:block --cpu-bind=cores \
      python "$@" --outputDir "runs/$tag" > "logs/${tag}_$J.log" 2>&1
    echo "END   $tag rc=$? $(date) :: $(grep -h '^RESULT' "logs/${tag}_$J.log" | tail -1)"
  ) &
  sleep 5      # sruns start in order, so the lanes occupy consecutive core blocks
}

# ---- 4 x 32 ranks = 128: the three models and plain SA on the NLF(1)-0416 at alpha 4 ------------
lane nlf0416_gr_a4     32 3h run_airfoil.py --airfoil nlf0416 --model gr     --alpha 4 --writeVolume
lane nlf0416_sgamma_a4 32 3h run_airfoil.py --airfoil nlf0416 --model sgamma --alpha 4 --nuinf 2.79e-7
lane nlf0416_bcmo_a4   32 3h run_airfoil.py --airfoil nlf0416 --model bcmo   --alpha 4
lane nlf0416_sa_a4     32 3h run_airfoil.py --airfoil nlf0416 --model sa     --alpha 4

# Other fillings of one node (replace the block above):
#   swept wing, 7 Reynolds numbers (7 x 18 = 126):
#     for re in 1.796238e6 2.0e6 2.203762e6 2.37e6 2.498432e6 3.0e6 3.498432e6; do
#       lane swept_re$re 18 3h run_swept_nlf0415.py --reynolds $re; done
#   sickle wing + AR5 L2 (64 + 64):
#     lane sickle 64 23h run_sickle_wing.py --writeVolume
#     lane ar5_L2 64 3h  run_ar5_wing.py --model gr --level L2 --writeVolume
#   two single-point optimisations (64 + 64, 48 h):
#     lane c2_gr     64 47h opt_airfoil.py --model gr --points 0.3:0.5
#     lane c2_gr_mtx 64 47h opt_airfoil.py --model gr --points 0.3:0.5 --dissipation matrix
wait
echo "job $J finished $(date)"
