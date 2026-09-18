#!/bin/bash
#SBATCH --job-name=bcm_regtests
#SBATCH --account=f202500002hpcvlabistulx
#SBATCH --partition=normal-x86
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=logs/regtests_%j.out
# SA-BCM derivative regression suite on Deucalion (machV4 real + complex builds), 2026-09-18.
# Stages: blockette (block == blockette residual), real (dot products fwd<->rev, _b<->_fast_b, AD/FD),
# cs (AD vs complex step), adjoint (total derivatives vs refs / CS). All three variants
# (smooth, hard, pg) run per file. Tiny 2-rank tests: this job cannot fill a node (8 tasks asked).
# Refs for the pg variant do not exist yet -> those RegTest cases fail on "missing ref" until
# `./run_bcm_tests.sh train` is run; the self-contained FD/CS/dot-product classes are the point.
set -uo pipefail
BASE=/projects/F202500002HPCVLABISTUL/pedrocardoso/MDO_Lab_IST
source $BASE/env_mdoV4.sh
cd $BASE/MDOLab/adflow_sa-bcm/tests/reg_tests; mkdir -p logs
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PY=python NP=2 NP_CS=1
echo "head $(git rev-parse --short HEAD) $(date)"
for mode in ${MODES:-blockette real cs adjoint}; do
  echo "=== $mode $(date)"; timeout 100m ./run_bcm_tests.sh $mode > logs/regtests_${mode}_$SLURM_JOB_ID.log 2>&1; echo "rc=$?"
  grep -E "^(OK|FAIL|Passed|Failed|Skipped|Ran|.*tests? (passed|failed))" logs/regtests_${mode}_$SLURM_JOB_ID.log | tail -5
done
date
