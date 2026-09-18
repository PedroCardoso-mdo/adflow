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
# Real stages under env_mdoV4.sh; complex stages under env_mdo_cs.sh (libcomplexify + complex PETSc
# on the path — the real library must NOT see that path, hence separate sub-shells). The pg variant
# refs are trained first (BCM_VARIANTS=pg), smooth/hard refs untouched. Tiny 2-rank tests: this job
# cannot fill a node (8 tasks asked).
set -uo pipefail
BASE=/projects/F202500002HPCVLABISTUL/pedrocardoso/MDO_Lab_IST
T=$BASE/MDOLab/adflow_sa-bcm/tests/reg_tests
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
J=$SLURM_JOB_ID; cd $T; mkdir -p logs; echo "head $(git rev-parse --short HEAD) $(date)"
real () { ( source $BASE/env_mdoV4.sh; cd $T; python -m testflo -n 2 "$@" -v ); }
cplx () { ( source $BASE/env_mdo_cs.sh; cd $T; python -m testflo -n 1 "$@" -v ); }
summ () { grep -E "^(Passed|Failed|Skipped)" $1 | tr '\n' ' '; echo; }
if [[ "${TRAIN_PG:-1}" == 1 && ! -f refs/jacvecfwd_bcm_pg_tut_wing.json ]]; then
  echo "=== train pg refs $(date)"; BCM_VARIANTS=pg real test_jacVecProdFWD_bcm.py test_jacVecProdBWDFast_bcm.py test_adjoint_bcm.py -m "train*" > logs/regtests_trainpg_$J.log 2>&1; summ logs/regtests_trainpg_$J.log; ls refs | grep pg
fi
echo "=== blockette $(date)"; real test_blockette_bcm.py > logs/regtests_blockette_$J.log 2>&1; summ logs/regtests_blockette_$J.log
echo "=== real $(date)";      real test_jacVecProdFWD_bcm.py test_jacVecProdBWDFast_bcm.py > logs/regtests_real_$J.log 2>&1; summ logs/regtests_real_$J.log
echo "=== cs $(date)";        cplx test_jacVecProdFWD_bcm.py -m "cmplx_test_*" > logs/regtests_cs_$J.log 2>&1; summ logs/regtests_cs_$J.log
echo "=== adjoint real $(date)"; real test_adjoint_bcm.py > logs/regtests_adjoint_$J.log 2>&1; summ logs/regtests_adjoint_$J.log
echo "=== adjoint cs $(date)";   cplx test_adjoint_bcm.py -m "cmplx_test_*" > logs/regtests_adjointcs_$J.log 2>&1; summ logs/regtests_adjointcs_$J.log
for f in logs/regtests_*_$J.log; do echo "--- $f"; grep -E "\.\.\. (FAIL|OK)" $f | sed -E 's/\(mpi [0-9]\) //; s/ \(00:.*//' | awk '{print $NF, $1}' | sort | uniq -c | awk '{print $2, $1}' | sort | uniq -c | awk '{printf "%s x%s  ", $2, $1} END {print ""}'; grep -E "\.\.\. FAIL" $f | sed -E 's/\(mpi [0-9]\) //; s/ +\.\.\..*//'; done
date
