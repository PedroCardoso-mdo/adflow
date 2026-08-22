#!/usr/bin/env bash
set -euo pipefail

AMDUPROF="/opt/AMDuProf_5.2-606/bin/AMDuProfPcm"
MACH_VENV="/home/mdo/packages_v2/mach/bin/activate"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PYTHON_SCRIPT="${SCRIPT_DIR}/analysis.py"
BLOCKETTE_FILE="${REPO_ROOT}/src/NKSolver/blockette.F90"
RESULTS_ROOT="${SCRIPT_DIR}/results_bs_sweep"
PYTHON_BIN=""

cd "${SCRIPT_DIR}"

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

ORIGINAL_BS=""

get_current_bs() {
    sed -n -E 's/^[[:space:]]*integer\(kind=intType\), parameter :: BS = ([0-9]+)/\1/p' "${BLOCKETTE_FILE}" | head -n1
}

set_bs() {
    local bs="$1"
    sed -i -E "s/^([[:space:]]*integer\(kind=intType\), parameter :: BS = )[0-9]+/\1${bs}/" "${BLOCKETTE_FILE}"

    local updated_bs
    updated_bs="$(get_current_bs)"
    if [[ "${updated_bs}" != "${bs}" ]]; then
        echo "Failed to set BS=${bs} in ${BLOCKETTE_FILE}" >&2
        exit 1
    fi

    echo "Set BS=${bs} in ${BLOCKETTE_FILE}"
}

restore_original_bs() {
    if [[ -n "${ORIGINAL_BS}" ]]; then
        set_bs "${ORIGINAL_BS}"
        echo "Restored BS=${ORIGINAL_BS}"
    fi
}

activate_env_and_install() {
    # shellcheck disable=SC1090
    source "${MACH_VENV}"

    if [[ -z "${VIRTUAL_ENV:-}" ]]; then
        echo "Virtual environment is not active after sourcing ${MACH_VENV}" >&2
        exit 1
    fi

    PYTHON_BIN="${VIRTUAL_ENV}/bin/python"
    if [[ ! -x "${PYTHON_BIN}" ]]; then
        echo "Python executable not found in virtual environment: ${PYTHON_BIN}" >&2
        exit 1
    fi

    pushd "${REPO_ROOT}" >/dev/null
    make
    "${PYTHON_BIN}" -m pip install .
    popd >/dev/null
}

run_case() {
    local base_dir="$1"
    local np="$2"
    local omp="$3"

    local case_dir="${base_dir}/np_${np}_omp_${omp}"
    mkdir -p "${case_dir}"

    echo "Running ${case_dir}"

    export OMP_NUM_THREADS="${omp}"
    export OMP_PROC_BIND=true
    export OMP_PLACES=cores

    ${AMDUPROF} \
        -a -A system -C \
        -o "${case_dir}/mem_np${np}_omp${omp}.csv" \
        -- \
        mpirun --bind-to core --report-bindings -np "${np}" "${PYTHON_BIN}" "${PYTHON_SCRIPT}" \
        > "${case_dir}/output_np${np}_omp${omp}.log" 2>&1
}

run_suite() {
    local suite_root="$1"

    mkdir -p "${suite_root}"
    echo "=== 1A: 1 to 16 ranks, 1 OpenMP thread ==="
    run_pure_mpi "${suite_root}/results_smt_off_pure_mpi"
}

run_pure_mpi() {
    local base_dir="$1"
    for np in $(seq 1 16); do
        run_case "${base_dir}" "${np}" 1
    done
}

run_hybrid() {
    local base_dir="$1"
    for np in $(seq 1 16); do
        for omp in $(seq 2 15); do
            if (( np * omp < 16 )); then
                run_case "${base_dir}" "${np}" "${omp}"
            fi
        done
    done
}

offline_smt() {
    echo "Turning OFF logical CPUs 16-31"
    for i in /sys/devices/system/cpu/cpu{16..31}/online; do
        echo 0 | sudo tee "$i" >/dev/null
    done
}

online_smt() {
    echo "Turning ON logical CPUs 16-31"
    for i in /sys/devices/system/cpu/cpu{16..31}/online; do
        echo 1 | sudo tee "$i" >/dev/null
    done
}

mkdir -p "${RESULTS_ROOT}"

ORIGINAL_BS="$(get_current_bs)"
if [[ -z "${ORIGINAL_BS}" ]]; then
    echo "Could not determine original BS value from ${BLOCKETTE_FILE}" >&2
    exit 1
fi

trap restore_original_bs EXIT

echo "Using original BS=${ORIGINAL_BS}"

echo "Running baseline case with BS=${ORIGINAL_BS}"
activate_env_and_install
run_suite "${RESULTS_ROOT}/bs_${ORIGINAL_BS}"

for bs in 2 4 10 12 16; do
    echo "Running sweep case with BS=${bs}"
    set_bs "${bs}"
    activate_env_and_install
    run_suite "${RESULTS_ROOT}/bs_${bs}"
done



