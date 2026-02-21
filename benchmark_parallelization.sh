#!/usr/bin/env bash
################################################################################
# 3DPolyS-LE Parallelization Methods Benchmark Script
################################################################################
#
# This script benchmarks different parallelization methods for 3DPolyS-LE
# simulations to compare performance and verify reproducibility.
#
# Usage:
#   ./benchmark_parallelization.sh [options] [config_file] [seed] [output_dir]
#
# Options:
#   --tests=TYPE         Select which tests to run:
#                        all      - Run all tests (default)
#                        mpi      - Run both MPI tests (with GPU and without GPU/OpenACC)
#                        no-mpi   - Run only non-MPI tests (single process: with GPU or without GPU)
#                        single   - Same as no-mpi: single-process (no MPI) tests
#                        gpu      - Run "with GPU/OpenACC" tests only (MPI and no-MPI)
#                        cpu      - Run "without GPU (--no-gpu)" tests only (MPI and no-MPI)
#                        gpu-vs-no-gpu - Run only the two single-process tests (with vs without GPU)
#
#   --max-process=NUM    Maximum number of MPI processes to use (default: all available CPUs)
#
# Arguments:
#   config_file  - Path to simulation config file (default: ./test/test_tads_init_benchmark.cfg)
#   seed         - Random seed for reproducibility (default: 42)
#   output_dir   - Output directory for results (default: ./benchmark_results)
#
# Examples:
#   ./benchmark_parallelization.sh
#   ./benchmark_parallelization.sh --tests=mpi          # both MPI tests (with GPU, without GPU)
#   ./benchmark_parallelization.sh --tests=no-mpi       # single process only (with GPU, without GPU)
#   ./benchmark_parallelization.sh --max-process=8
#   ./benchmark_parallelization.sh --tests=gpu-vs-no-gpu # only single-process GPU vs --no-gpu comparison
#
################################################################################

set -e  # Exit on error

# Show help and exit (used by -h/--help)
show_help() {
    cat << 'HELP'
3DPolyS-LE Parallelization Methods Benchmark Script
====================================================

Benchmarks different parallelization methods for 3DPolyS-LE simulations
to compare performance and verify reproducibility.

Usage:
  ./benchmark_parallelization.sh [options] [config_file] [seed] [output_dir]

Options:
  -h, --help             Show this help and exit (no run).

  --tests=TYPE           Select which tests to run:
                         all      - Run all tests (default)
                         mpi      - Run both MPI tests (with GPU, without GPU/OpenACC)
                         no-mpi   - Run only non-MPI tests (single process)
                         single   - Same as no-mpi
                         gpu      - Run "with GPU/OpenACC" tests only
                         cpu      - Run "without GPU (--no-gpu)" tests only
                         gpu-vs-no-gpu - Run only single-process with GPU vs without GPU

  --max-process=NUM      Maximum number of MPI processes to use (default: all available CPUs)

Arguments:
  config_file  - Path to simulation config file (default: ./test/test_tads_init_benchmark.cfg)
  seed         - Random seed for reproducibility (default: 42)
  output_dir   - Output directory for results (default: ./benchmark_results)

Examples:
  ./benchmark_parallelization.sh -h
  ./benchmark_parallelization.sh --tests=mpi
  ./benchmark_parallelization.sh --tests=no-mpi
  ./benchmark_parallelization.sh --max-process=8
  ./benchmark_parallelization.sh --tests=gpu --max-process=16 ./test/test_tads_init_benchmark.cfg 42 ./my_benchmarks

Test plan (when --tests=all; N = --max-process or number of CPUs):
  MPI (multi-process) tests:
    [1] MPI N procs, with GPU/OpenACC
    [2] MPI N procs, without GPU/OpenACC (--no-gpu, OpenMP per process)
  No-MPI (single process) tests:
    [3] With GPU/OpenACC, no MPI
    [4] Without GPU (--no-gpu), no MPI — uses OpenMP for replica parallelism
  Reproducibility: with the same seed, all four runs must produce identical
  output (GPU, OpenMP, MPI+GPU, MPI+no-GPU). Checksums verify this.
HELP
}

# Handle -h / --help (anywhere in args)
for arg in "$@"; do
    if [[ "$arg" == "-h" || "$arg" == "--help" ]]; then
        show_help
        exit 0
    fi
done

# Parse options
TEST_SELECTION="all"
MAX_PROCESS=""
while [[ "$1" =~ ^-- ]]; do
    case "$1" in
        --tests=*)
            TEST_SELECTION="${1#*=}"
            shift
            ;;
        --max-process=*)
            MAX_PROCESS="${1#*=}"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Default parameters
CONFIG_FILE="${1:-./test/test_tads_init_benchmark.cfg}"
SEED="${2:-42}"
OUTPUT_DIR="${3:-./benchmark_results}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
BENCH_DIR="${OUTPUT_DIR}/${TIMESTAMP}"
LOG_DIR="${BENCH_DIR}/logs"
RESULTS_DIR="${BENCH_DIR}/results"

# Create directories
mkdir -p "${LOG_DIR}"
mkdir -p "${RESULTS_DIR}"

# Executable path
EXECUTABLE="./bin/3dpolys_le"

# Check if executable exists
if [ ! -f "${EXECUTABLE}" ]; then
    echo "ERROR: Executable not found at ${EXECUTABLE}"
    echo "Please build the project first: cmake --build cmake-build"
    exit 1
fi

# Check if config file exists
if [ ! -f "${CONFIG_FILE}" ]; then
    echo "ERROR: Config file not found: ${CONFIG_FILE}"
    exit 1
fi

# Detect hardware configuration
detect_hardware() {
    echo "=============================================================================="
    echo "HARDWARE CONFIGURATION"
    echo "=============================================================================="

    # System information
    echo "System Information:"
    echo "  Hostname:        $(hostname)"
    echo "  OS:              $(uname -s) $(uname -r)"
    echo "  Architecture:    $(uname -m)"
    echo "  Date:            $(date)"
    echo ""

    # CPU information
    echo "CPU Information:"
    if [[ "$(uname -s)" == "Darwin" ]]; then
        # macOS
        CPU_MODEL=$(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo "Unknown")
        CPU_CORES=$(sysctl -n hw.physicalcpu 2>/dev/null || echo "Unknown")
        CPU_LOGICAL=$(sysctl -n hw.logicalcpu 2>/dev/null || echo "Unknown")
        echo "  Model:           ${CPU_MODEL}"
        echo "  Physical Cores:  ${CPU_CORES}"
        echo "  Logical Cores:   ${CPU_LOGICAL}"
    else
        # Linux
        CPU_MODEL=$(grep "model name" /proc/cpuinfo | head -1 | cut -d: -f2 | xargs)
        CPU_CORES=$(grep "^cpu cores" /proc/cpuinfo | head -1 | cut -d: -f2 | xargs)
        CPU_LOGICAL=$(nproc)
        echo "  Model:           ${CPU_MODEL}"
        echo "  Physical Cores:  ${CPU_CORES}"
        echo "  Logical Cores:   ${CPU_LOGICAL}"
    fi
    echo ""

    # Memory information
    echo "Memory Information:"
    if [[ "$(uname -s)" == "Darwin" ]]; then
        MEM_TOTAL=$(sysctl -n hw.memsize | awk '{printf "%.1f GB", $1/1024/1024/1024}')
        echo "  Total Memory:    ${MEM_TOTAL}"
    else
        MEM_TOTAL=$(free -h | grep "Mem:" | awk '{print $2}')
        echo "  Total Memory:    ${MEM_TOTAL}"
    fi
    echo ""

    # GPU information
    echo "GPU Information:"
    if command -v nvidia-smi &> /dev/null; then
        GPU_INFO=$(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "Unable to query")
        echo "  NVIDIA GPU:      ${GPU_INFO}"
    elif [[ "$(uname -s)" == "Darwin" ]]; then
        GPU_INFO=$(system_profiler SPDisplaysDataType 2>/dev/null | grep "Chipset Model:" | cut -d: -f2 | xargs || echo "Unknown")
        echo "  GPU:             ${GPU_INFO}"
    else
        echo "  GPU:             No NVIDIA GPU detected"
    fi
    echo ""

    # MPI configuration
    echo "MPI Configuration:"
    if command -v mpirun &> /dev/null; then
        MPI_VERSION=$(mpirun --version 2>&1 | head -1 || echo "Unknown")
        echo "  MPI Available:   Yes"
        echo "  Version:         ${MPI_VERSION}"
    else
        echo "  MPI Available:   No"
    fi
    echo ""

    echo "=============================================================================="
    echo ""
}

# Detect available CPUs
if [[ "$(uname -s)" == "Darwin" ]]; then
    AVAILABLE_CPUS=$(sysctl -n hw.logicalcpu 2>/dev/null || echo "8")
else
    AVAILABLE_CPUS=$(nproc 2>/dev/null || echo "8")
fi

# Set MAX_CPUS from parameter or use all available
if [ -n "${MAX_PROCESS}" ]; then
    MAX_CPUS="${MAX_PROCESS}"
    # Validate that MAX_CPUS is a positive integer
    if ! [[ "${MAX_CPUS}" =~ ^[0-9]+$ ]] || [ "${MAX_CPUS}" -lt 1 ]; then
        echo "ERROR: --max-process must be a positive integer"
        exit 1
    fi
    # Warn if exceeding available CPUs
    if [ "${MAX_CPUS}" -gt "${AVAILABLE_CPUS}" ]; then
        echo "WARNING: --max-process=${MAX_CPUS} exceeds available CPUs (${AVAILABLE_CPUS})"
        echo "         This may lead to oversubscription and reduced performance."
        echo ""
    fi
else
    MAX_CPUS="${AVAILABLE_CPUS}"
fi

echo "=============================================================================="
echo "3DPolyS-LE Parallelization Benchmark Suite"
echo "=============================================================================="
echo "Configuration:"
echo "  Config file:     ${CONFIG_FILE}"
echo "  Seed:            ${SEED}"
echo "  Output dir:      ${BENCH_DIR}"
echo "  Test selection:  ${TEST_SELECTION}"
echo "  Available CPUs:  ${AVAILABLE_CPUS}"
echo "  Max MPI procs:   ${MAX_CPUS}$([ "${MAX_CPUS}" != "${AVAILABLE_CPUS}" ] && echo " (custom)" || echo " (default)")"
echo "=============================================================================="
echo ""

# Display hardware configuration
detect_hardware

# Format seconds as human-readable time (round to min/h when appropriate)
format_elapsed() {
    local secs="$1"
    if [[ ! "$secs" =~ ^[0-9]*\.?[0-9]+$ ]]; then
        echo "N/A"
        return
    fi
    awk -v s="$secs" 'BEGIN {
        if (s >= 3600) printf "%.1f h\n", s/3600
        else if (s >= 60) printf "%.1f min\n", s/60
        else printf "%.1f s\n", s
    }'
}

# Output files used for reproducibility verification
REPRO_OUTPUT_FILES="config.out dr.out contact.out process.out Nlef.out"

# Function to run benchmark and capture timing
run_benchmark() {
    local test_name=$1
    local mpi_procs=$2
    local extra_flags=$3
    local output_subdir="${RESULTS_DIR}/${test_name}"
    local log_file="${LOG_DIR}/${test_name}.log"

    # Single-process runs still need mpirun -np 1 (the program calls MPI_Init)
    local np=1
    [ ${mpi_procs} -gt 0 ] && np=${mpi_procs}

    echo "Running: ${test_name}"
    echo "  MPI processes: ${np}"
    echo "  Extra flags: ${extra_flags}"
    echo "  Output: ${output_subdir}"

    # Use fixed OMP_NUM_THREADS for no-GPU (OpenMP) runs so OpenMP path is deterministic and comparable
    local env_prefix=""
    #if [[ "${extra_flags}" == *"--no-gpu"* ]]; then
    #    env_prefix="OMP_NUM_THREADS=4 "
    #    echo "  OpenMP: OMP_NUM_THREADS=4 (reproducibility)"
    #fi

    RUN_CMD="${env_prefix}/usr/bin/time -p mpirun -np ${np} ${EXECUTABLE} --seed:${SEED} -o:${output_subdir} ${extra_flags} ${CONFIG_FILE}"
    echo "  Command: ${RUN_CMD}"

    # Run the benchmark (always via mpirun so MPI_Init succeeds)
    eval "${env_prefix}/usr/bin/time -p mpirun -np ${np} ${EXECUTABLE} \
        --seed:${SEED} \
        -o:${output_subdir} \
        ${extra_flags} \
        ${CONFIG_FILE} \
        > \"${log_file}\" 2>&1"

    local exit_code=$?

    # Show elapsed time (from time -p "real" line), rounded to s/min/h
    local real_secs=""
    if [ -f "${log_file}" ]; then
        real_secs=$(grep '^real ' "${log_file}" 2>/dev/null | awk '{print $2}')
    fi
    if [ -n "${real_secs}" ]; then
        echo "  Time: $(format_elapsed "${real_secs}")"
    else
        echo "  Time: N/A"
    fi

    if [ ${exit_code} -eq 0 ]; then
        echo "  ✓ Completed successfully"

        # Extract timing information
        grep "Elapsed time" "${log_file}" | tail -2 >> "${LOG_DIR}/timing_summary.txt"
        echo "---" >> "${LOG_DIR}/timing_summary.txt"

        # Calculate checksums for all simulation outputs (reproducibility: GPU vs non-GPU must match)
        for f in config.out dr.out contact.out process.out Nlef.out; do
            if [ -f "${output_subdir}/${f}" ]; then
                printf "%s  %s  %s\n" "${test_name}" "${f}" "$(md5sum "${output_subdir}/${f}" | awk '{print $1}')" >> "${LOG_DIR}/checksums.txt"
            fi
        done
    else
        echo "  ✗ Failed with exit code ${exit_code}"
        echo "  See log: ${log_file}"
    fi

    echo ""
}

# Function to check if a test should run
# test_type: mpi-gpu, mpi-cpu, no-mpi-gpu, no-mpi-cpu
should_run_test() {
    local test_type=$1

    case "${TEST_SELECTION}" in
        all)
            return 0
            ;;
        mpi)
            # All tests that use MPI (multi-process: with GPU or without GPU)
            [[ "${test_type}" == "mpi-gpu" || "${test_type}" == "mpi-cpu" ]]
            ;;
        no-mpi|single)
            # Single process, no MPI (with GPU or without GPU)
            [[ "${test_type}" == "no-mpi-gpu" || "${test_type}" == "no-mpi-cpu" ]]
            ;;
        gpu)
            [[ "${test_type}" == *"gpu"* ]]
            ;;
        cpu)
            [[ "${test_type}" == *"cpu"* ]]
            ;;
        gpu-vs-no-gpu)
            # Only the two single-process tests that compare with vs without GPU
            [[ "${test_type}" == "no-mpi-gpu" || "${test_type}" == "no-mpi-cpu" ]]
            ;;
        *)
            return 0
            ;;
    esac
}

# Build test plan
echo "=============================================================================="
echo "TEST PLAN"
echo "=============================================================================="
echo "The following tests will be executed:"
echo ""

TEST_COUNT=0

# MPI tests (all use MPI; multi-process)
echo "  MPI (multi-process) tests:"
if should_run_test "mpi-gpu"; then
    echo "  [$((++TEST_COUNT))] MPI ${MAX_CPUS} procs, with GPU/OpenACC"
fi
if should_run_test "mpi-cpu"; then
    echo "  [$((++TEST_COUNT))] MPI ${MAX_CPUS} procs, without GPU/OpenACC (--no-gpu)"
fi

echo "  No-MPI (single process) tests:"
if should_run_test "no-mpi-gpu"; then
    echo "  [$((++TEST_COUNT))] With GPU/OpenACC, no MPI"
fi
if should_run_test "no-mpi-cpu"; then
    echo "  [$((++TEST_COUNT))] Without GPU (--no-gpu), no MPI"
fi

echo ""
echo "Total tests to run: ${TEST_COUNT}"
echo "=============================================================================="
echo ""

# Initialize timing summary file
echo "=== Timing Summary ===" > "${LOG_DIR}/timing_summary.txt"
echo "Benchmark run: ${TIMESTAMP}" >> "${LOG_DIR}/timing_summary.txt"
echo "Seed: ${SEED}" >> "${LOG_DIR}/timing_summary.txt"
echo "Test selection: ${TEST_SELECTION}" >> "${LOG_DIR}/timing_summary.txt"
echo "Max CPUs: ${MAX_CPUS}" >> "${LOG_DIR}/timing_summary.txt"
echo "" >> "${LOG_DIR}/timing_summary.txt"

# Initialize checksums file (all output files for reproducibility verification)
echo "=== MD5 Checksums (test_name  file  checksum) ===" > "${LOG_DIR}/checksums.txt"
echo "Format: each line is: <test_name>  <output_file>  <md5>" >> "${LOG_DIR}/checksums.txt"
echo "" >> "${LOG_DIR}/checksums.txt"

echo "Starting benchmark tests..."
echo ""

################################################################################
# MPI tests (multi-process)
################################################################################
# Test 1: MPI, with GPU/OpenACC
if should_run_test "mpi-gpu"; then
run_benchmark "mpi_${MAX_CPUS}proc_gpu" ${MAX_CPUS} ""
fi

# Test 2: MPI, without GPU/OpenACC (--no-gpu)
if should_run_test "mpi-cpu"; then
run_benchmark "mpi_${MAX_CPUS}proc_no_gpu" ${MAX_CPUS} "--no-gpu"
fi

################################################################################
# No-MPI tests (single process) — direct GPU vs no-GPU comparison
################################################################################
# Test 3: With GPU/OpenACC, no MPI
if should_run_test "no-mpi-gpu"; then
run_benchmark "gpu_no_mpi" 0 ""
fi

# Test 4: Without GPU (--no-gpu), no MPI
if should_run_test "no-mpi-cpu"; then
run_benchmark "no_gpu_no_mpi" 0 "--no-gpu"
fi

################################################################################
# Summary and verification
################################################################################
echo "=============================================================================="
echo "Benchmark completed!"
echo "=============================================================================="
echo ""
echo "Results location: ${BENCH_DIR}"
echo ""

# Check reproducibility (all output files; same seed => identical across GPU, OpenMP, MPI, no-MPI)
echo "Verifying reproducibility (MD5 checksums):"
cat "${LOG_DIR}/checksums.txt"
echo ""

# Helper: compare two result dirs; return 0 if all REPRO_OUTPUT_FILES match
compare_result_dirs() {
    local dir1="$1"
    local dir2="$2"
    local same=1
    for outfile in ${REPRO_OUTPUT_FILES}; do
        [ -f "${dir1}/${outfile}" ] && [ -f "${dir2}/${outfile}" ] || continue
        c1=$(md5sum "${dir1}/${outfile}" 2>/dev/null | awk '{print $1}')
        c2=$(md5sum "${dir2}/${outfile}" 2>/dev/null | awk '{print $1}')
        if [ -n "${c1}" ] && [ -n "${c2}" ] && [ "${c1}" != "${c2}" ]; then
            same=0
            break
        fi
    done
    return ${same}
}

REPRO_FAIL=0
for outfile in ${REPRO_OUTPUT_FILES}; do
    unique=$(awk -v f="${outfile}" '$2 == f {print $3}' "${LOG_DIR}/checksums.txt" | sort -u | wc -l)
    total=$(awk -v f="${outfile}" '$2 == f {print $3}' "${LOG_DIR}/checksums.txt" | wc -l)
    if [ "${total}" -eq 0 ]; then
        continue
    fi
    if [ "${unique}" -eq 1 ]; then
        echo "  ${outfile}: ✓ identical across all runs (${total} runs)"
    else
        echo "  ${outfile}: ⚠ ${unique} distinct checksums (expected 1) - reproducibility NOT verified"
        REPRO_FAIL=1
    fi
done

# Pairwise reproducibility checks (when both sides exist)
# GPU vs no-GPU (OpenMP) — single process
if [ -f "${RESULTS_DIR}/gpu_no_mpi/config.out" ] && [ -f "${RESULTS_DIR}/no_gpu_no_mpi/config.out" ]; then
    if compare_result_dirs "${RESULTS_DIR}/gpu_no_mpi" "${RESULTS_DIR}/no_gpu_no_mpi"; then
        echo "  GPU vs no-GPU (OpenMP, same seed): ✓ all output files byte-identical"
    else
        echo "  GPU vs no-GPU (OpenMP): ⚠ some output files differ"
        REPRO_FAIL=1
    fi
fi

# MPI vs no-MPI (GPU): multi-process GPU vs single-process GPU
mpi_gpu_dir="${RESULTS_DIR}/mpi_${MAX_CPUS}proc_gpu"
if [ -f "${mpi_gpu_dir}/config.out" ] && [ -f "${RESULTS_DIR}/gpu_no_mpi/config.out" ]; then
    if compare_result_dirs "${mpi_gpu_dir}" "${RESULTS_DIR}/gpu_no_mpi"; then
        echo "  MPI vs no-MPI (GPU, same seed): ✓ all output files byte-identical"
    else
        echo "  MPI vs no-MPI (GPU): ⚠ some output files differ"
        REPRO_FAIL=1
    fi
fi

# MPI vs no-MPI (no-GPU / OpenMP): multi-process CPU vs single-process OpenMP
mpi_no_gpu_dir="${RESULTS_DIR}/mpi_${MAX_CPUS}proc_no_gpu"
if [ -f "${mpi_no_gpu_dir}/config.out" ] && [ -f "${RESULTS_DIR}/no_gpu_no_mpi/config.out" ]; then
    if compare_result_dirs "${mpi_no_gpu_dir}" "${RESULTS_DIR}/no_gpu_no_mpi"; then
        echo "  MPI vs no-MPI (OpenMP/CPU, same seed): ✓ all output files byte-identical"
    else
        echo "  MPI vs no-MPI (OpenMP/CPU): ⚠ some output files differ"
        REPRO_FAIL=1
    fi
fi

if [ ${REPRO_FAIL} -eq 0 ]; then
    echo "✓ REPRODUCIBILITY VERIFIED: All runs produced identical results for the same seed"
    echo "  (GPU, OpenMP/no-GPU, MPI+GPU, MPI+no-GPU)"
else
    echo "⚠ WARNING: Some outputs differ across runs (GPU/OpenMP/MPI/no-MPI)"
fi
echo ""

# Display timing summary
echo "Timing summary:"
cat "${LOG_DIR}/timing_summary.txt"
echo ""

echo "=============================================================================="
echo "Generate report:"
echo ""
echo "  ./generate_benchmark_report.sh ${BENCH_DIR}"
echo ""
echo "=============================================================================="

./generate_benchmark_report.sh "${BENCH_DIR}"