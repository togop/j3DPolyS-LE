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
#                        no-mpi   - Run only non-MPI tests (GPU/CPU only)
#                        single   - Run single MPI process tests
#                        multi    - Run multi-process MPI tests
#                        gpu      - Run GPU-based tests only
#                        cpu      - Run CPU-based tests only
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
#   ./benchmark_parallelization.sh --tests=no-mpi
#   ./benchmark_parallelization.sh --max-process=8
#   ./benchmark_parallelization.sh --tests=gpu --max-process=16 ./test/test_tads_init_benchmark.cfg 42 ./my_benchmarks
#
################################################################################

set -e  # Exit on error

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

# Function to run benchmark and capture timing
run_benchmark() {
    local test_name=$1
    local mpi_procs=$2
    local extra_flags=$3
    local threads=$4  # New parameter for thread count (optional, 0 = auto)
    local output_subdir="${RESULTS_DIR}/${test_name}"
    local log_file="${LOG_DIR}/${test_name}.log"

    # Add thread flag if specified and > 0
    local thread_flag=""
    if [ -n "${threads}" ] && [ ${threads} -gt 0 ]; then
        thread_flag="--threads:${threads}"
    fi

    echo "Running: ${test_name}"
    echo "  MPI processes: ${mpi_procs}"
    echo "  Threads: ${threads:-auto}"
    echo "  Extra flags: ${extra_flags} ${thread_flag}"
    echo "  Output: ${output_subdir}"

    # Run the benchmark
    if [ ${mpi_procs} -eq 0 ]; then
        # Run without MPI (direct execution)
        /usr/bin/time -p ${EXECUTABLE} \
            --seed:${SEED} \
            -o:${output_subdir} \
            ${extra_flags} \
            ${thread_flag} \
            ${CONFIG_FILE} \
            > "${log_file}" 2>&1
    else
        # Run with MPI
        /usr/bin/time -p mpirun -np ${mpi_procs} ${EXECUTABLE} \
            --seed:${SEED} \
            -o:${output_subdir} \
            ${extra_flags} \
            ${thread_flag} \
            ${CONFIG_FILE} \
            > "${log_file}" 2>&1
    fi

    local exit_code=$?

    if [ ${exit_code} -eq 0 ]; then
        echo "  ✓ Completed successfully"

        # Extract timing information
        grep "Elapsed time" "${log_file}" | tail -2 >> "${LOG_DIR}/timing_summary.txt"
        echo "---" >> "${LOG_DIR}/timing_summary.txt"

        # Calculate checksums for verification
        if [ -f "${output_subdir}/config.out" ]; then
            md5sum "${output_subdir}/config.out" >> "${LOG_DIR}/checksums.txt"
        fi
    else
        echo "  ✗ Failed with exit code ${exit_code}"
        echo "  See log: ${log_file}"
    fi

    echo ""
}

# Function to check if a test should run
should_run_test() {
    local test_type=$1  # no-mpi, single, multi, gpu, cpu

    case "${TEST_SELECTION}" in
        all)
            return 0
            ;;
        no-mpi)
            [[ "${test_type}" == "no-mpi" ]]
            ;;
        single)
            [[ "${test_type}" == "single" ]]
            ;;
        multi)
            [[ "${test_type}" == "multi" ]]
            ;;
        gpu)
            [[ "${test_type}" == *"gpu"* ]]
            ;;
        cpu)
            [[ "${test_type}" == *"cpu"* ]]
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

# Multi-process MPI tests - GPU
if should_run_test "multi-gpu"; then
    echo "  [$((++TEST_COUNT))] MPI ${MAX_CPUS} processes (GPU/OpenACC) - 1 thread"
    echo "  [$((++TEST_COUNT))] MPI ${MAX_CPUS} processes (GPU/OpenACC) - auto threads"
fi

# Multi-process MPI tests - CPU
if should_run_test "multi-cpu"; then
    echo "  [$((++TEST_COUNT))] MPI ${MAX_CPUS} processes (OpenMP/CPU) - 1 thread"
    echo "  [$((++TEST_COUNT))] MPI ${MAX_CPUS} processes (OpenMP/CPU) - auto threads"
fi

# Multi-process MPI tests - No acceleration
if should_run_test "multi"; then
    echo "  [$((++TEST_COUNT))] MPI ${MAX_CPUS} processes (no OpenACC, no OpenMP)"
fi

# No-MPI tests (moved to end)
if should_run_test "no-mpi-gpu"; then
    echo "  [$((++TEST_COUNT))] GPU/OpenACC (no MPI) - auto threads"
fi
if should_run_test "no-mpi-cpu"; then
    echo "  [$((++TEST_COUNT))] CPU only (no MPI) - auto threads"
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

# Initialize checksums file
echo "=== MD5 Checksums (config.out) ===" > "${LOG_DIR}/checksums.txt"
echo "" >> "${LOG_DIR}/checksums.txt"

echo "Starting benchmark tests..."
echo ""

################################################################################
# Test 1: MPI MAX_CPUS processes (GPU/OpenACC) - 1 thread
################################################################################
if should_run_test "multi-gpu"; then
run_benchmark "mpi_${MAX_CPUS}proc_gpu_t1" ${MAX_CPUS} "" 1

################################################################################
# Test 2: MPI MAX_CPUS processes (GPU/OpenACC) - auto threads
################################################################################
run_benchmark "mpi_${MAX_CPUS}proc_gpu" ${MAX_CPUS} "" 0
fi

################################################################################
# Test 3: MPI MAX_CPUS processes (OpenMP/CPU) - 1 thread
################################################################################
if should_run_test "multi-cpu"; then
run_benchmark "mpi_${MAX_CPUS}proc_openmp_t1" ${MAX_CPUS} "--no-gpu-prefer" 1

################################################################################
# Test 4: MPI MAX_CPUS processes (OpenMP/CPU) - auto threads
################################################################################
run_benchmark "mpi_${MAX_CPUS}proc_openmp" ${MAX_CPUS} "--no-gpu-prefer" 0
fi

################################################################################
# Test 5: MPI MAX_CPUS processes (no OpenACC, no OpenMP)
################################################################################
if should_run_test "multi"; then
run_benchmark "mpi_${MAX_CPUS}proc_plain" ${MAX_CPUS} "--no-gpu-prefer --no-openmp" 0
fi

################################################################################
# Test 6: GPU/OpenACC (no MPI) - auto threads
################################################################################
if should_run_test "no-mpi-gpu"; then
run_benchmark "gpu_no_mpi" 0 "" 0
fi

################################################################################
# Test 7: CPU only (no MPI) - auto threads
################################################################################
if should_run_test "no-mpi-cpu"; then
run_benchmark "cpu_no_mpi" 0 "--no-gpu-prefer" 0
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

# Check reproducibility
echo "Verifying reproducibility (MD5 checksums):"
cat "${LOG_DIR}/checksums.txt"
echo ""

# Count unique checksums (excluding header lines)
unique_checksums=$(grep -v "===" "${LOG_DIR}/checksums.txt" | grep "config.out" | awk '{print $1}' | sort -u | wc -l)
echo "Unique checksums: ${unique_checksums}"

if [ ${unique_checksums} -eq 1 ]; then
    echo "✓ REPRODUCIBILITY VERIFIED: All runs produced identical results"
else
    echo "⚠ WARNING: Different checksums detected - results may not be reproducible"
fi
echo ""

# Display timing summary
echo "Timing summary:"
cat "${LOG_DIR}/timing_summary.txt"
echo ""

echo "=============================================================================="
echo "To generate the report, run:"
echo "  ./generate_benchmark_report.sh ${BENCH_DIR}"
echo "=============================================================================="
