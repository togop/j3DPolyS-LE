#!/bin/bash
################################################################################
# 3DPolyS-LE Parallelization Methods Benchmark Script
################################################################################
#
# This script benchmarks different parallelization methods for 3DPolyS-LE
# simulations to compare performance and verify reproducibility.
#
# Usage:
#   ./benchmark_parallelization.sh [config_file] [seed] [output_dir]
#
# Arguments:
#   config_file  - Path to simulation config file (default: ./test/test_tads_init_benchmark.cfg)
#   seed         - Random seed for reproducibility (default: 42)
#   output_dir   - Output directory for results (default: ./benchmark_results)
#
# Example:
#   ./benchmark_parallelization.sh ./test/test_tads_init_benchmark.cfg 42 ./my_benchmarks
#
################################################################################

set -e  # Exit on error

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

echo "=============================================================================="
echo "3DPolyS-LE Parallelization Benchmark Suite"
echo "=============================================================================="
echo "Configuration:"
echo "  Config file: ${CONFIG_FILE}"
echo "  Seed:        ${SEED}"
echo "  Output dir:  ${BENCH_DIR}"
echo "  Date:        $(date)"
echo "=============================================================================="
echo ""

# Function to run benchmark and capture timing
run_benchmark() {
    local test_name=$1
    local mpi_procs=$2
    local extra_flags=$3
    local output_subdir="${RESULTS_DIR}/${test_name}"
    local log_file="${LOG_DIR}/${test_name}.log"

    echo "Running: ${test_name}"
    echo "  MPI processes: ${mpi_procs}"
    echo "  Extra flags: ${extra_flags}"
    echo "  Output: ${output_subdir}"

    # Run the benchmark
    /usr/bin/time -p mpirun -np ${mpi_procs} ${EXECUTABLE} \
        --seed:${SEED} \
        -o:${output_subdir} \
        ${extra_flags} \
        ${CONFIG_FILE} \
        > "${log_file}" 2>&1

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

# Start benchmarking
echo "Starting benchmark tests..."
echo ""

# Initialize timing summary file
echo "=== Timing Summary ===" > "${LOG_DIR}/timing_summary.txt"
echo "Benchmark run: ${TIMESTAMP}" >> "${LOG_DIR}/timing_summary.txt"
echo "Seed: ${SEED}" >> "${LOG_DIR}/timing_summary.txt"
echo "" >> "${LOG_DIR}/timing_summary.txt"

# Initialize checksums file
echo "=== MD5 Checksums (config.out) ===" > "${LOG_DIR}/checksums.txt"
echo "" >> "${LOG_DIR}/checksums.txt"

################################################################################
# Test 1: GPU/OpenACC (default) - Single MPI process
################################################################################
run_benchmark "gpu_default" 1 ""

################################################################################
# Test 2: OpenMP (CPU) - Single MPI process
################################################################################
run_benchmark "openmp_cpu" 1 "--no-gpu-prefer"

################################################################################
# Test 3: MPI 2 processes (GPU)
################################################################################
run_benchmark "mpi_2proc_gpu" 2 ""

################################################################################
# Test 4: MPI 4 processes (GPU)
################################################################################
run_benchmark "mpi_4proc_gpu" 4 ""

################################################################################
# Test 5: MPI 8 processes (GPU)
################################################################################
run_benchmark "mpi_8proc_gpu" 8 ""

################################################################################
# Test 6: MPI 2 processes (OpenMP/CPU)
################################################################################
run_benchmark "mpi_2proc_openmp" 2 "--no-gpu-prefer"

################################################################################
# Test 7: MPI 4 processes (OpenMP/CPU)
################################################################################
run_benchmark "mpi_4proc_openmp" 4 "--no-gpu-prefer"

################################################################################
# Test 8: MPI 8 processes (OpenMP/CPU)
################################################################################
run_benchmark "mpi_8proc_openmp" 8 "--no-gpu-prefer"

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
