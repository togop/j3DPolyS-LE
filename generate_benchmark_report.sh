#!/bin/bash
################################################################################
# Generate Benchmark Report from Results
################################################################################
#
# This script generates a markdown report from benchmark results.
#
# Usage:
#   ./generate_benchmark_report.sh <benchmark_dir>
#
# Arguments:
#   benchmark_dir - Directory containing benchmark results
#
################################################################################

set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 <benchmark_dir>"
    echo "Example: $0 ./benchmark_results/20260207_123456"
    exit 1
fi

BENCH_DIR="$1"
LOG_DIR="${BENCH_DIR}/logs"
RESULTS_DIR="${BENCH_DIR}/results"
REPORT_FILE="${BENCH_DIR}/BENCHMARK_REPORT.md"

if [ ! -d "${BENCH_DIR}" ]; then
    echo "ERROR: Benchmark directory not found: ${BENCH_DIR}"
    exit 1
fi

echo "Generating benchmark report..."
echo "  Input: ${BENCH_DIR}"
echo "  Output: ${REPORT_FILE}"

# Extract system information
HOSTNAME=$(hostname)
OS=$(uname -s)
OS_VERSION=$(uname -r)
ARCH=$(uname -m)
NPROC=$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo "Unknown")

# Extract timing data
extract_timing() {
    local test_name=$1
    local log_file="${LOG_DIR}/${test_name}.log"

    if [ -f "${log_file}" ]; then
        sim_time=$(grep "END Running similations Elapsed time" "${log_file}" | tail -1 | awk '{print $NF}')
        merge_time=$(grep "Merged files Elapsed time" "${log_file}" | tail -1 | awk '{print $NF}')
        method=$(grep "Method:" "${log_file}" | tail -1 | awk -F'Method:' '{print $2}' | tr -d ' ')

        echo "${sim_time:-N/A}|${merge_time:-N/A}|${method:-N/A}"
    else
        echo "N/A|N/A|N/A"
    fi
}

# Generate report
cat > "${REPORT_FILE}" << 'EOF_HEADER'
# 3DPolyS-LE Parallelization Methods Benchmark Report

## Executive Summary

This report presents performance benchmarks for different parallelization methods in 3DPolyS-LE, comparing GPU/OpenACC, OpenMP/CPU, and MPI scaling with various process counts. All tests used the same seed value to verify reproducibility.

---

## System Information

EOF_HEADER

# Add system information
cat >> "${REPORT_FILE}" << EOF
- **Hostname:** ${HOSTNAME}
- **Operating System:** ${OS} ${OS_VERSION}
- **Architecture:** ${ARCH}
- **CPU Cores:** ${NPROC}
- **Test Date:** $(date)
- **Benchmark Directory:** ${BENCH_DIR}

---

## Test Configuration

EOF

# Extract config parameters from log
if [ -f "${LOG_DIR}/gpu_default.log" ]; then
    CONFIG_FILE=$(grep "Load initial parameters" "${LOG_DIR}/gpu_default.log" | head -1 | awk -F'from ' '{print $2}' | awk '{print $1}')
    SEED=$(grep "Random seed set to:" "${LOG_DIR}/gpu_default.log" | head -1 | awk '{print $NF}')

    cat >> "${REPORT_FILE}" << EOF
**Configuration File:** \`${CONFIG_FILE}\`

**Seed Value:** ${SEED}

**Parameters** (extracted from simulation logs):
- Seed: ${SEED}
- All other parameters from configuration file

---

## Benchmark Results

### Performance Summary Table

| Test Configuration | MPI Procs | Parallelization | Simulation Time (s) | Merge Time (s) | Method |
|-------------------|-----------|-----------------|---------------------|----------------|--------|
EOF

    # Extract timing for each test
    for test in gpu_default openmp_cpu mpi_2proc_gpu mpi_4proc_gpu mpi_8proc_gpu mpi_2proc_openmp mpi_4proc_openmp mpi_8proc_openmp; do
        timing=$(extract_timing "${test}")
        sim_time=$(echo "${timing}" | cut -d'|' -f1)
        merge_time=$(echo "${timing}" | cut -d'|' -f2)
        method=$(echo "${timing}" | cut -d'|' -f3)

        case ${test} in
            gpu_default)
                desc="GPU/OpenACC (default)"
                procs="1"
                para="GPU/OpenACC"
                ;;
            openmp_cpu)
                desc="OpenMP (CPU)"
                procs="1"
                para="OpenMP/CPU"
                ;;
            mpi_2proc_gpu)
                desc="MPI 2 processes (GPU)"
                procs="2"
                para="GPU/OpenACC"
                ;;
            mpi_4proc_gpu)
                desc="MPI 4 processes (GPU)"
                procs="4"
                para="GPU/OpenACC"
                ;;
            mpi_8proc_gpu)
                desc="MPI 8 processes (GPU)"
                procs="8"
                para="GPU/OpenACC"
                ;;
            mpi_2proc_openmp)
                desc="MPI 2 processes (OpenMP)"
                procs="2"
                para="OpenMP/CPU"
                ;;
            mpi_4proc_openmp)
                desc="MPI 4 processes (OpenMP)"
                procs="4"
                para="OpenMP/CPU"
                ;;
            mpi_8proc_openmp)
                desc="MPI 8 processes (OpenMP)"
                procs="8"
                para="OpenMP/CPU"
                ;;
        esac

        echo "| ${desc} | ${procs} | ${para} | ${sim_time} | ${merge_time} | ${method} |" >> "${REPORT_FILE}"
    done
fi

# Add reproducibility verification
cat >> "${REPORT_FILE}" << 'EOF'

---

## Reproducibility Verification

### MD5 Checksums

EOF

cat "${LOG_DIR}/checksums.txt" >> "${REPORT_FILE}"

cat >> "${REPORT_FILE}" << 'EOF'

### Verification Result

EOF

# Check checksums (excluding header lines)
unique_checksums=$(grep -v "===" "${LOG_DIR}/checksums.txt" | grep "config.out" | awk '{print $1}' | sort -u | wc -l)

if [ ${unique_checksums} -eq 1 ]; then
    cat >> "${REPORT_FILE}" << 'EOF'
✓ **REPRODUCIBILITY VERIFIED**

All test runs produced byte-for-byte identical results. The MD5 checksums for all output files are identical, confirming perfect reproducibility across all parallelization methods when using the same seed.

EOF
else
    cat >> "${REPORT_FILE}" << 'EOF'
⚠ **WARNING: Different Checksums Detected**

The test runs produced different checksums. This may indicate:
- Non-deterministic behavior in the simulation
- Different random number sequences
- Hardware-specific variations

Please review the logs for more details.

EOF
fi

# Add timing summary
cat >> "${REPORT_FILE}" << 'EOF'
---

## Detailed Timing Information

EOF

cat "${LOG_DIR}/timing_summary.txt" >> "${REPORT_FILE}"

# Add footer
cat >> "${REPORT_FILE}" << 'EOF'

---

## How to Use This Report

This benchmark report provides:
1. **Performance comparison** across different parallelization methods
2. **Reproducibility verification** via MD5 checksums
3. **Detailed timing information** for each test configuration

### Interpreting Results

- **Simulation Time:** Core computation time (Monte Carlo + LEF dynamics)
- **Merge Time:** File I/O time for combining MPI process outputs
- **Method:** Parallelization method used (OpenACC/GPU or OpenMP/CPU)

### Recommendations

Based on your results:
- Compare single-process GPU vs OpenMP performance
- Evaluate MPI scaling efficiency (speedup vs number of processes)
- Use the configuration with best performance for your problem size

---

## Files Included

- `logs/` - Complete log files for each test
- `results/` - Simulation output files for each test
- `logs/checksums.txt` - MD5 checksums for reproducibility verification
- `logs/timing_summary.txt` - Extracted timing information

---

**Report generated:** $(date)

EOF

echo "✓ Report generated successfully: ${REPORT_FILE}"
echo ""
echo "To view the report:"
echo "  cat ${REPORT_FILE}"
echo "  # or"
echo "  open ${REPORT_FILE}  # on macOS"
