#!/usr/bin/env bash
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

# Convert seconds to human-readable format
format_time() {
    local seconds=$1

    # Handle N/A or invalid input
    if [ "${seconds}" = "N/A" ] || [ -z "${seconds}" ]; then
        echo "N/A"
        return
    fi

    # Convert to float for bc
    local time_val=$(echo "${seconds}" | awk '{printf "%.3f", $1}')

    # Check if >= 3600 seconds (1 hour)
    if (( $(echo "${time_val} >= 3600" | bc -l) )); then
        local hours=$(echo "${time_val} / 3600" | bc -l | awk '{printf "%.2f", $1}')
        echo "${hours} h"
    # Check if >= 60 seconds (1 minute)
    elif (( $(echo "${time_val} >= 60" | bc -l) )); then
        local minutes=$(echo "${time_val} / 60" | bc -l | awk '{printf "%.2f", $1}')
        echo "${minutes} min"
    else
        echo "${time_val} s"
    fi
}

# Extract timing data
extract_timing() {
    local test_name=$1
    local log_file="${LOG_DIR}/${test_name}.log"

    if [ -f "${log_file}" ]; then
        sim_time=$(grep "END Running similations Elapsed time" "${log_file}" | tail -1 | awk '{print $NF}')
        merge_time=$(grep "Merged files Elapsed time" "${log_file}" | tail -1 | awk '{print $NF}')
        method=$(grep "Method:" "${log_file}" | tail -1 | awk -F'Method:' '{print $2}' | tr -d ' ')

        # Extract thread count information
        threads="auto"
        if grep -q "OpenMP threads set to:" "${log_file}"; then
            threads=$(grep "OpenMP threads set to:" "${log_file}" | head -1 | awk '{print $NF}')
        elif grep -q "OpenMP enabled with default thread count:" "${log_file}"; then
            threads=$(grep "OpenMP enabled with default thread count:" "${log_file}" | head -1 | awk '{print $NF}')
            threads="${threads} (auto)"
        elif [ "${method}" = "OpenACC" ]; then
            threads="N/A (GPU)"
        fi

        echo "${sim_time:-N/A}|${merge_time:-N/A}|${method:-N/A}|${threads:-auto}"
    else
        echo "N/A|N/A|N/A|N/A"
    fi
}

# Generate report
cat > "${REPORT_FILE}" << 'EOF_HEADER'
# 3DPolyS-LE Parallelization Methods Benchmark Report

## Executive Summary

This report presents performance benchmarks for different parallelization methods in 3DPolyS-LE, comparing:
- **No MPI**: Direct execution using GPU or CPU only
- **Single MPI Process**: GPU/OpenACC and OpenMP/CPU with single process
- **Multi-process MPI**: Scaling with 2, 4, and 8 MPI processes

All tests used the same seed value to verify reproducibility.

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

# Extract config parameters from log (check both no-MPI and MPI versions)
LOG_FILE=""
if [ -f "${LOG_DIR}/gpu_no_mpi.log" ]; then
    LOG_FILE="${LOG_DIR}/gpu_no_mpi.log"
elif [ -f "${LOG_DIR}/gpu_default.log" ]; then
    LOG_FILE="${LOG_DIR}/gpu_default.log"
fi

if [ -n "${LOG_FILE}" ]; then
    CONFIG_FILE=$(grep "Load initial parameters" "${LOG_FILE}" | head -1 | awk -F'from ' '{print $2}' | awk '{print $1}')
    SEED=$(grep "Random seed set to:" "${LOG_FILE}" | head -1 | awk '{print $NF}')

    cat >> "${REPORT_FILE}" << EOF
**Configuration File:** \`${CONFIG_FILE}\`

**Seed Value:** ${SEED}

**Parameters** (extracted from simulation logs):
- Seed: ${SEED}
- All other parameters from configuration file

---

## Benchmark Results

### Performance Summary Table

| Test Configuration | MPI Procs | Parallelization | Threads | Simulation Time | Method |
|-------------------|-----------|-----------------|---------|----------------|--------|
EOF

    # Store timing data to temp files for scaling analysis
    TEMP_DIR=$(mktemp -d)
    trap "rm -rf ${TEMP_DIR}" EXIT

    # Extract timing for each test - find all test directories dynamically
    test_list=$(find "${RESULTS_DIR}" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort)

    # Extract timing for each test
    for test in ${test_list}; do
        timing=$(extract_timing "${test}")
        sim_time_raw=$(echo "${timing}" | cut -d'|' -f1)
        merge_time_raw=$(echo "${timing}" | cut -d'|' -f2)
        method=$(echo "${timing}" | cut -d'|' -f3)
        threads=$(echo "${timing}" | cut -d'|' -f4)

        # Store raw times for scaling calculations
        echo "${sim_time_raw}" > "${TEMP_DIR}/${test}.time"

        # Store baselines (prefer no-MPI versions if available)
        if [ "${test}" = "gpu_no_mpi" ]; then
            echo "${sim_time_raw}" > "${TEMP_DIR}/gpu_baseline.time"
        elif [ "${test}" = "gpu_default" ] && [ ! -f "${TEMP_DIR}/gpu_baseline.time" ]; then
            echo "${sim_time_raw}" > "${TEMP_DIR}/gpu_baseline.time"
        fi

        if [ "${test}" = "cpu_no_mpi" ]; then
            echo "${sim_time_raw}" > "${TEMP_DIR}/openmp_baseline.time"
        elif [ "${test}" = "openmp_cpu" ] && [ ! -f "${TEMP_DIR}/openmp_baseline.time" ]; then
            echo "${sim_time_raw}" > "${TEMP_DIR}/openmp_baseline.time"
        fi

        # Format times for display
        sim_time=$(format_time "${sim_time_raw}")

        # Parse test name to extract description and process count
        desc="${test}"
        procs="1"
        para="${method}"

        # Determine parallelization type from method
        if [ "${method}" = "OpenACC" ]; then
            para="GPU/OpenACC"
        elif [ "${method}" = "OpenMP" ]; then
            para="OpenMP/CPU"
        fi

        # Extract MPI process count from test name
        if [[ "${test}" =~ mpi_([0-9]+)proc ]]; then
            procs="${BASH_REMATCH[1]}"
        elif [[ "${test}" =~ no_mpi ]]; then
            procs="0"
        fi

        # Create human-readable description
        if [ "${test}" = "gpu_no_mpi" ]; then
            desc="GPU only (no MPI)"
        elif [ "${test}" = "cpu_no_mpi_t8" ]; then
            desc="CPU only (no MPI) - 8 threads"
        elif [[ "${test}" =~ cpu_no_mpi_t([0-9]+) ]]; then
            desc="CPU only (no MPI) - ${BASH_REMATCH[1]} threads"
        elif [ "${test}" = "cpu_no_mpi" ]; then
            desc="CPU only (no MPI) - auto threads"
        elif [ "${test}" = "gpu_default" ]; then
            desc="GPU/OpenACC (default)"
        elif [[ "${test}" =~ openmp_cpu_t([0-9]+) ]]; then
            desc="OpenMP (CPU) - ${BASH_REMATCH[1]} threads"
        elif [ "${test}" = "openmp_cpu" ]; then
            desc="OpenMP (CPU) - auto threads"
        elif [[ "${test}" =~ mpi_([0-9]+)proc_openmp_t([0-9]+) ]]; then
            desc="MPI ${BASH_REMATCH[1]} procs (OpenMP) - ${BASH_REMATCH[2]} threads"
        elif [[ "${test}" =~ mpi_([0-9]+)proc_openmp ]]; then
            desc="MPI ${BASH_REMATCH[1]} procs (OpenMP)"
        elif [[ "${test}" =~ mpi_([0-9]+)proc_gpu ]]; then
            desc="MPI ${BASH_REMATCH[1]} procs (GPU)"
        fi

        # Display MPI process count (show "None" for non-MPI tests)
        mpi_display="${procs}"
        if [ "${procs}" = "0" ]; then
            mpi_display="None"
        fi

        echo "| ${desc} | ${mpi_display} | ${para} | ${threads} | ${sim_time} | ${method} |" >> "${REPORT_FILE}"
    done

    # Add visual performance summary
    cat >> "${REPORT_FILE}" << 'EOF'

---

## Visual Performance Summary

### GPU/OpenACC Performance
EOF

    # GPU scaling chart
    gpu_baseline=$(cat "${TEMP_DIR}/gpu_baseline.time" 2>/dev/null || echo "N/A")
    if [ "${gpu_baseline}" != "N/A" ] && [ -n "${gpu_baseline}" ]; then
        echo "" >> "${REPORT_FILE}"
        echo "\`\`\`" >> "${REPORT_FILE}"
        echo "Baseline (1 proc):  $(format_time ${gpu_baseline})" >> "${REPORT_FILE}"

        for procs in 2 4 8; do
            test_name="mpi_${procs}proc_gpu"
            time_raw=$(cat "${TEMP_DIR}/${test_name}.time" 2>/dev/null || echo "N/A")
            if [ "${time_raw}" != "N/A" ] && [ -n "${time_raw}" ]; then
                speedup=$(echo "scale=2; ${gpu_baseline} / ${time_raw}" | bc -l)
                efficiency=$(echo "scale=1; (${speedup} / ${procs}) * 100" | bc -l)
                bar_length=$(echo "${speedup} * 10" | bc -l | awk '{printf "%d", $1}')
                bar=$(printf '█%.0s' $(seq 1 ${bar_length} 2>/dev/null || echo ""))
                printf "MPI %d procs:      %s  %s (speedup: %.2fx, efficiency: %.1f%%)\n" \
                    ${procs} "$(format_time ${time_raw})" "${bar}" ${speedup} ${efficiency} >> "${REPORT_FILE}"
            fi
        done
        echo "\`\`\`" >> "${REPORT_FILE}"
    fi

    cat >> "${REPORT_FILE}" << 'EOF'

### OpenMP/CPU Performance
EOF

    # OpenMP scaling chart
    openmp_baseline=$(cat "${TEMP_DIR}/openmp_baseline.time" 2>/dev/null || echo "N/A")
    if [ "${openmp_baseline}" != "N/A" ] && [ -n "${openmp_baseline}" ]; then
        echo "" >> "${REPORT_FILE}"
        echo "\`\`\`" >> "${REPORT_FILE}"
        echo "Baseline (1 proc):  $(format_time ${openmp_baseline})" >> "${REPORT_FILE}"

        for procs in 2 4 8; do
            test_name="mpi_${procs}proc_openmp"
            time_raw=$(cat "${TEMP_DIR}/${test_name}.time" 2>/dev/null || echo "N/A")
            if [ "${time_raw}" != "N/A" ] && [ -n "${time_raw}" ]; then
                speedup=$(echo "scale=2; ${openmp_baseline} / ${time_raw}" | bc -l)
                efficiency=$(echo "scale=1; (${speedup} / ${procs}) * 100" | bc -l)
                bar_length=$(echo "${speedup} * 10" | bc -l | awk '{printf "%d", $1}')
                bar=$(printf '█%.0s' $(seq 1 ${bar_length} 2>/dev/null || echo ""))
                printf "MPI %d procs:      %s  %s (speedup: %.2fx, efficiency: %.1f%%)\n" \
                    ${procs} "$(format_time ${time_raw})" "${bar}" ${speedup} ${efficiency} >> "${REPORT_FILE}"
            fi
        done
        echo "\`\`\`" >> "${REPORT_FILE}"
    fi

    # Add scaling efficiency table
    cat >> "${REPORT_FILE}" << 'EOF'

---

## Scaling Efficiency Analysis

### MPI Scaling Metrics

| Configuration | Processes | Simulation Time | Speedup | Parallel Efficiency | Scaling Category |
|--------------|-----------|----------------|---------|-------------------|-----------------|
EOF

    # GPU scaling efficiency
    gpu_baseline=$(cat "${TEMP_DIR}/gpu_baseline.time" 2>/dev/null || echo "N/A")
    if [ "${gpu_baseline}" != "N/A" ] && [ -n "${gpu_baseline}" ]; then
        for procs in 1 2 4 8; do
            if [ ${procs} -eq 1 ]; then
                test_name="gpu_default"
            else
                test_name="mpi_${procs}proc_gpu"
            fi

            time_raw=$(cat "${TEMP_DIR}/${test_name}.time" 2>/dev/null || echo "N/A")
            if [ "${time_raw}" != "N/A" ] && [ -n "${time_raw}" ]; then
                if [ ${procs} -eq 1 ]; then
                    speedup="1.00"
                    efficiency="100.0"
                    category="Baseline"
                else
                    speedup=$(echo "scale=2; ${gpu_baseline} / ${time_raw}" | bc -l)
                    efficiency=$(echo "scale=1; (${speedup} / ${procs}) * 100" | bc -l)

                    # Categorize scaling
                    if (( $(echo "${efficiency} >= 90" | bc -l) )); then
                        category="Excellent (≥90%)"
                    elif (( $(echo "${efficiency} >= 75" | bc -l) )); then
                        category="Good (75-90%)"
                    elif (( $(echo "${efficiency} >= 50" | bc -l) )); then
                        category="Moderate (50-75%)"
                    else
                        category="Poor (<50%)"
                    fi
                fi

                time_fmt=$(format_time "${time_raw}")
                echo "| GPU/OpenACC | ${procs} | ${time_fmt} | ${speedup}x | ${efficiency}% | ${category} |" >> "${REPORT_FILE}"
            fi
        done
    fi

    # OpenMP scaling efficiency
    openmp_baseline=$(cat "${TEMP_DIR}/openmp_baseline.time" 2>/dev/null || echo "N/A")
    if [ "${openmp_baseline}" != "N/A" ] && [ -n "${openmp_baseline}" ]; then
        for procs in 1 2 4 8; do
            if [ ${procs} -eq 1 ]; then
                test_name="openmp_cpu"
            else
                test_name="mpi_${procs}proc_openmp"
            fi

            time_raw=$(cat "${TEMP_DIR}/${test_name}.time" 2>/dev/null || echo "N/A")
            if [ "${time_raw}" != "N/A" ] && [ -n "${time_raw}" ]; then
                if [ ${procs} -eq 1 ]; then
                    speedup="1.00"
                    efficiency="100.0"
                    category="Baseline"
                else
                    speedup=$(echo "scale=2; ${openmp_baseline} / ${time_raw}" | bc -l)
                    efficiency=$(echo "scale=1; (${speedup} / ${procs}) * 100" | bc -l)

                    # Categorize scaling
                    if (( $(echo "${efficiency} >= 90" | bc -l) )); then
                        category="Excellent (≥90%)"
                    elif (( $(echo "${efficiency} >= 75" | bc -l) )); then
                        category="Good (75-90%)"
                    elif (( $(echo "${efficiency} >= 50" | bc -l) )); then
                        category="Moderate (50-75%)"
                    else
                        category="Poor (<50%)"
                    fi
                fi

                time_fmt=$(format_time "${time_raw}")
                echo "| OpenMP/CPU | ${procs} | ${time_fmt} | ${speedup}x | ${efficiency}% | ${category} |" >> "${REPORT_FILE}"
            fi
        done
    fi

    cat >> "${REPORT_FILE}" << 'EOF'

### Performance Metrics Definitions

#### Simulation Time
The wall-clock time (in seconds, minutes, or hours) required to complete the Monte Carlo simulation and LEF dynamics calculations. This excludes file I/O and merging operations.

#### Speedup
Measures how much faster the parallel execution is compared to the baseline single-process run:

```
Speedup = Baseline Simulation Time / Parallel Simulation Time
```

- **Speedup = 1.0x**: No improvement (same speed as baseline)
- **Speedup = 2.0x**: Twice as fast as baseline
- **Speedup < 1.0x**: Slower than baseline (negative scaling due to overhead)

#### Parallel Efficiency
Measures how effectively additional processors are being utilized:

```
Parallel Efficiency = (Speedup / Number of Processes) × 100%
```

Where **Number of Processes** is the MPI process count (1, 2, 4, or 8).

**Example Calculation:**
- Baseline (1 proc): 0.116 s
- MPI 2 procs: 0.115 s
- Speedup = 0.116 / 0.115 = 1.00x
- Efficiency = (1.00 / 2) × 100% = 50.0%

#### Efficiency Categories

- **Excellent (≥90%)**: Nearly linear scaling - adding processors provides proportional speedup
- **Good (75-90%)**: Strong scaling - efficient use of additional processors
- **Moderate (50-75%)**: Acceptable scaling - some overhead but still beneficial
- **Poor (<50%)**: Weak scaling - overhead dominates, consider fewer processors

**Note:** For single-trajectory simulations (Niter=1), poor MPI scaling is expected since work cannot be effectively divided. MPI scaling improves significantly with multiple trajectories (Niter ≥ number of processes).

EOF
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
- **Method:** Parallelization method used (OpenACC/GPU or OpenMP/CPU)
- **Threads:** Number of OpenMP threads (N/A for GPU tests, "auto" for system default)

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

**Report generated:**
EOF

echo "$(date)" >> "${REPORT_FILE}"
echo "" >> "${REPORT_FILE}"

echo "✓ Report generated successfully: ${REPORT_FILE}"
echo ""
echo "To view the report:"
echo "  cat ${REPORT_FILE}"
echo "  # or"
echo "  open ${REPORT_FILE}  # on macOS"
