#!/usr/bin/env bash
################################################################################
# Generate Benchmark Report from Results
################################################################################
#
# Generates a markdown report from benchmark_parallelization.sh results.
# Compares simulation modes: with GPU/OpenACC vs without GPU (--no-gpu).
# Supports MPI and single-process tests.
#
# Usage:
#   ./generate_benchmark_report.sh <benchmark_dir>
#
# Arguments:
#   benchmark_dir - Directory containing benchmark results (e.g. ./benchmark_results/YYYYMMDD_HHMMSS)
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

        # If log shows GPU disabled by user, treat as CPU (no GPU) run
        if grep -q "GPU (OpenACC) disabled by user" "${log_file}"; then
            method="CPU (no GPU)"
        fi

        # Extract thread count information
        threads="auto"
        if grep -q "OpenMP threads set to:" "${log_file}"; then
            threads=$(grep "OpenMP threads set to:" "${log_file}" | head -1 | awk '{print $NF}')
        elif grep -q "OpenMP enabled with default thread count:" "${log_file}"; then
            threads=$(grep "OpenMP enabled with default thread count:" "${log_file}" | head -1 | awk '{print $NF}')
            threads="${threads} (auto)"
        elif [ "${method}" = "OpenACC" ] || [ "${method}" = "GPU/OpenACC" ]; then
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
- **With vs without GPU/OpenACC**: Same binary run with default (GPU when available) vs \`--no-gpu\` (CPU only)
- **No MPI**: Direct execution — with GPU/OpenACC or without GPU (\`--no-gpu\`)
- **Single MPI Process**: With GPU/OpenACC or without GPU (\`--no-gpu\`)
- **Multi-process MPI**: Scaling with N MPI processes (with GPU or without GPU per process)

All tests used the same seed value to verify reproducibility.

---

## System Information

EOF_HEADER

# Detect GPU information - first try from simulation logs
GPU_INFO="Not detected"
GPU_LOG_INFO=""

# Check logs for GPU information
for log_file in "${LOG_DIR}"/*.log; do
    if [ -f "${log_file}" ]; then
        # Try to extract GPU info from log (common patterns)
        if grep -q "GPU" "${log_file}" 2>/dev/null; then
            GPU_LOG_INFO=$(grep -i "GPU\|OpenACC\|device" "${log_file}" | grep -v "no-gpu" | head -3)
        fi
        break
    fi
done

# Try system-level detection
if command -v nvidia-smi &> /dev/null; then
    GPU_INFO=$(nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader 2>/dev/null | head -1 || echo "NVIDIA GPU detected but unable to query")
elif [[ "$(uname -s)" == "Darwin" ]]; then
    GPU_MODEL=$(system_profiler SPDisplaysDataType 2>/dev/null | grep "Chipset Model:" | head -1 | cut -d: -f2 | xargs || echo "Unknown")
    if [ -n "${GPU_MODEL}" ] && [ "${GPU_MODEL}" != "Unknown" ]; then
        GPU_VRAM=$(system_profiler SPDisplaysDataType 2>/dev/null | grep "VRAM" | head -1 | cut -d: -f2 | xargs || echo "Unknown")
        if [ -n "${GPU_VRAM}" ] && [ "${GPU_VRAM}" != "Unknown" ]; then
            GPU_INFO="${GPU_MODEL}, ${GPU_VRAM}"
        else
            GPU_INFO="${GPU_MODEL}"
        fi
    fi
fi

# Add system information
cat >> "${REPORT_FILE}" << EOF
- **Hostname:** ${HOSTNAME}
- **Operating System:** ${OS} ${OS_VERSION}
- **Architecture:** ${ARCH}
- **CPU Cores:** ${NPROC}
- **GPU:** ${GPU_INFO}
- **Test Date:** $(date)
- **Benchmark Directory:** ${BENCH_DIR}

EOF

# Add GPU details from log if available
if [ -n "${GPU_LOG_INFO}" ]; then
    cat >> "${REPORT_FILE}" << EOF
### GPU Information from Simulation Logs

\`\`\`
${GPU_LOG_INFO}
\`\`\`

EOF
fi

cat >> "${REPORT_FILE}" << EOF
---

## Test Configuration

EOF

# Extract config parameters from log
LOG_FILE=""
if [ -f "${LOG_DIR}/gpu_no_mpi.log" ]; then
    LOG_FILE="${LOG_DIR}/gpu_no_mpi.log"
elif [ -f "${LOG_DIR}/no_gpu_no_mpi.log" ]; then
    LOG_FILE="${LOG_DIR}/no_gpu_no_mpi.log"
elif [ -f "${LOG_DIR}/cpu_no_mpi.log" ]; then
    LOG_FILE="${LOG_DIR}/cpu_no_mpi.log"
else
    # Find first available log file
    LOG_FILE=$(find "${LOG_DIR}" -name "*.log" -type f | head -1)
fi

if [ -n "${LOG_FILE}" ]; then
    CONFIG_FILE=$(grep "Load initial parameters" "${LOG_FILE}" | head -1 | awk -F'from ' '{print $2}' | awk '{print $1}')
    SEED=$(grep "Random seed set to:" "${LOG_FILE}" | head -1 | awk '{print $NF}')

    # Extract simulation parameters from log (format: "Nchain=3000" or "km=1.3e-03")
    get_param() {
        local name="$1"
        grep "${name}=" "${LOG_FILE}" 2>/dev/null | head -1 | sed -n "s/.*${name}=\([^[:space:]]*\).*/\1/p"
    }

    NCHAIN=$(get_param "Nchain")
    L=$(get_param "L")
    NITER=$(get_param "Niter")
    NMEAS=$(get_param "Nmeas")
    NINTER=$(get_param "Ninter")
    NLEF=$(get_param "Nlef")
    KB=$(get_param "kb")
    KU=$(get_param "ku")
    KM=$(get_param "km")
    EA=$(get_param "Ea")
    EI=$(get_param "Ei")
    BURNIN=$(get_param "burnin")
    BURNOUT=$(get_param "burnout")
    BURNOUTM=$(get_param "burnoutM")
    INIT_MODE=$(get_param "init_mode")
    Z_LOOP=$(get_param "z_loop")
    UNIDIRECTIONAL=$(get_param "unidirectional")

    cat >> "${REPORT_FILE}" << EOF
**Configuration File:** \`${CONFIG_FILE}\`

**Seed Value:** ${SEED}

### Simulation Parameters

Parameters used for the benchmark (from simulation logs):

| Parameter | Value | Description |
|-----------|-------|-------------|
| Nchain | ${NCHAIN:-—} | Polymer chain length (monomers of 2 kb) |
| L | ${L:-—} | Compartment box size |
| Niter | ${NITER:-—} | Number of independent trajectories |
| Nmeas | ${NMEAS:-—} | Number of measurement snapshots |
| Ninter | ${NINTER:-—} | Monte Carlo steps between snapshots |
| Nlef | ${NLEF:-—} | Max bound LEFs |
| kb | ${KB:-—} | LEF binding rate |
| ku | ${KU:-—} | LEF half-unbinding rate |
| km | ${KM:-—} | LEF movement rate |
| Ea | ${EA:-—} | Extrusion energy |
| Ei | ${EI:-—} | Interaction energy |
| burnin | ${BURNIN:-—} | Steps before introducing LEFs |
| burnout | ${BURNOUT:-—} | Steps after last measurement (LEFs removed) |
| burnoutM | ${BURNOUTM:-—} | Measurements at end with LEFs removed |
| init_mode | ${INIT_MODE:-—} | Initial folding (z=zigzag, h=helices, s=continue) |
| z_loop | ${Z_LOOP:-—} | Allow LEFs to traverse |
| unidirectional | ${UNIDIRECTIONAL:-—} | LEF direction mode |

---

## Benchmark Results

### Performance Summary Table

| Test Configuration | MPI Procs | Parallelization | Threads | Simulation Time | Method |
|-------------------|-----------|-----------------|---------|----------------|--------|
EOF

    # Store timing data to temp files for scaling analysis
    TEMP_DIR=$(mktemp -d)
    trap "rm -rf ${TEMP_DIR}" EXIT

    # Test selection from benchmark run (mpi = use mpi_*_no_gpu as baseline)
    TEST_SELECTION=$(grep "^Test selection:" "${LOG_DIR}/timing_summary.txt" 2>/dev/null | sed 's/^Test selection:[[:space:]]*//' | tr -d '\r' || echo "all")

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

        # Store baselines (prefer no-MPI versions) for GPU vs no-GPU comparison
        if [ "${test}" = "gpu_no_mpi" ]; then
            echo "${sim_time_raw}" > "${TEMP_DIR}/gpu_baseline.time"
        fi

        if [ "${test}" = "no_gpu_no_mpi" ] || [ "${test}" = "cpu_no_mpi" ]; then
            echo "${sim_time_raw}" > "${TEMP_DIR}/cpu_baseline.time"
        fi

        # Format times for display
        sim_time=$(format_time "${sim_time_raw}")

        # Parse test name to extract description and process count
        desc="${test}"
        procs="1"
        para="${method}"

        # Determine parallelization type from method
        if [ "${method}" = "OpenACC" ] || [ "${method}" = "GPU/OpenACC" ]; then
            para="GPU/OpenACC"
        elif [ "${method}" = "OpenMP" ]; then
            para="OpenMP/CPU"
        elif [ "${method}" = "CPU (no GPU)" ]; then
            para="CPU (--no-gpu)"
        fi

        # Extract MPI process count from test name
        if [[ "${test}" =~ mpi_([0-9]+)proc ]]; then
            procs="${BASH_REMATCH[1]}"
        elif [[ "${test}" =~ no_mpi ]]; then
            procs="0"
        fi

        # Create human-readable description
        if [ "${test}" = "gpu_no_mpi" ]; then
            desc="With GPU/OpenACC (no MPI)"
        elif [ "${test}" = "no_gpu_no_mpi" ]; then
            desc="Without GPU (--no-gpu, no MPI)"
        elif [ "${test}" = "cpu_no_mpi" ]; then
            desc="Without GPU (--no-gpu, no MPI)"
        elif [[ "${test}" =~ mpi_([0-9]+)proc_gpu ]]; then
            desc="MPI ${BASH_REMATCH[1]} procs, with GPU/OpenACC"
        elif [[ "${test}" =~ mpi_([0-9]+)proc_no_gpu ]]; then
            desc="MPI ${BASH_REMATCH[1]} procs, without GPU/OpenACC (--no-gpu)"
        fi

        # Display MPI process count (show "None" for non-MPI tests)
        mpi_display="${procs}"
        if [ "${procs}" = "0" ]; then
            mpi_display="None"
        fi

        echo "| ${desc} | ${mpi_display} | ${para} | ${threads} | ${sim_time} | ${method} |" >> "${REPORT_FILE}"
    done

    # Ensure MPI no-GPU row is present if missing (for backwards compatibility with older runs)
    max_procs=""
    for test in ${test_list}; do
        if [[ "${test}" =~ mpi_([0-9]+)proc ]]; then
            max_procs="${BASH_REMATCH[1]}"
            break
        fi
    done
    if [ -n "${max_procs}" ] && ! echo "${test_list}" | grep -q "mpi_${max_procs}proc_no_gpu"; then
        echo "| MPI ${max_procs} procs, without GPU/OpenACC (--no-gpu) | ${max_procs} | — | — | N/A | Not run |" >> "${REPORT_FILE}"
    fi

    # Add GPU vs no-GPU comparison section when both single-process tests exist
    gpu_time=$(cat "${TEMP_DIR}/gpu_no_mpi.time" 2>/dev/null || echo "N/A")
    no_gpu_time=$(cat "${TEMP_DIR}/no_gpu_no_mpi.time" 2>/dev/null || cat "${TEMP_DIR}/cpu_no_mpi.time" 2>/dev/null || echo "N/A")
    if [ -n "${gpu_time}" ] && [ "${gpu_time}" != "N/A" ] && [ -n "${no_gpu_time}" ] && [ "${no_gpu_time}" != "N/A" ]; then
        gpu_vs_speedup=$(echo "scale=2; ${no_gpu_time} / ${gpu_time}" | bc -l 2>/dev/null || echo "N/A")
        cat >> "${REPORT_FILE}" << EOF

---

## GPU vs no-GPU Comparison (single process)

Same binary and config; only difference is \`--no-gpu\` for the CPU-only run.

| Mode | Simulation time | Relative |
|------|-----------------|----------|
| **With GPU/OpenACC** | $(format_time ${gpu_time}) | baseline |
| **Without GPU (--no-gpu)** | $(format_time ${no_gpu_time}) | $(echo "scale=2; ${no_gpu_time} / ${gpu_time}" | bc -l 2>/dev/null || echo "N/A")x slower |

EOF
        if [ "${gpu_vs_speedup}" != "N/A" ] && [ -n "${gpu_vs_speedup}" ]; then
            echo "With GPU is **${gpu_vs_speedup}x** faster than with \`--no-gpu\` (single process)." >> "${REPORT_FILE}"
            echo "" >> "${REPORT_FILE}"
        fi
    fi

    # Add visual performance summary
    cat >> "${REPORT_FILE}" << 'EOF'

---

## Visual Performance Summary

EOF

    # Choose baseline: when only MPI tests were run, use mpi_*_no_gpu as baseline
    cpu_baseline=$(cat "${TEMP_DIR}/cpu_baseline.time" 2>/dev/null || echo "N/A")
    gpu_baseline=$(cat "${TEMP_DIR}/gpu_baseline.time" 2>/dev/null || echo "N/A")
    baseline_time=""
    baseline_label=""
    baseline_skip=""

    if [ "${TEST_SELECTION}" = "mpi" ] && [ -n "${max_procs}" ] && [ -f "${TEMP_DIR}/mpi_${max_procs}proc_no_gpu.time" ]; then
        baseline_time=$(cat "${TEMP_DIR}/mpi_${max_procs}proc_no_gpu.time" 2>/dev/null)
        baseline_label="MPI ${max_procs} procs, without GPU/OpenACC"
        baseline_skip="mpi_${max_procs}proc_no_gpu"
    elif [ "${cpu_baseline}" != "N/A" ] && [ -n "${cpu_baseline}" ]; then
        baseline_time="${cpu_baseline}"
        baseline_label="Without GPU (no MPI)"
        if [ -f "${TEMP_DIR}/no_gpu_no_mpi.time" ]; then
            baseline_skip="no_gpu_no_mpi"
        else
            baseline_skip="cpu_no_mpi"
        fi
    fi

    if [ -n "${baseline_time}" ] && [ "${baseline_time}" != "N/A" ]; then
        echo "### Performance Comparison (Relative to ${baseline_label} baseline)" >> "${REPORT_FILE}"
        echo "" >> "${REPORT_FILE}"
        echo "\`\`\`" >> "${REPORT_FILE}"
        echo "Baseline (${baseline_label}):  $(format_time ${baseline_time})" >> "${REPORT_FILE}"
        echo "" >> "${REPORT_FILE}"

        # Find all test results
        for test_file in "${TEMP_DIR}"/*.time; do
            test_name=$(basename "${test_file}" .time)

            # Skip synthetic baseline files and the baseline test itself
            if [ "${test_name}" = "cpu_baseline" ] || [ "${test_name}" = "gpu_baseline" ]; then
                continue
            fi
            if [ -n "${baseline_skip}" ] && [ "${test_name}" = "${baseline_skip}" ]; then
                continue
            fi

            time_raw=$(cat "${test_file}" 2>/dev/null || echo "N/A")
            if [ "${time_raw}" != "N/A" ] && [ -n "${time_raw}" ]; then
                speedup=$(echo "scale=2; ${baseline_time} / ${time_raw}" | bc -l)
                bar_length=$(echo "${speedup} * 10" | bc -l | awk '{printf "%d", $1}')
                bar=$(printf '█%.0s' $(seq 1 ${bar_length} 2>/dev/null || echo ""))

                # Format test name for display
                display_name="${test_name}"

                printf "%-40s %s  %s (%.2fx)\n" \
                    "${display_name}:" "$(format_time ${time_raw})" "${bar}" ${speedup} >> "${REPORT_FILE}"
            fi
        done

        # If MPI no-GPU test was not run, list it as Not run
        no_gpu_mpi_test=""
        for t in ${test_list}; do
            if [[ "${t}" =~ mpi_([0-9]+)proc ]]; then
                no_gpu_mpi_test="mpi_${BASH_REMATCH[1]}proc_no_gpu"
                break
            fi
        done
        if [ -n "${no_gpu_mpi_test}" ] && [ ! -f "${TEMP_DIR}/${no_gpu_mpi_test}.time" ]; then
            printf "%-40s %s\n" "${no_gpu_mpi_test}:" "Not run" >> "${REPORT_FILE}"
        fi

        echo "\`\`\`" >> "${REPORT_FILE}"
    fi

    # Add analysis notes
    cat >> "${REPORT_FILE}" << 'EOF'

---

## Performance Analysis

### Performance Metrics Definitions

#### Simulation Time
The wall-clock time (in seconds, minutes, or hours) required to complete the Monte Carlo simulation and LEF dynamics calculations. This excludes file I/O and merging operations.

#### Speedup
Measures how much faster each configuration is compared to the baseline:
- When only MPI tests were run, baseline = MPI N procs, without GPU/OpenACC.
- Otherwise baseline = without GPU (no MPI).

```
Speedup = Baseline Time / Configuration Time
```

- **Speedup = 1.0x**: Same speed as baseline
- **Speedup = 2.0x**: Twice as fast as baseline
- **Speedup > 1.0x**: Faster than baseline
- **Speedup < 1.0x**: Slower than baseline

#### Test Modes
- **With GPU/OpenACC**: Simulation uses GPU when the program was built with OpenACC.
- **Without GPU (--no-gpu)**: Simulation runs on CPU only; use for comparison or when no GPU is available.

EOF
fi

# Add reproducibility verification
REPRO_OUTPUT_FILES="config.out dr.out contact.out process.out Nlef.out"

cat >> "${REPORT_FILE}" << 'EOF'

---

## Reproducibility Verification

Reproducibility is checked across all run types (same seed required):
- **GPU (OpenACC)** — single process or MPI
- **OpenMP (no-GPU)** — single process with `--no-gpu`, uses OpenMP for replica parallelism
- **MPI + GPU** — multiple processes with GPU
- **MPI + no-GPU** — multiple processes, CPU/OpenMP per process

When all runs use the same seed, output files (config.out, dr.out, contact.out, process.out, Nlef.out) must be byte-identical across GPU, OpenMP, MPI, and no-MPI.

### MD5 Checksums

EOF

cat "${LOG_DIR}/checksums.txt" >> "${REPORT_FILE}"

cat >> "${REPORT_FILE}" << 'EOF'

### Verification Result

EOF

# Check reproducibility per output file (checksums format: test_name  file  md5)
REPRO_FAIL=0
for outfile in ${REPRO_OUTPUT_FILES}; do
    unique=$(awk -v f="${outfile}" '$2 == f {print $3}' "${LOG_DIR}/checksums.txt" 2>/dev/null | sort -u | wc -l)
    total=$(awk -v f="${outfile}" '$2 == f {print $3}' "${LOG_DIR}/checksums.txt" 2>/dev/null | wc -l)
    if [ "${total}" -eq 0 ]; then
        continue
    fi
    if [ "${unique}" -eq 1 ]; then
        echo "- **${outfile}:** ✓ identical across all runs (${total} runs)" >> "${REPORT_FILE}"
    else
        echo "- **${outfile}:** ⚠ ${unique} distinct checksums (expected 1)" >> "${REPORT_FILE}"
        REPRO_FAIL=1
    fi
done

# Pairwise checks (GPU vs OpenMP, MPI vs no-MPI)
compare_result_dirs_report() {
    local dir1="$1"
    local dir2="$2"
    local label="$3"
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
    if [ "${same}" -eq 1 ]; then
        echo "- **${label}:** ✓ byte-identical" >> "${REPORT_FILE}"
    else
        echo "- **${label}:** ⚠ outputs differ" >> "${REPORT_FILE}"
        return 1
    fi
}

echo "" >> "${REPORT_FILE}"
echo "**Pairwise checks (same seed):**" >> "${REPORT_FILE}"
if [ -f "${RESULTS_DIR}/gpu_no_mpi/config.out" ] && [ -f "${RESULTS_DIR}/no_gpu_no_mpi/config.out" ]; then
    compare_result_dirs_report "${RESULTS_DIR}/gpu_no_mpi" "${RESULTS_DIR}/no_gpu_no_mpi" "GPU vs OpenMP (no-GPU)" || REPRO_FAIL=1
fi
# MPI vs no-MPI (discover MAX_CPUS from existing mpi_* dirs)
for mpi_gpu in "${RESULTS_DIR}"/mpi_*proc_gpu; do
    [ -d "${mpi_gpu}" ] || continue
    base=$(basename "${mpi_gpu}" _gpu)
    nprocs="${base#mpi_}"
    nprocs="${nprocs%proc}"
    if [ -f "${mpi_gpu}/config.out" ] && [ -f "${RESULTS_DIR}/gpu_no_mpi/config.out" ]; then
        compare_result_dirs_report "${mpi_gpu}" "${RESULTS_DIR}/gpu_no_mpi" "MPI (${nprocs} proc) GPU vs no-MPI GPU" || REPRO_FAIL=1
        break
    fi
done
for mpi_cpu in "${RESULTS_DIR}"/mpi_*proc_no_gpu; do
    [ -d "${mpi_cpu}" ] || continue
    base=$(basename "${mpi_cpu}" _no_gpu)
    nprocs="${base#mpi_}"
    nprocs="${nprocs%proc}"
    if [ -f "${mpi_cpu}/config.out" ] && [ -f "${RESULTS_DIR}/no_gpu_no_mpi/config.out" ]; then
        compare_result_dirs_report "${mpi_cpu}" "${RESULTS_DIR}/no_gpu_no_mpi" "MPI (${nprocs} proc) OpenMP vs no-MPI OpenMP" || REPRO_FAIL=1
        break
    fi
done

echo "" >> "${REPORT_FILE}"
if [ ${REPRO_FAIL} -eq 0 ]; then
    cat >> "${REPORT_FILE}" << 'EOF'
✓ **REPRODUCIBILITY VERIFIED**

All test runs produced byte-for-byte identical results for the same seed across:
GPU (OpenACC), OpenMP (no-GPU), MPI+GPU, and MPI+no-GPU. Determinism is confirmed for all execution modes.

EOF
else
    cat >> "${REPORT_FILE}" << 'EOF'
⚠ **WARNING: Different Checksums Detected**

Some output files differ across runs. This may indicate:
- Non-deterministic behavior in the simulation
- Different random number sequences between GPU and CPU/OpenMP paths
- MPI vs no-MPI ordering or seeding differences
- Hardware-specific or compiler-specific variations

Ensure the same seed is used (e.g. \`--seed:42\`) and review the logs for details.

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
- **Method:** With GPU/OpenACC or without GPU (CPU only, e.g. when run with \`--no-gpu\`)
- **Threads:** Shown where applicable (N/A for GPU runs)

### Recommendations

Based on your results:
- Compare single-process runs: with GPU/OpenACC vs without GPU (\`--no-gpu\`) for a direct speedup.
- Compare MPI runs: with GPU vs without GPU at the same process count.
- Use the configuration with best performance for your problem size and hardware.

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
