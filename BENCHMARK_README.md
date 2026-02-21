# 3DPolyS-LE Parallelization Benchmark Suite

## Overview

This benchmark suite evaluates the performance of different parallelization methods in 3DPolyS-LE, including:
- GPU/OpenACC acceleration
- OpenMP/CPU parallelization
- MPI scaling (1, 2, 4, 8 processes)
- Combined MPI + OpenMP configurations

The suite automatically:
- ✓ Runs simulations with multiple parallelization configurations
- ✓ Uses consistent seed values to ensure reproducibility
- ✓ Verifies result consistency via MD5 checksums
- ✓ Generates detailed performance reports
- ✓ Compares execution times across methods

---

## Quick Start

### Prerequisites

1. **Build the project:**
   ```bash
   mkdir -p cmake-build
   cd cmake-build
   cmake ..
   cmake --build .
   cd ..
   ```

2. **Ensure you have:**
   - Compiled `bin/3dpolys_le` executable
   - MPI installed (`mpirun` available)
   - Test configuration file (default: `./test/test_tads_init_benchmark.cfg`)

### Running the Benchmark

**Option 1: Use defaults (recommended for first run)**
```bash
./benchmark_parallelization.sh
```

**Option 2: Specify custom parameters**
```bash
./benchmark_parallelization.sh <config_file> <seed> <output_dir>
```

**Example:**
```bash
./benchmark_parallelization.sh ./test/test_tads_init_benchmark.cfg 42 ./my_benchmarks
```

### Generating the Report

After the benchmark completes:
```bash
./generate_benchmark_report.sh <benchmark_directory>
```

**Example:**
```bash
./generate_benchmark_report.sh ./benchmark_results/20260207_123456
```

The report will be saved as `BENCHMARK_REPORT.md` in the benchmark directory.

---

## Detailed Usage

### benchmark_parallelization.sh

**Purpose:** Run comprehensive parallelization benchmarks

**Synopsis:**
```bash
./benchmark_parallelization.sh [config_file] [seed] [output_dir]
```

**Arguments:**
| Argument | Description | Default |
|----------|-------------|---------|
| `config_file` | Path to simulation configuration file | `./test/test_tads_init_benchmark.cfg` |
| `seed` | Random seed for reproducibility | `42` |
| `output_dir` | Base output directory for results | `./benchmark_results` |

**What it does:**
1. Creates timestamped output directory: `output_dir/YYYYMMDD_HHMMSS/`
2. Runs 8 different benchmark tests:
   - GPU/OpenACC (1 MPI process)
   - OpenMP/CPU (1 MPI process)
   - MPI 2, 4, 8 processes with GPU
   - MPI 2, 4, 8 processes with OpenMP
3. Saves results and logs for each test
4. Generates checksums for reproducibility verification
5. Creates timing summary

**Output structure:**
```
benchmark_results/
└── YYYYMMDD_HHMMSS/
    ├── logs/
    │   ├── gpu_default.log
    │   ├── openmp_cpu.log
    │   ├── mpi_2proc_gpu.log
    │   ├── mpi_4proc_gpu.log
    │   ├── mpi_8proc_gpu.log
    │   ├── mpi_2proc_openmp.log
    │   ├── mpi_4proc_openmp.log
    │   ├── mpi_8proc_openmp.log
    │   ├── checksums.txt
    │   └── timing_summary.txt
    └── results/
        ├── gpu_default/
        ├── openmp_cpu/
        ├── mpi_2proc_gpu/
        ├── mpi_4proc_gpu/
        ├── mpi_8proc_gpu/
        ├── mpi_2proc_openmp/
        ├── mpi_4proc_openmp/
        └── mpi_8proc_openmp/
```

### generate_benchmark_report.sh

**Purpose:** Generate markdown report from benchmark results

**Synopsis:**
```bash
./generate_benchmark_report.sh <benchmark_directory>
```

**Arguments:**
| Argument | Description | Required |
|----------|-------------|----------|
| `benchmark_directory` | Path to benchmark results directory | Yes |

**What it does:**
1. Extracts timing information from logs
2. Compiles performance comparison table
3. Verifies reproducibility via checksums
4. Generates comprehensive markdown report

**Output:**
- Creates `BENCHMARK_REPORT.md` in the benchmark directory

---

## Understanding the Results

### Performance Metrics

**Simulation Time:**
- Core computation time (Monte Carlo steps + LEF dynamics)
- Excludes I/O and initialization overhead
- Primary metric for comparing parallelization efficiency

**Merge Time:**
- Time spent combining output files from MPI processes
- Increases with number of MPI processes
- Generally small compared to simulation time

**Total Time:**
- Wall-clock time from start to finish
- Includes all overhead (startup, I/O, synchronization)
- Best metric for real-world performance

### Reproducibility Verification

The benchmark suite verifies reproducibility by:
1. Using the same seed value for all runs
2. Computing MD5 checksums of all simulation output files (config.out, dr.out, contact.out, process.out, Nlef.out)
3. Comparing checksums across all test configurations
4. Pairwise checks: **GPU vs OpenMP (no-GPU)**, **MPI vs no-MPI (GPU)**, **MPI vs no-MPI (OpenMP)**

**Expected result:** All checksums should be identical, confirming that:
- Random number generation is deterministic with fixed seed
- GPU (OpenACC) and CPU/OpenMP implementations produce identical results
- MPI and no-MPI runs produce identical merged output for the same seed
- OpenMP (single process with `--no-gpu`, `OMP_NUM_THREADS` set in benchmark) matches GPU and MPI results
- Results are reproducible across execution modes (GPU, OpenMP, MPI+GPU, MPI+no-GPU)

### Interpreting Scaling

**Strong Scaling (Fixed Problem Size):**
- How performance improves with more processes for the same problem
- Ideal: Linear speedup (2x processes = 2x faster)
- Reality: Diminishing returns due to communication overhead

**Efficiency:**
- Efficiency = Speedup / Number of Processes
- 100% = perfect scaling
- >80% = excellent scaling
- <50% = poor scaling (overhead dominates)

**When MPI helps:**
- Multiple trajectories (Niter > number of processes)
- Large problem size (high computation/communication ratio)
- Independent simulations that can run in parallel

**When MPI hurts:**
- Single or few trajectories
- Small problem size
- High communication/computation ratio

---

## Example Workflow

### Complete Benchmark Run

```bash
# 1. Build the project (if not already built)
mkdir -p cmake-build
cd cmake-build
cmake ..
cmake --build .
cd ..

# 2. Run the benchmark suite (this may take a while)
./benchmark_parallelization.sh

# The script will output something like:
# Results location: ./benchmark_results/20260207_210000

# 3. Generate the report
./generate_benchmark_report.sh ./benchmark_results/20260207_210000

# 4. View the report
cat ./benchmark_results/20260207_210000/BENCHMARK_REPORT.md

# Or on macOS:
open ./benchmark_results/20260207_210000/BENCHMARK_REPORT.md
```

### Custom Benchmark with Smaller Test

```bash
# Use a smaller test configuration for faster benchmarking
./benchmark_parallelization.sh \
    ./test/test_tads_shell_input.cfg \
    12345 \
    ./quick_bench

# Generate report
./generate_benchmark_report.sh ./quick_bench/$(ls -t ./quick_bench | head -1)
```

### Testing on Different Hardware

```bash
# On a GPU system
./benchmark_parallelization.sh ./test/test_tads_init_benchmark.cfg 42 ./bench_gpu_server

# On a CPU-only system
./benchmark_parallelization.sh ./test/test_tads_init_benchmark.cfg 42 ./bench_cpu_only

# Compare the reports to see hardware-specific performance
```

---

## Configuration Parameters

### Default Test Configuration

The default benchmark uses `./test/test_tads_init_benchmark.cfg` with:
- **Nchain:** 8860 monomers
- **Niter:** 50 trajectories
- **Nmeas:** 3 measurements
- **Ninter:** 120000 steps between measurements
- **Nlef:** 200 LEFs
- **Boundaries:** Loaded from data file

This represents a realistic production-scale simulation.

### Customizing the Benchmark

**For faster testing** (reduce simulation time):
- Decrease `Niter` (fewer trajectories)
- Decrease `Ninter` (fewer steps between measurements)
- Decrease `Nmeas` (fewer measurement points)

**For comprehensive benchmarking** (better statistics):
- Increase `Niter` (more trajectories for better MPI scaling)
- Increase `Ninter` (longer trajectories, more representative)

Edit the configuration file to adjust these parameters, or create a custom config file.

---

## Troubleshooting

### "Executable not found"

**Problem:** Script can't find `./bin/3dpolys_le`

**Solution:** Build the project first:
```bash
mkdir -p cmake-build
cd cmake-build
cmake ..
cmake --build .
cd ..
```

### "Config file not found"

**Problem:** Specified config file doesn't exist

**Solution:**
- Check the file path: `ls ./test/test_tads_init_benchmark.cfg`
- Use absolute path: `./benchmark_parallelization.sh /full/path/to/config.cfg`
- Create custom config file based on existing examples

### "mpirun: command not found"

**Problem:** MPI is not installed or not in PATH

**Solution:**
- Install MPI (e.g., `brew install open-mpi` on macOS)
- Or use system-specific package manager
- Verify installation: `which mpirun`

### Different checksums detected

**Problem:** Output files have different checksums across runs

**Possible causes:**
1. Seed not being used correctly
2. Non-deterministic code paths
3. Hardware-specific floating point differences
4. Bug in parallelization

**Solution:**
- Verify `--seed` parameter is being passed
- Check logs for seed confirmation
- Compare specific output files to identify differences
- Review recent code changes for non-deterministic operations

### Benchmark takes too long

**Problem:** Benchmark is running for hours

**Solution:**
- Use a smaller configuration file
- Reduce `Niter`, `Nmeas`, or `Ninter` in config
- Run subset of tests by modifying the script
- Use faster hardware (GPU vs CPU)

---

## Advanced Usage

### Running Specific Tests Only

Edit `benchmark_parallelization.sh` and comment out tests you don't need:

```bash
# Comment out tests you want to skip
# run_benchmark "mpi_8proc_gpu" 8 ""
# run_benchmark "mpi_8proc_openmp" 8 "--no-gpu-prefer"
```

### Adding New Test Configurations

Add custom test configurations to the script:

```bash
# Add after existing tests
run_benchmark "my_custom_test" 4 "--no-gpu-prefer --threads:16"
```

### Batch Testing Across Configurations

```bash
#!/bin/bash
# Test multiple configurations
for config in ./test/*.cfg; do
    echo "Testing $config"
    ./benchmark_parallelization.sh "$config" 42 "./batch_results"
done
```

### Automated Performance Regression Testing

```bash
#!/bin/bash
# Run benchmark before and after code changes
./benchmark_parallelization.sh ./test/test_tads_init_benchmark.cfg 42 ./baseline

# Make code changes...

./benchmark_parallelization.sh ./test/test_tads_init_benchmark.cfg 42 ./modified

# Compare reports to detect performance regressions
```

---

## Best Practices

### For Reproducible Research

1. **Always use a fixed seed:**
   ```bash
   ./benchmark_parallelization.sh config.cfg 42 ./results
   ```

2. **Document your hardware:**
   - The report automatically captures system info
   - Keep reports with benchmark results
   - Note any special hardware configurations (GPU model, CPU specs)

3. **Version control:**
   - Commit configuration files with your code
   - Tag releases with benchmark results
   - Track performance changes over time

### For Performance Optimization

1. **Start small:**
   - Use small config for initial testing
   - Identify bottlenecks quickly
   - Scale up once optimized

2. **Test on target hardware:**
   - Benchmark on the system where you'll run production simulations
   - Different hardware may show different optimal configurations

3. **Consider problem size:**
   - Small problems: Prefer single-process OpenMP
   - Large problems: Use MPI with GPU
   - Match MPI processes to trajectory count

### For Continuous Integration

```bash
# Add to CI pipeline
- name: Run performance benchmark
  run: |
    ./benchmark_parallelization.sh ./test/quick_benchmark.cfg 42 ./ci_bench
    ./generate_benchmark_report.sh ./ci_bench/$(ls -t ./ci_bench | head -1)

- name: Upload benchmark results
  uses: actions/upload-artifact@v2
  with:
    name: benchmark-results
    path: ./ci_bench/**/BENCHMARK_REPORT.md
```

---

## File Reference

| File | Purpose |
|------|---------|
| `benchmark_parallelization.sh` | Main benchmark execution script |
| `generate_benchmark_report.sh` | Report generation script |
| `BENCHMARK_README.md` | This documentation file |
| `BENCHMARK_PARALLELIZATION_METHODS.md` | Previous comprehensive benchmark report |

---

## Support

For questions or issues:
1. Check this README for common solutions
2. Review log files in the benchmark output directory
3. Verify your build and environment setup
4. Check for known issues in the project repository

---

## Changelog

**2026-02-07:**
- Initial release of automated benchmark suite
- Support for GPU, OpenMP, and MPI configurations
- Automatic reproducibility verification
- Markdown report generation

---

**Last Updated:** 2026-02-07
