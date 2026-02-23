# GPU Compilation Guide for 3DPolyS-LE

This guide explains how to compile 3DPolyS-LE with GPU support using NVIDIA's HPC SDK (NVHPC compiler).

## Requirements

### Hardware
- **NVIDIA GPU** with compute capability ≥ 6.0 (Pascal or newer recommended)
  - A100 (cc80), V100 (cc70), P100 (cc60), H100 (cc90)
  - Check your GPU: `nvidia-smi --query-gpu=compute_cap --format=csv`

### Software
- **Linux** (Ubuntu 20.04+, RHEL 8+, or equivalent)
- **NVIDIA Drivers** (version ≥ 450.80.02)
- **NVIDIA HPC SDK** (version ≥ 23.1)
- **OpenMPI** compiled with NVHPC support
- **HDF5** (version ≥ 1.10)

### Why Not macOS?
- macOS with Apple Silicon (M1/M2/M3/M4) does NOT support OpenACC
- Apple GPUs use Metal API, not CUDA/OpenACC
- Use CPU-only build on macOS: `make build-gcc-noacc`

## Installation

### 1. Install NVIDIA HPC SDK

#### Via Package Manager (Recommended)
```bash
# Ubuntu/Debian
wget https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK
sudo apt-key add DEB-GPG-KEY-NVIDIA-HPC-SDK
echo "deb [trusted=yes] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /" | \
    sudo tee /etc/apt/sources.list.d/nvhpc.list
sudo apt-get update
sudo apt-get install nvhpc-24-1  # or latest version

# Set environment variables
export NVARCH=`uname -s`_`uname -m`
export NVCOMPILERS=/opt/nvidia/hpc_sdk
export PATH=$NVCOMPILERS/$NVARCH/24.1/compilers/bin:$PATH
export LD_LIBRARY_PATH=$NVCOMPILERS/$NVARCH/24.1/compilers/lib:$LD_LIBRARY_PATH
export MANPATH=$NVCOMPILERS/$NVARCH/24.1/compilers/man:$MANPATH
```

#### On HPC Clusters
```bash
# Check available modules
module avail nvhpc

# Load NVHPC and dependencies
module load nvhpc/23.11
module load openmpi/4.1.4  # or latest compatible version
module load hdf5/1.14.0
```

### 2. Verify Installation
```bash
# Check nvfortran compiler
nvfortran --version

# Check GPU visibility
nvidia-smi

# Check CUDA runtime
nvaccelinfo
```

## Compilation

### Quick Start

**For NVIDIA A100 (most common):**
```bash
make build-nvhpc-a100
```

**For NVIDIA V100:**
```bash
make build-nvhpc-v100
```

**Auto-detect GPU architecture:**
```bash
make build-nvhpc
```

### Advanced Compilation Options

**Custom compute capability:**
```bash
make build-nvhpc GPU_ARCH=cc75
```

**With specific OpenMPI:**
```bash
cmake -Bcmake-build -S . \
    -DENABLE_OPENACC=ON \
    -DCMAKE_Fortran_COMPILER=nvfortran \
    -DCMAKE_C_COMPILER=nvc \
    -DGPU_ARCH=cc80 \
    -DMPI_Fortran_COMPILER=mpifort
cmake --build cmake-build --target 3dpolys_le -j 6
```

**Debug build with GPU support:**
```bash
cmake -Bcmake-build-debug -S . \
    -DENABLE_OPENACC=ON \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_Fortran_COMPILER=nvfortran \
    -DCMAKE_C_COMPILER=nvc \
    -DGPU_ARCH=cc80
cmake --build cmake-build-debug --target 3dpolys_le -j 6
```

## GPU Compute Capabilities

| GPU Model | Compute Capability | Makefile Target |
|-----------|-------------------|-----------------|
| H100      | 9.0 (cc90)        | `build-nvhpc-h100` |
| A100      | 8.0 (cc80)        | `build-nvhpc-a100` |
| V100      | 7.0 (cc70)        | `build-nvhpc-v100` |
| P100      | 6.0 (cc60)        | `build-nvhpc-p100` |
| RTX 3090  | 8.6 (cc86)        | `build-nvhpc GPU_ARCH=cc86` |
| RTX 4090  | 8.9 (cc89)        | `build-nvhpc GPU_ARCH=cc89` |

## Verification

### Check Binary GPU Support
```bash
# Check for CUDA/OpenACC symbols
ldd py3dpolys_le/bin/3dpolys_le | grep -i cuda

# Check OpenACC information (requires NVHPC compiler)
nvdisasm py3dpolys_le/bin/3dpolys_le | grep -i "\.target"
```

### Run Benchmarks
```bash
# GPU benchmark
./benchmark_parallelization.sh --tests=gpu

# Compare GPU vs CPU
./benchmark_parallelization.sh --tests=gpu-vs-no-gpu
```

## Performance Tuning

### Environment Variables

**Set number of OpenACC gangs (GPU thread blocks):**
```bash
export ACC_NUM_GANGS=256
```

**Enable profiling:**
```bash
export ACC_PROFILER=1
export NVCOMPILER_ACC_TIME=1
```

**Force GPU device:**
```bash
export CUDA_VISIBLE_DEVICES=0  # Use GPU 0
```

### Compile-Time Options

The GPU optimizations include:
- **Batched measurements**: Reduces host↔device transfers by 10x
- **Optimized transfer size**: Skips unnecessary array transfers (70% reduction)
- **Async operations**: Overlaps GPU computation with host I/O
- **Gang/vector tuning**: Optimized thread organization (`vector_length(128)`)
- **RNG buffering**: Pre-generates random numbers in batches

Adjust batch size in configuration:
```fortran
! In sim_src/3dpolys_le.f03:882
call do_simulation_replicas_flat(replicas, rank_Niter, Ninter, Nmeas, &
                                  burnin, burnout, burnoutM, batch_size=20)
```

## Troubleshooting

### Compilation Errors

**Error: `nvfortran: command not found`**
```bash
# Add NVHPC to PATH
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.1/compilers/bin:$PATH
```

**Error: `HDF5 not found`**
```bash
# Install HDF5 or set HDF5_ROOT
export HDF5_ROOT=/path/to/hdf5
```

**Error: MPI compiler mismatch**
```bash
# Ensure OpenMPI was compiled with NVHPC
module load openmpi/nvhpc
```

### Runtime Errors

**Error: `CUDA_ERROR_NO_DEVICE`**
```bash
# Check GPU visibility
nvidia-smi
export CUDA_VISIBLE_DEVICES=0
```

**Error: `Out of memory`**
```bash
# Reduce number of replicas (Niter) in config file
# Or use managed memory (already enabled by default)
```

**Slow GPU performance (worse than CPU)**
- Check GPU is actually being used: `nvidia-smi` during run
- Ensure `-acc` flag was used during compilation
- Try increasing `Niter` (number of replicas) for better GPU saturation
- Minimum recommended: `Niter >= 128`

### Verification Tests

**Check if GPU code is running:**
```bash
# This should show GPU activity during run
watch -n 0.5 nvidia-smi

# Check for OpenACC kernels in log
./bin/3dpolys_le --config test/test_tads_init_benchmark.cfg 2>&1 | grep "gpu_flat_mod"
```

**Expected log output:**
```
[INFO] gpu_flat_mod - GPU batch size (measurements per transfer): 10
[INFO] gpu_flat_mod - Replica batch (flat GPU), Measurement: 10 (batched 10 intervals)
```

## Expected Performance

With proper GPU compilation and hardware:

| Configuration | CPU (16 cores) | GPU (A100) | Speedup |
|---------------|----------------|------------|---------|
| Niter=128, Nmeas=10 | ~20 min | ~2-3 min | 7-10x |
| Niter=512, Nmeas=10 | ~80 min | ~6-8 min | 10-13x |

**Note:** GPU performance scales with number of replicas (`Niter`). For optimal GPU utilization, use `Niter >= 256`.

## Support

- **NVIDIA HPC SDK docs**: https://docs.nvidia.com/hpc-sdk/
- **OpenACC specification**: https://www.openacc.org/specification
- **Issue tracker**: https://github.com/yourrepo/3dpolys-le/issues

## Summary

```bash
# Complete workflow for GPU compilation on HPC cluster:
module load nvhpc openmpi hdf5
make build-nvhpc-a100
./benchmark_parallelization.sh --tests=gpu
```

For macOS or systems without NVIDIA GPUs, use CPU-only build:
```bash
make build-gcc-noacc
```
