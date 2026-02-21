# Setting Up Dependencies for j3DPolySLE

## System Dependencies

### HDF5

HDF5 requires system libraries to be installed:

**macOS:**
```bash
brew install hdf5
```

**Linux (Debian/Ubuntu):**
```bash
sudo apt-get install libhdf5-dev
```

**Linux (RHEL/CentOS):**
```bash
sudo yum install hdf5-devel
```

**Linux (Fedora):**
```bash
sudo dnf install hdf5-devel
```

### Python and Cooler (for PyCall)

The package uses PyCall to interface with Python's `cooler` library.

**Install Python cooler:**
```bash
pip install cooler
```

Or via conda:
```bash
conda install -c conda-forge cooler
```

## Building Julia Packages

After installing system dependencies, build the Julia packages:

```julia
using Pkg
Pkg.activate(".")
Pkg.build("HDF5")
Pkg.build("PyCall")
```

Or use the automated setup:

```bash
julia --project=. start_julia.jl
```

The script will attempt to build missing dependencies automatically.

## Verification

Test that everything works:

```julia
using Pkg
Pkg.activate(".")
using HDF5
using PyCall
pyimport("cooler")
println("✓ All dependencies working!")
```

## HPC / SLURM clusters (MPI)

If you see:
```text
Could NOT find MPI_Fortran (missing: MPI_Fortran_WORKS)
Could NOT find MPI (missing: MPI_Fortran_FOUND Fortran)
```
when building with `module load gcc openmpi nvhpc`, the usual cause is a **compiler mismatch**: the default Fortran becomes **nvfortran** (from nvhpc) while OpenMPI was built with **gfortran**, so the MPI Fortran test fails.

**Option 1 – Use MPI wrappers as compilers (recommended)**  
Configure and build using the MPI Fortran/C wrappers so the whole build matches the MPI stack:
```bash
module load gcc openmpi    # optional: omit nvhpc if you are not using OpenACC with NVHPC
# Use MPI wrappers as compilers (set before any cmake run)
export FC=mpifort
export CC=mpicc
make build
```
Or pass them to CMake explicitly:
```bash
cmake -G "CodeBlocks - Unix Makefiles" -Bcmake-build -S . \
  -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc
cmake --build cmake-build --target 3dpolys_le -- -j 6
```

**Option 2 – Build with gcc + OpenMPI only**  
If you do not need NVHPC for this build, load only gcc and openmpi so the default Fortran is gfortran (matches OpenMPI):
```bash
module purge
module load gcc openmpi
make build
```
Then load nvhpc only when you need to run or build with OpenACC and an NVHPC-built OpenMPI.

**Option 3 – Use an OpenMPI built for NVHPC**  
If your cluster provides something like `openmpi/3.x-nvhpc`, load that with nvhpc so the MPI and Fortran compiler match.

**“nvfortran not in PATH” when using mpifort (e.g. Compute Canada)**  
If you see:
```text
The Open MPI wrapper compiler was unable to find the specified compiler nvfortran in your PATH.
```
then the `mpifort` in use was **built for NVHPC** and calls **nvfortran**; either make `nvfortran` available or switch to the GCC toolchain.

**Option A – Use the NVHPC toolchain (for best OpenACC/GPU)**  
Load the compiler **before** the MPI module so `nvfortran` is on `PATH` when `mpifort` runs. In the **same** shell where you run `make build`:
```bash
module purge
module load nvhpc/23.11         # or nvhpc/23.x — must match the nvhpc used to build openmpi
which nvfortran                 # must print a path; if not, try another nvhpc module
module load openmpi/4.1.5       # OpenMPI built with nvhpc (path contains nvhpc23)
which mpifort                   # should be .../nvhpc23/openmpi/4.1.5/bin/mpifort
make build
```
If `which nvfortran` fails after loading nvhpc, try `module avail nvhpc` and load the exact compiler module that provides `nvfortran`.

**Option B – Use the GCC toolchain (no nvfortran needed)**  
Build with gfortran and the OpenMPI built for GCC so the NVHPC wrapper is never used. Load an HDF5 module for the same compiler too (see **HDF5 on HPC** below).
```bash
module purge
module load gcc openmpi hdf5     # do not load nvhpc — use gcc-built openmpi and hdf5
make build-gcc
```
This uses gfortran and the GCC-built OpenMPI; OpenACC will use gfortran’s implementation (less capable than nvfortran for GPU).

## Troubleshooting

### HDF5 Issues

If HDF5 still fails after installing system libraries:

1. Check that HDF5 is in your PATH:
   ```bash
   which h5cc  # Should show path to HDF5 compiler
   ```

2. Set HDF5_DIR environment variable:
   ```bash
   export HDF5_DIR=$(brew --prefix hdf5)  # macOS
   # or
   export HDF5_DIR=/usr  # Linux (adjust as needed)
   ```

**HDF5 on HPC (e.g. Compute Canada / Alliance)**  
CMake looks for HDF5 with **Fortran** support. On clusters, load an HDF5 module that matches your compiler **before** configuring:

```bash
module load gcc openmpi
module load hdf5                # or hdf5-mpi/1.12, or the gcc-built variant
module spider hdf5              # list available HDF5 modules if unsure
make build-gcc
```

If the module sets `HDF5_ROOT` or `HDF5_DIR`, CMake will use it. If HDF5 is still not found, set it explicitly after loading the module:

```bash
module load hdf5
export HDF5_DIR=$EBROOTHDF5     # EasyBuild sets this; or use the path from "module show hdf5"
make build-gcc
```

3. Rebuild HDF5:
   ```julia
   Pkg.build("HDF5", verbose=true)
   ```

### PyCall Issues

If PyCall fails:

1. Set Python path:
   ```julia
   ENV["PYTHON"] = ""  # Use default Python
   Pkg.build("PyCall")
   ```

2. Or specify Python explicitly:
   ```julia
   ENV["PYTHON"] = "/usr/bin/python3"  # Adjust path
   Pkg.build("PyCall")
   ```

3. Verify cooler is installed:
   ```bash
   python -c "import cooler; print(cooler.__version__)"
   ```

