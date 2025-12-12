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

