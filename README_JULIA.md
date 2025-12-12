# 3DPolyS-LE Julia Port

This directory contains the Julia port of the 3DPolyS-LE Python codebase.

## Structure

The Julia code is organized in the `src/` directory with the following modules:

- **j3DPolySLE.jl** - Main module that exports all submodules
- **Logging.jl** - Logging utilities
- **JobRunner.jl** - Job execution and management
- **HicAnalysis.jl** - Hi-C data analysis functions
- **PlotHic.jl** - Hi-C plotting functionality
- **PlotSimStats.jl** - Simulation statistics plotting
- **HicConverters.jl** - Hi-C file format converters
- **Stats3dpolysLe.jl** - Statistics computation module
- **Runner3dpolysLe.jl** - Main runner with CLI interface

## Dependencies

The Julia port requires the following packages (listed in `Project.toml`):

- ArgParse.jl - Command-line argument parsing
- CSV.jl - CSV file reading/writing
- DataFrames.jl - Data manipulation
- HDF5.jl - HDF5 file support
- Plots.jl - Plotting library
- Statistics.jl - Statistical functions
- Distributions.jl - Probability distributions
- TOML.jl - TOML configuration file parsing
- PyCall.jl - Python interop (for cooler library)
- Glob.jl - File globbing

## Installation

### As a Julia Package (Recommended)

Install the package for use from anywhere:

```julia
using Pkg
Pkg.develop(path="/path/to/3DPolyS-LE")  # Development mode
# or
Pkg.add(path="/path/to/3DPolyS-LE")      # Production mode
```

Or install from a Git repository:

```julia
Pkg.add(url="https://gitlab.com/togop/3DPolyS-LE.git", rev ="julia_port")
```

See `README_INSTALL.md` for detailed installation instructions.

### Local Development

For local development:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Usage

### Command-Line Interface

The project provides several command-line tools in the `bin/` directory:

- **`3dpolys_le_runner`** - Main runner with batch commands for simulations
- **`3dpolys_le_stats`** - Statistical analysis and comparison with experimental data
- **`plot_hic`** - Plot Hi-C contact maps
- **`plot_sim_stats`** - Plot simulation statistics as 3D plots or heatmaps
- **`hic_converters`** - Convert Hi-C files between formats
- **`hdf5_to_cooler`** - Convert HDF5 files to cooler format
- **`3dpolys_le`** - Wrapper for the Fortran binary (if available)

All scripts are executable and can be run directly:
```bash
./bin/3dpolys_le_runner run -i input.cfg -o output_folder
```

See `README_CLI.md` for detailed usage information.

**Full user guide (install without clone/build, pure Julia API, batch modes, analysis, plots):** see [JULIA_USER_GUIDE.md](JULIA_USER_GUIDE.md).

### Programmatic Usage

The main entry point is through `Runner3dpolysLe.main()`, which provides the same CLI interface as the Python version.

## Notes

- Some functionality (particularly cooler library operations) uses PyCall to interface with the Python `cooler` library, as there is no native Julia implementation.
- The code maintains compatibility with the original Python command-line interface.
- File locking functionality may need adjustment depending on available Julia packages.

## Differences from Python Version

1. **Module Structure**: Julia uses a different module system, so imports are structured differently
2. **Type System**: Julia's type system is more explicit, requiring type annotations in some places
3. **Error Handling**: Uses Julia's exception handling instead of Python's
4. **DataFrames**: Uses DataFrames.jl instead of pandas, with some API differences
5. **Plotting**: Uses Plots.jl instead of matplotlib, requiring some plot code adjustments

