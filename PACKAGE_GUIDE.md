# Julia Package Installation Guide

This guide explains how to set up j3DPolySLE as an installable Julia package.

## Package Structure

The package follows Julia's standard package structure:

```
3DPolyS-LE/
├── Project.toml          # Package metadata, dependencies, and version
├── Manifest.toml         # Locked dependency versions (auto-generated, gitignored)
├── src/
│   └── j3DPolySLE.jl     # Main package module
│   └── *.jl              # Submodules
├── test/
│   └── runtests.jl       # Test suite
└── bin/                   # CLI executable scripts
```

## Installation Methods

### 1. Development Installation (Recommended for Contributors)

Creates a symlink to your local copy. Changes are immediately available:

```julia
using Pkg
Pkg.develop(path="/path/to/3DPolyS-LE")
```

**After installation:**
```julia
using j3DPolySLE
j3DPolySLE.Runner3dpolysLe.main()
```

### 2. Production Installation

Copies the package to Julia's package directory:

```julia
using Pkg
Pkg.add(path="/path/to/3DPolyS-LE")
```

### 3. Install from Git Repository

```julia
using Pkg
Pkg.add(url="https://gitlab.com/togop/3DPolyS-LE.git")
```

For a specific branch:
```julia
Pkg.add(url="https://gitlab.com/togop/3DPolyS-LE.git", rev="main")
```

## Package Metadata

The `Project.toml` file contains:

- **name**: `j3DPolySLE`
- **uuid**: Unique identifier for the package
- **version**: `2025.2`
- **authors**: Package maintainers
- **description**: Package description
- **license**: MIT
- **repository**: Git repository URL
- **[deps]**: Required dependencies
- **[compat]**: Version compatibility constraints
- **[extras]**: Optional dependencies (e.g., Test)

## Using the Installed Package

Once installed, use it from any Julia session:

```julia
using j3DPolySLE

# Access modules
j3DPolySLE.Runner3dpolysLe.main()
j3DPolySLE.PlotHic.main()
j3DPolySLE.Stats3dpolysLe.main()
```

## CLI Scripts

The CLI scripts in `bin/` automatically detect if the package is installed:

- If installed: Uses the installed package
- If not installed: Falls back to local project activation

Make scripts available system-wide:

```bash
export PATH="$PATH:/path/to/3DPolyS-LE/bin"
```

## Dependencies

All dependencies are automatically installed with the package:

- **ArgParse.jl** - CLI argument parsing
- **CSV.jl** - CSV file handling  
- **DataFrames.jl** - Data manipulation
- **HDF5.jl** - HDF5 file support
- **Plots.jl** - Plotting library
- **Statistics.jl** - Statistical functions
- **Distributions.jl** - Probability distributions
- **TOML.jl** - Configuration parsing
- **PyCall.jl** - Python interop (for cooler)
- **Glob.jl** - File globbing

## Testing

Run the test suite:

```julia
using Pkg
Pkg.test("j3DPolySLE")
```

Or from project directory:

```julia
Pkg.activate(".")
Pkg.test()
```

## Updating

### Development Installation
Just pull latest changes - no reinstall needed.

### Production Installation
```julia
Pkg.rm("j3DPolySLE")
Pkg.add(path="/path/to/3DPolyS-LE")
```

## Uninstallation

```julia
using Pkg
Pkg.rm("j3DPolySLE")
```

## Package Location

After installation, the package is located at:

- **Development**: `~/.julia/dev/j3DPolySLE/` (symlink)
- **Production**: `~/.julia/packages/j3DPolySLE/<version>/`

## Troubleshooting

### Package Not Found
```julia
Pkg.status()  # Check if installed
Pkg.add(path="/path/to/3DPolyS-LE")  # Reinstall
```

### Dependency Issues
```julia
Pkg.instantiate()  # Install all dependencies
```

### PyCall Setup
```julia
ENV["PYTHON"] = ""
Pkg.build("PyCall")
```

Then install cooler:
```bash
pip install cooler
```

## Publishing to Julia Registry

To publish to the Julia General Registry:

1. Push to a public Git repository (GitHub/GitLab)
2. Register at [JuliaRegistries/General](https://github.com/JuliaRegistries/General)
3. Follow the [package registration guide](https://julialang.github.io/Pkg.jl/v1/creating-packages/)

For now, install directly from Git:
```julia
Pkg.add(url="https://gitlab.com/togop/3DPolyS-LE.git")
```


