# Installing j3DPolySLE as a Julia Package

This guide explains how to install j3DPolySLE as a proper Julia package that can be used from anywhere.

## Quick Start

### From Local Directory

```julia
using Pkg
Pkg.develop(path="/path/to/3DPolyS-LE")
```

### From Git Repository

```julia
using Pkg
Pkg.add(url="https://gitlab.com/togop/3DPolyS-LE.git")
```

## Installation Methods

### Method 1: Development Installation (Recommended for Development)

This creates a symlink to your local copy, so changes are immediately available:

```julia
using Pkg
Pkg.develop(path="/path/to/3DPolyS-LE")
```

**Advantages:**
- Changes to source code are immediately available
- Easy to update and develop
- No need to reinstall after changes

### Method 2: Production Installation

This copies the package to Julia's package directory:

```julia
using Pkg
Pkg.add(path="/path/to/3DPolyS-LE")
```

**Advantages:**
- Stable installation
- Independent of source location
- Good for production use

### Method 3: Install from Git

If the package is in a Git repository:

```julia
using Pkg
Pkg.add(url="https://gitlab.com/togop/3DPolyS-LE.git")
```

Or for a specific branch/tag:

```julia
Pkg.add(url="https://gitlab.com/togop/3DPolyS-LE.git", rev="main")
```

## Using the Installed Package

Once installed, you can use it from any Julia session:

```julia
using j3DPolySLE

# Use the modules
j3DPolySLE.Runner3dpolysLe.main()
j3DPolySLE.PlotHic.main()
# etc.
```

## CLI Scripts After Installation

The CLI scripts in `bin/` will automatically detect if the package is installed and use it. If not installed, they fall back to activating the local project.

To make CLI scripts available system-wide:

1. **Add bin/ to PATH:**
   ```bash
   export PATH="$PATH:/path/to/3DPolyS-LE/bin"
   ```

2. **Or create symlinks:**
   ```bash
   ln -s /path/to/3DPolyS-LE/bin/* ~/bin/
   ```

3. **Or use Julia's package system:**
   After installation, you can create wrapper scripts that use the installed package.

## Package Structure

The package follows Julia's standard structure:

```
3DPolyS-LE/
├── Project.toml          # Package metadata and dependencies
├── Manifest.toml         # Locked dependency versions (gitignored)
├── src/
│   └── j3DPolySLE.jl     # Main package file
├── test/
│   └── runtests.jl       # Test suite
└── bin/                   # CLI scripts (optional)
```

## Dependencies

All dependencies are automatically installed when you install the package:

- ArgParse.jl - CLI argument parsing
- CSV.jl - CSV file handling
- DataFrames.jl - Data manipulation
- HDF5.jl - HDF5 file support
- Plots.jl - Plotting
- Statistics.jl - Statistical functions
- Distributions.jl - Probability distributions
- TOML.jl - Configuration file parsing
- PyCall.jl - Python interop (for cooler library)
- Glob.jl - File globbing

## Updating the Package

### Development Installation

Just pull the latest changes:
```bash
cd /path/to/3DPolyS-LE
git pull
```

### Production Installation

Reinstall:
```julia
using Pkg
Pkg.rm("j3DPolySLE")
Pkg.add(path="/path/to/3DPolyS-LE")
```

## Uninstallation

```julia
using Pkg
Pkg.rm("j3DPolySLE")
```

## Troubleshooting

### Package Not Found

If you get `LoadError: ArgumentError: Package j3DPolySLE not found`:

1. Check that you've installed it: `Pkg.status()`
2. Try reinstalling: `Pkg.add(path="/path/to/3DPolyS-LE")`

### Dependency Issues

If dependencies fail to install:

```julia
using Pkg
Pkg.instantiate()
```

### PyCall Issues

The package uses PyCall for the cooler library. Ensure Python and cooler are set up:

```julia
using Pkg
ENV["PYTHON"] = ""  # Use default Python
Pkg.build("PyCall")
```

Then install cooler in Python:
```bash
pip install cooler
```

## Testing

Run the test suite:

```julia
using Pkg
Pkg.test("j3DPolySLE")
```

Or from the project directory:

```julia
using Pkg
Pkg.activate(".")
Pkg.test()
```

## Publishing to Julia Package Registry (Future)

To publish to the Julia General Registry:

1. Create a GitHub/GitLab repository
2. Register with [JuliaRegistries](https://github.com/JuliaRegistries/General)
3. Follow the [package registration guide](https://julialang.github.io/Pkg.jl/v1/creating-packages/)

For now, the package can be installed directly from the Git repository.


