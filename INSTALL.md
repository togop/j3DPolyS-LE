# Installation Guide for j3DPolySLE

This guide explains how to install the j3DPolySLE Julia package.

## Prerequisites

- Julia 1.6 or higher
- Git (for installation from repository)

## Installation Methods

### Method 1: Install from Local Directory (Development)

If you have cloned the repository:

```julia
using Pkg
Pkg.develop(path="/path/to/3DPolyS-LE")
```

Or from within the project directory:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

### Method 2: Install from Git Repository

```julia
using Pkg
Pkg.add(url="https://gitlab.com/togop/3DPolyS-LE.git", subdir=".")
```

Note: Replace the URL with your actual repository URL if different.

### Method 3: Install from Local Path (Production)

```julia
using Pkg
Pkg.add(path="/path/to/3DPolyS-LE")
```

## Post-Installation Setup

### 1. Install Dependencies

After installation, all dependencies should be automatically installed. If not:

```julia
using Pkg
Pkg.instantiate()
```

### 2. Set Up PyCall for Cooler Library

The package uses PyCall to interface with Python's `cooler` library. You may need to configure PyCall:

```julia
using Pkg
ENV["PYTHON"] = ""  # Use default Python
Pkg.build("PyCall")
```

Or if you have a specific Python environment:

```julia
ENV["PYTHON"] = "/path/to/python"
Pkg.build("PyCall")
```

Then install the cooler library in Python:

```bash
pip install cooler
```

Or via conda:

```bash
conda install -c conda-forge cooler
```

### 3. Make CLI Scripts Available

The CLI scripts are located in the `bin/` directory. To use them:

**Option A: Add to PATH**

```bash
export PATH="$PATH:/path/to/3DPolyS-LE/bin"
```

Add this to your `~/.bashrc` or `~/.zshrc` for persistence.

**Option B: Create Symlinks**

```bash
ln -s /path/to/3DPolyS-LE/bin/3dpolys_le_runner ~/bin/
ln -s /path/to/3DPolyS-LE/bin/3dpolys_le_stats ~/bin/
# ... etc
```

**Option C: Use Julia's Package Binaries**

After installation, you can call the functions directly from Julia:

```julia
using j3DPolySLE
j3DPolySLE.Runner3dpolysLe.main()
```

## Verification

Test the installation:

```julia
using j3DPolySLE
# Should load without errors
```

Test CLI (if added to PATH):

```bash
3dpolys_le_runner --help
```

## Troubleshooting

### PyCall Issues

If you encounter issues with PyCall:

1. Ensure Python is installed and accessible
2. Rebuild PyCall: `Pkg.build("PyCall")`
3. Verify cooler is installed: `python -c "import cooler; print(cooler.__version__)"`

### Missing Dependencies

If you get import errors:

```julia
using Pkg
Pkg.instantiate()
```

### CLI Scripts Not Found

- Ensure scripts are executable: `chmod +x bin/*`
- Check that the path is correct in the scripts
- Verify Julia is in your PATH: `which julia`

## Uninstallation

To remove the package:

```julia
using Pkg
Pkg.rm("j3DPolySLE")
```

## Development Installation

For development work:

```julia
using Pkg
Pkg.develop(path="/path/to/3DPolyS-LE")
```

This creates a symlink, so changes to the source code are immediately available.


