# CLI Implementation Summary

## Created CLI Entry Points

All command-line interfaces have been implemented matching the Python `pyproject.toml` entry points:

| Python Entry Point | Julia Script | Module Function |
|-------------------|--------------|-----------------|
| `3dpolys_le` | `bin/3dpolys_le` | Wrapper for Fortran binary |
| `3dpolys_le_runner` | `bin/3dpolys_le_runner` | `Runner3dpolysLe.main()` |
| `3dpolys_le_stats` | `bin/3dpolys_le_stats` | `Stats3dpolysLe.main()` |
| `plot_hic` | `bin/plot_hic` | `PlotHic.main()` |
| `plot_sim_stats` | `bin/plot_sim_stats` | `PlotSimStats.main()` |
| `hic_converters` | `bin/hic_converters` | `HicConverters.main()` |
| `hdf5_to_cooler` | `bin/hdf5_to_cooler` | `HicConverters.hdf5_to_cooler()` |

## Implementation Details

### Script Structure

All scripts follow this pattern:
```julia
#!/usr/bin/env julia
# Get the project root directory
project_root = abspath(joinpath(@__DIR__, ".."))
using Pkg
Pkg.activate(project_root)
using j3DPolySLE
# Call the appropriate main function
```

### Main Functions Added

1. **PlotHic.main()** - Added CLI argument parsing using ArgParse.jl
2. **PlotSimStats.main()** - Added CLI argument parsing and fixed function signatures
3. **HicConverters.main()** - Added CLI argument parsing
4. **HicConverters.hdf5_to_cooler()** - Alias function for compatibility

### Dependencies

All modules now import `ArgParse` for command-line argument parsing, matching the Python `argparse` functionality.

## Usage

All scripts are executable and can be run from the command line:

```bash
# From project root
./bin/3dpolys_le_runner run -i input.cfg -o output_folder

# Or if added to PATH
3dpolys_le_runner run -i input.cfg -o output_folder
```

## Compatibility

The CLI maintains compatibility with the Python version's command-line interface:
- Same argument names (with `--` or `-` prefixes)
- Same default values
- Same behavior where possible

## Notes

- The `3dpolys_le` script is a bash wrapper that searches for the Fortran binary
- All other scripts are Julia scripts that activate the project and call the appropriate module's `main()` function
- Scripts automatically activate the project environment, so they work from any directory


