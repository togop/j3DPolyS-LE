# Command-Line Interface for j3DPolySLE

This document describes the command-line interface (CLI) tools available in the Julia port of 3DPolyS-LE.

## Available Commands

All commands are located in the `bin/` directory and can be executed directly:

### 1. `3dpolys_le_runner`

Main runner with batch commands for running simulations and analysis.

**Usage:**
```bash
./bin/3dpolys_le_runner <run_command> [options]
```

**Run Commands:**
- `grid_nlef_km` - Perform simulations over a range of Nlef and km parameter values
- `new_stats` - Perform comparative statistical analysis on simulations
- `decay_plots` - Produce distance-contact-decay plots
- `chip_seq_plots` - Produce ChIP-seq plots
- `contact_radius_analysis` - Perform contact analysis for new parameters
- `multi_decay_plot` - Produce multi distance-contact-decay plot for a single simulation
- `multi_decay_exps_plot` - Produce multi distance-contact-decay plot for multiple datasets
- `run` - Run a single simulation including extracting Hi-C matrices and statistics

**Example:**
```bash
./bin/3dpolys_le_runner run -i input.cfg -o output_folder --nlef 200 --km 2.7e-3
```

### 2. `3dpolys_le_stats`

Perform statistical analysis on simulation results and compare with experimental data.

**Usage:**
```bash
./bin/3dpolys_le_stats [options]
```

**Options:**
- `-o, --output_folder` - Simulation's output folder
- `-a, --analyse` - Analysis step output folder
- `-b, --boundary` - Boundary file in TSV format
- `-e, --exp_cool` - Experimental cooler file for comparison
- `-f, --stats_file` - Statistics repository file (default: ./sim_stats.csv)
- `-r, --radius_contact` - Contact radius in lattice units
- `--replace` - Replace existing files

**Example:**
```bash
./bin/3dpolys_le_stats -o output_folder -a analyse_folder -e exp_data.cool
```

### 3. `plot_hic`

Plot Hi-C contact maps from simulation or experimental data.

**Usage:**
```bash
./bin/plot_hic [options]
```

**Options:**
- `-o, --output_folder` - Output folder containing Hi-C files (default: .)
- `-w, --hic_wildcard` - Hi-C files wildcard (default: hic*.hdf5)
- `-r, --resolution` - Hi-C data resolution in bp (default: 2000)
- `-b, --balanced` - Use balanced matrices if available
- `-c, --cmap` - Color map (default: hot_r)
- `-f, --plot_format` - Image format: png, tif, svg (default: png)
- `--hic_chrs` - Chromosome synonyms to plot

**Example:**
```bash
./bin/plot_hic -o analyse_folder -r 10000 -c YlOrRd
```

### 4. `plot_sim_stats`

Plot simulation statistics as 3D scatter plots or 2D heatmaps.

**Usage:**
```bash
./bin/plot_sim_stats [options]
```

**Options:**
- `-f, --stats_file` - Statistics repository file (default: ./sim_stats.csv)
- `-o, --output_folder` - Output folder for plots (default: .)
- `-e, --file_extension` - Image format: png, tif, svg (default: png)
- `-z, --z_column` - Column to plot on z-axis (default: chi2_log)
- `-p, --plot_mode` - Plot method: 3d or hmap (default: hmap)
- `-c, --cmap` - Color map (default: YlGnBu_r)
- `-lr, --list_contact_radii` - List of contact radii to plot
- `--list_nlef` - Filter by list of Nlef values

**Example:**
```bash
./bin/plot_sim_stats -f sim_stats.csv -z chi2_log -p hmap -lr 2.84 3.55
```

### 5. `hic_converters`

Convert Hi-C files between different formats (.hdf5, .cool, .mcool).

**Usage:**
```bash
./bin/hic_converters -i <input_file> -o <output_file> [options]
```

**Options:**
- `-i, --input_file` - Input file (.hdf5 or .mat format)
- `-o, --output_file` - Output file (.cool or .mcool format)
- `--chr` - Chromosome name (default: chrS)
- `-r, --resolutions` - List of resolutions for .mcool (default: 2000)

**Example:**
```bash
./bin/hic_converters -i hic_81.hdf5 -o output.cool --chr chrX
```

### 6. `hdf5_to_cooler`

Alias for `hic_converters` - converts HDF5 files to cooler format.

**Usage:**
```bash
./bin/hdf5_to_cooler -i <input.hdf5> -o <output.cool> [options]
```

### 7. `3dpolys_le`

Wrapper script to call the Fortran 3dpolys_le binary (if available).

**Usage:**
```bash
./bin/3dpolys_le [options] [config_file]
```

This script searches for the Fortran binary in common locations and executes it with the provided arguments.

## Installation

To use these commands, ensure:

1. Julia is installed and in your PATH
2. The project dependencies are installed:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

3. The scripts are executable:
   ```bash
   chmod +x bin/*
   ```

## Adding to PATH

To use these commands from anywhere, add the `bin/` directory to your PATH:

```bash
export PATH="$PATH:/path/to/3DPolyS-LE/bin"
```

Or create symlinks in a directory already in your PATH (e.g., `~/bin/` or `/usr/local/bin/`).

## Compatibility

These CLI tools are designed to be compatible with the Python version's command-line interface, maintaining the same argument names and behavior where possible.


