# j3DPolySLE — Julia implementation user guide

This guide explains how to use the **Julia implementation** of 3DPolyS-LE **without cloning the repository or building** the Fortran binary, and how to use the **pure Julia interface**, run **batch simulations**, perform **analysis**, and generate **plots** with the package tools.

---

## 1. Using the Julia implementation without cloning or building

You can run the Julia tools (runner, stats, plotting, converters) from anywhere by installing the package. The **Fortran simulation binary** (`3dpolys_le`) is optional: the Julia package orchestrates runs and analysis; if you do not have the binary, use `cmd_run=stdout` or `cmd_run=file:path` in your config to **emit commands** and run the simulation elsewhere (e.g. on a cluster).

### 1.1 Install from Git (no clone needed on your machine)

In a Julia session:

```julia
using Pkg
Pkg.add(url="https://gitlab.com/togop/3DPolyS-LE.git", rev="julia_port")
```

Or from the default branch:

```julia
Pkg.add(url="https://gitlab.com/togop/3DPolyS-LE.git")
```

This installs **j3DPolySLE** and its dependencies. No repository clone or build step is required for the Julia part.

### 1.2 Install from a local clone (optional)

If you already have a clone:

```julia
using Pkg
Pkg.develop(path="/path/to/3DPolyS-LE")   # development (symlink)
# or
Pkg.add(path="/path/to/3DPolyS-LE")       # copy into Julia’s packages
```

### 1.3 Use the CLI without the repo in your working directory

After installation, call the tools via Julia so the installed package is used:

```bash
julia -e 'using j3DPolySLE; j3DPolySLE.Runner3dpolysLe.main()' -- run -i input.cfg -o ./out
julia -e 'using j3DPolySLE; j3DPolySLE.Stats3dpolysLe.main()' -- -o ./out -a ./out/r2.84 -i input.cfg -f sim_stats.csv
julia -e 'using j3DPolySLE; j3DPolySLE.PlotHic.main()' -- -o ./out/r2.84 -r 10000
julia -e 'using j3DPolySLE; j3DPolySLE.PlotSimStats.main()' -- -f sim_stats.csv -p hmap
julia -e 'using j3DPolySLE; j3DPolySLE.HicConverters.main()' -- -i hic_001.hdf5 -o out.cool --chr chrX
```

Everything after `--` is passed as command-line arguments to the corresponding tool.

To use the `bin/` scripts **without** having the repo as your working directory, either:

- Add the **clone’s** `bin/` to your `PATH` (scripts will still load the **installed** package if found), or  
- Create small wrappers that call `julia -e 'using j3DPolySLE; ...' -- ...` as above.

### 1.4 Dependencies

The package installs: ArgParse, CSV, DataFrames, HDF5, Plots, Statistics, Distributions, TOML, PyCall, Glob, FileIO, LinearAlgebra. For cooler-based comparison and some I/O, **PyCall** is used; ensure a Python with **cooler** available (e.g. `pip install cooler`).

---

## 2. Pure Julia interface (non-CLI)

You can drive the same logic from Julia without using the shell CLI by passing argument vectors or by using the internal API.

### 2.1 Runner: CLI-style arguments as a vector

`Runner3dpolysLe.main` accepts an optional vector of strings (same as CLI):

```julia
using j3DPolySLE

# Same as: 3dpolys_le_runner run -i input.cfg -o ./out -l 200 -m 2.7e-3 -r 2.84 -e exp.cool
j3DPolySLE.Runner3dpolysLe.main([
    "run",
    "--input_cfg", "input.cfg",
    "--output_folder", "./out",
    "--nlef", "200",
    "--km", "2.7e-3",
    "--radius_contact", "2.84",
    "--exp_cool", "exp.cool",
])
```

All CLI options from the runner are supported; use the long names and pass them as alternating `"key", "value"` (or `"flag"` for booleans). Examples:

- **Input/output:** `--input_cfg`, `--output_folder`, `--analysis_folder` (`-a`), `--stats_file` (`-f`), `--cmd_run_file`
- **Simulation:** `--nlef` (`-l`), `--km` (`-m`), `--boundary` (`-b`), `--lef_loading_sites` (`-s`), `--tads_boundary` (`-t`), `--boundary_direction` (`-d`), `--z_loop` (`-z`), `--unidirectional` (`-u`), `--init_mode` (`-n`), `--radius_contact` (`-r`), `--list_contact_radii` (`-k`), `--contact_probability` (`-y`), `--replace`
- **Comparison:** `--exp_cool` (`-e`), `--cmp_chrs`, `--resolution` (`-q`)
- **Batch:** `--nlef_list`, `--km_list`
- **Stats only:** `--stats`, `--all_stats`

Example for **multi_decay_plot** and **grid_nlef_km**:

```julia
# multi_decay_plot
j3DPolySLE.Runner3dpolysLe.main([
    "multi_decay_plot",
    "--output_folder", "./out",
    "--analysis_folder", "./out/r2.84",
    "--input_cfg", "input.cfg",
    "--exp_cool", "exp.cool",
    "--replace",
])

# grid_nlef_km (parameter grid)
j3DPolySLE.Runner3dpolysLe.main([
    "grid_nlef_km",
    "--input_cfg", "input.cfg",
    "--output_folder", "./out_grid",
    "--exp_cool", "exp.cool",
    "--nlef_list", "100", "200", "300",
    "--km_list", "1e-3", "2.7e-3", "5e-3",
    "--list_contact_radii", "2.84", "3.55",
    "--replace",
])
```

### 2.2 Providing all parameters as in the CLI

Parameters come from (1) the config file `input.cfg` and (2) CLI overrides. In pure Julia you can:

- Use a config file and override only what you need via the argument vector (as above), or  
- Omit or minimalize the config and pass everything via the vector.

Config file format (TOML) under `[3dpolys_le]` includes for example:

- **Polymer:** `Nchain`, `L`, `Ea`, `init_mode`
- **Measurements:** `Niter`, `Nmeas`, `Ninter`, `burnin`, `burnout`, `burnoutM`
- **LEFs:** `kb`, `ku`, `km`, `Nlef`, `lef_loading_sites`, `basal_loading_factor`
- **Boundaries:** `boundary`, `boundary_factor`, `boundary_score`, `boundary_direction`
- **Options:** `z_loop`, `unidirectional`
- **Analysis:** `radius_contact`, `chrom`, `cmp_chrs`, `exp_cool`, `tads_boundary`

Runner CLI overrides (and thus the same names in the argument vector) include: `--input_cfg`, `--output_folder`, `--analysis_folder`, `--nlef`, `--km`, `--boundary`, `--lef_loading_sites`, `--tads_boundary`, `--boundary_direction`, `--exp_cool`, `--radius_contact`, `--contact_probability`, `--z_loop`, `--unidirectional`, `--init_mode`, `--cmp_chrs`, `--resolution`, `--stats_file`, `--replace`, `--nlef_list`, `--km_list`, `--list_contact_radii`, `--stats`, `--all_stats`, `--cmd_run_file`, etc.

So “all parameters as in the CLI” means: use an `input.cfg` that matches your run and pass any overrides with the same option names in the string vector.

### 2.3 Other tools (stats, plot_hic, plot_sim_stats, hic_converters)

These modules expose `main()` that read from `ARGS`. So from Julia you can either:

- **Set `ARGS` and call `main()`:**

```julia
using j3DPolySLE

# 3dpolys_le_stats
original_ARGS = copy(ARGS)
ARGS[:] = ["-o", "./out", "-a", "./out/r2.84", "-i", "input.cfg", "-f", "sim_stats.csv", "-e", "exp.cool"]
try
    j3DPolySLE.Stats3dpolysLe.main()
finally
    ARGS[:] = original_ARGS
end

# plot_hic
ARGS[:] = ["-o", "./out/r2.84", "-r", "10000", "-c", "YlOrRd", "-f", "png"]
j3DPolySLE.PlotHic.main()

# plot_sim_stats
ARGS[:] = ["-f", "sim_stats.csv", "-o", ".", "-z", "chi2_log", "-p", "hmap", "-k", "2.84", "3.55"]
j3DPolySLE.PlotSimStats.main()

# hic_converters
ARGS[:] = ["-i", "hic_001.hdf5", "-o", "out.cool", "--chr", "chrX"]
j3DPolySLE.HicConverters.main()
```

- **Call programmatic APIs where available:**
  - **PlotHic:** `PlotHic.run(output_folder; hic_wildcard="hic*.hdf5", resolution=2000, balanced=false, hic_chrs=PlotHic.CHR_X_SYNONYMS, cmap="hot_r", clim=[-2.75, 0], title=PlotHic.DEFAULT_TITLE, plot_format="png")`
  - **PlotSimStats:** `PlotSimStats.main(args_dict)` where `args_dict` is a Dict with keys: `"stats_file"`, `"output_folder"` (default `"."`), `"file_extension"` (default `"png"`), `"z_column"` (default `"chi2_log"`), `"plot_mode"` (`"hmap"` or `"3d"`), `"cmap"`, `"list_nlef"` (optional vector of Int or `nothing`), `"list_contact_radii"` (vector of strings, e.g. `["2.84", "3.55"]`). Example: `j3DPolySLE.PlotSimStats.main(Dict("stats_file" => "sim_stats.csv", "output_folder" => ".", "plot_mode" => "hmap", "z_column" => "chi2_log", "list_contact_radii" => ["2.84", "3.55"]))`.
  - **HicConverters:** `HicConverters.main(input_file, output_file, chr="chrS", resolutions=[2000])` for HDF5 → .cool/.mcool.

---

## 3. Where to find results

- **Simulation output folder** (`-o` / `output_folder`):  
  Contains raw simulation outputs: `config.out`, `contact.out`, `Nlef.out`, `process.out`, `dr.out`, and optionally a generated `cmd.sh` when using a container prefix.

- **Analysis subfolder** (e.g. `output_folder/r2.84` for radius 2.84):  
  Contains: `hic_001.hdf5`, `hic_002.hdf5`, … (Hi-C matrices), `chip_lef.out`, `chip_lef.bedGraph`, `xyzconfig_001.out`, and a `plots/` subfolder with comparison/decay plots when experimental data is provided.

- **Statistics repository:**  
  The file given by `-f` / `--stats_file` (default `./sim_stats.csv`): one row per (simulation, analysis) with columns such as `sim_hic_file`, `sim_out_folder`, `exp_cool`, `resolution`, `boundary`, `boundary_direction`, `tads_boundary`, `input.cfg`, `nlef`, `km`, `radius_contact`, `chi2_log`, `alpha_log`, `chi2_lin`, `alpha_lin`, `chr`.

- **Hi-C plots:**  
  Next to each `hic_*.hdf5`, e.g. `hic_001_hot_r.png` (or other `cmap`/format from config or CLI).

- **Multi-decay plots:**  
  In the **simulation output folder** (not inside analysis subfolders), produced by `multi_decay_plot`.

- **plot_sim_stats output:**  
  In the folder given by `-o` for that command (default `.`): heatmaps or 3D plots named from the stats file and parameters.

---

## 4. Batch simulation modes (3dpolys_le_runner)

All of these can be run via CLI or via `Runner3dpolysLe.main([...])` with the same options.

### 4.1 `run` — single run (simulate + analyse + stats)

Runs one simulation (unless output already exists and `--replace` is not set), then analysis for the contact radii in `--list_contact_radii`, then stats.

```bash
3dpolys_le_runner run -i input.cfg -o ./out --nlef 200 --km 2.7e-3 -r 2.84 -e exp.cool
```

Optional: `-a` analysis folder, `-b` boundary file, `-t` TADs boundary, `-d` boundary_direction, `-z` z_loop, `-u` unidirectional, `-n` init_mode, `-k` list_contact_radii, `-y` contact_probability, `--cmp_chrs`, `--replace`, `--stats`, `--all_stats`, `--cmd_run_file`, etc.

### 4.2 `grid_nlef_km` — parameter grid

Sweeps over lists of `Nlef` and `km`, running the same pipeline (sim + analysis + stats) for each pair. Contact radii are taken from `--list_contact_radii`.

```bash
3dpolys_le_runner grid_nlef_km -i input.cfg -o ./out_grid -e exp.cool \
  --nlef_list 50 100 200 300 --km_list 5.4e-4 2.7e-3 5e-3 \
  --list_contact_radii 2.84 3.55 --replace
```

Defaults (if not overridden): `--nlef_list` and `--km_list` have built-in default vectors; see `--help`.

### 4.3 `multi_decay_plot` — distance–contact decay plots

Builds distance–contact probability decay plots (log and linear) from existing analysis folders and optional experimental cooler.

```bash
3dpolys_le_runner multi_decay_plot -o ./out -a ./out/r2.84 -i input.cfg -e exp.cool
```

If `-a` is omitted, all subfolders of `-o` (except `plots`) are treated as analysis folders. Outputs go in the simulation output folder.

### 4.4 Generating commands without executing (no Fortran binary)

To only **print** or **write** the underlying commands (e.g. to run on another machine):

- In `input.cfg`, under `[job_runner]`: set `cmd_run=stdout` (print) or `cmd_run=file:/path/to/script.sh` (append to file).
- Or pass `--cmd_run_file /path/to/script.sh` to the runner.

Then run the runner as usual; it will emit the `3dpolys_le` / `3dpolys_le_stats` / `multi_decay_plot` commands instead of executing them.

---

## 5. Analysis possible from simulation results

- **Chi² comparison with experimental Hi-C:** For each analysis folder, compare the last `hic_*.hdf5` to an experimental .cool (or .mcool) on specified chromosomes; compute chi² and scaling (alpha) in **log** and **linear** modes. Done by `3dpolys_le_stats` (and internally by the runner after each run).
- **Distance–contact probability decay:** Compare simulation vs experimental contact probability vs genomic distance; chi² in log/linear. Done by `multi_decay_plot`.
- **Statistics repository:** Append one row per analysis to `sim_stats.csv` (or the file set by `-f`), so you can later filter and plot by `nlef`, `km`, `radius_contact`, etc.
- **Hi-C contact maps:** From `hic_*.hdf5` (and optionally .cool/.mcool) via `plot_hic`.
- **Parameter sweeps:** Use `grid_nlef_km` then `plot_sim_stats` on the resulting stats file to get heatmaps or 3D plots of chi² (or other columns) over `nlef` and `km`.
- **Format conversion:** HDF5 → .cool/.mcool via `hic_converters` / `hdf5_to_cooler` for use in other tools.

---

## 6. Plots and diagrams from results

### 6.1 plot_hic — Hi-C contact maps

Reads `hic_*.hdf5` (and optionally .cool/.mcool) from an analysis folder and saves log-scaled heatmaps.

**Usage:**

```bash
plot_hic -o <analyse_folder> [options]
```

**Options:**

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--output_folder` | `-o` | `.` | Folder containing `hic_*.hdf5` (and optionally .cool) |
| `--hic_wildcard` | `-w` | `hic*.hdf5` | Glob for Hi-C files |
| `--resolution` | `-r` | `2000` | Resolution in bp |
| `--balanced` | `-b` | false | Use balanced matrices for .cool/.mcool |
| `--cmap` | `-c` | `hot_r` | Colormap (e.g. hot_r, YlOrRd, viridis) |
| `--clim` | | `[-2.75, 0]` | Color limits (two numbers) |
| `--plot_format` | `-f` | `png` | Output format: png, tif, svg |
| `--title` | | (template) | Title template |
| `--hic_chrs` | | chr synonyms | Chromosome names to plot |

**Examples:**

```bash
plot_hic -o ./out/r2.84
plot_hic -o ./out/r2.84 -r 10000 -c YlOrRd -f png
```

---

### 6.2 plot_sim_stats — Stats over parameters (heatmaps / 3D)

Plots columns from the stats CSV (e.g. `chi2_log`, `chi2_lin`, `alpha_log`, `alpha_lin`) as 2D heatmaps or 3D scatter over `nlef` and `km`, optionally filtered by contact radius and Nlef.

**Usage:**

```bash
plot_sim_stats -f <stats_file> [options]
```

**Options:**

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--stats_file` | `-f` | `./sim_stats.csv` | Statistics CSV from runs |
| `--output_folder` | `-o` | `.` | Where to save plots |
| `--file_extension` | `-e` | `png` | png, tif, svg |
| `--z_column` | `-z` | `chi2_log` | Column for z-axis / heatmap value |
| `--plot_mode` | `-p` | `hmap` | `hmap` (2D heatmap) or `3d` (3D scatter) |
| `--cmap` | `-c` | `YlGnBu_r` | Colormap |
| `--list_nlef` | | (all) | Restrict to these Nlef values |
| `--list_contact_radii` | `-k` | `2.84 3.55` | Restrict to these contact radii |

**Examples:**

```bash
plot_sim_stats -f sim_stats.csv -z chi2_log -p hmap -k 2.84 3.55
plot_sim_stats -f sim_stats.csv -z chi2_lin -p 3d -o ./plots
```

---

### 6.3 hic_converters — HDF5 ↔ cooler

Converts Hi-C matrices between formats. Input: .hdf5 (3DPolyS-LE style). Output: .cool or .mcool.

**Usage:**

```bash
hic_converters -i <input_file> -o <output_file> [options]
```

**Options:**

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--input_file` | `-i` | (required) | Input .hdf5 (or .mat if supported) |
| `--output_file` | `-o` | (required) | Output .cool or .mcool |
| `--chr` | | `chrS` | Chromosome name in cooler |
| `--resolutions` | `-r` | `2000` | Resolutions for .mcool (space-separated list) |

**Examples:**

```bash
hic_converters -i ./out/r2.84/hic_003.hdf5 -o sim.cool --chr chrX
hic_converters -i hic_003.hdf5 -o sim.mcool -r 2000 5000 10000
```

---

### 6.4 hdf5_to_cooler

Alias for the same conversion: `hic_converters` with the same `-i`, `-o`, `--chr`, `-r` options. Example:

```bash
hdf5_to_cooler -i hic_003.hdf5 -o sim.cool --chr chrX
```

---

## 7. Optional running commands summary

| Tool | Purpose | Typical usage |
|------|----------|----------------|
| **3dpolys_le_runner** | Run simulations and analysis pipelines | `3dpolys_le_runner run -i input.cfg -o ./out` |
| | | `3dpolys_le_runner grid_nlef_km -i input.cfg -o ./out -e exp.cool` |
| | | `3dpolys_le_runner multi_decay_plot -o ./out -a ./out/r2.84 -i input.cfg -e exp.cool` |
| **3dpolys_le_stats** | Compute chi² vs experimental data, append to stats CSV, trigger Hi-C plot | `3dpolys_le_stats -o ./out -a ./out/r2.84 -i input.cfg -f sim_stats.csv -e exp.cool` |
| **plot_hic** | Plot Hi-C contact maps from analysis folder | `plot_hic -o ./out/r2.84 -r 10000 -c YlOrRd` |
| **plot_sim_stats** | Plot stats CSV as heatmaps or 3D | `plot_sim_stats -f sim_stats.csv -z chi2_log -p hmap -k 2.84 3.55` |
| **hic_converters** | HDF5 → .cool / .mcool | `hic_converters -i hic.hdf5 -o out.cool --chr chrX` |
| **hdf5_to_cooler** | Same as hic_converters | `hdf5_to_cooler -i hic.hdf5 -o out.cool` |

For full option lists, run:

```bash
3dpolys_le_runner --help
3dpolys_le_stats --help
plot_hic --help
plot_sim_stats --help
hic_converters --help
hdf5_to_cooler --help
```

Or from Julia:

```julia
julia -e 'using j3DPolySLE; j3DPolySLE.Runner3dpolysLe.main(["--help"])'
# etc.
```

---

## 8. Quick reference: from zero to plots (Julia only, no build)

1. **Install:** `Pkg.add(url="https://gitlab.com/togop/3DPolyS-LE.git", rev="julia_port")`
2. **Get a config and (optional) experimental .cool** (e.g. from the repo’s `test/` or your own).
3. **Emit commands** (if you don’t have the Fortran binary): in config set `cmd_run=stdout` or `cmd_run=file:run.sh`, then run the runner; execute the printed/written commands where the binary is available.
4. **Or run locally** (if you have the Fortran binary):  
   `julia -e 'using j3DPolySLE; j3DPolySLE.Runner3dpolysLe.main()' -- run -i input.cfg -o ./out -e exp.cool`
5. **Stats:**  
   `julia -e 'using j3DPolySLE; j3DPolySLE.Stats3dpolysLe.main()' -- -o ./out -a ./out/r2.84 -i input.cfg -f sim_stats.csv -e exp.cool`
6. **Decay plots:**  
   `julia -e 'using j3DPolySLE; j3DPolySLE.Runner3dpolysLe.main()' -- multi_decay_plot -o ./out -a ./out/r2.84 -i input.cfg -e exp.cool`
7. **Hi-C and parameter plots:**  
   `plot_hic -o ./out/r2.84` and `plot_sim_stats -f sim_stats.csv -p hmap`

Results: simulation outputs in `-o`, analysis (including Hi-C and `plots/`) in `-o/r<radius>`, stats in `sim_stats.csv`, multi-decay and heatmaps in the locations described in section 3.
