# Yggdrasil / BinaryBuilder: Deployable Fortran Package

This document describes how to build and use **3DPolyS_LE** as a deployable binary artifact (macOS and Linux) via [Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil) and [BinaryBuilder](https://github.com/JuliaPackaging/BinaryBuilder.jl). The recipe compiles the Fortran simulator `3dpolys_le` with HDF5 and MPI (MPICH).

## Artifacts and platforms

The recipe produces:

- **Executable:** `3dpolys_le` (Fortran binary)
- **Platforms:** Linux (x86_64, aarch64) and macOS (x86_64, aarch64), with MPI (MPICH) and multiple `libgfortran` ABI versions for compatibility.

## Prerequisites

- [Julia](https://julialang.org/) 1.6+
- [BinaryBuilder.jl](https://github.com/JuliaPackaging/BinaryBuilder.jl)
- For **local builds:** a clone of [Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil) (needed for the MPI platform helpers)

```bash
git clone https://github.com/JuliaPackaging/Yggdrasil.git
```

## Building the tarballs

### Option A: Build from the Yggdrasil repo (recommended)

1. Clone Yggdrasil and add the recipe:

   ```bash
   git clone https://github.com/JuliaPackaging/Yggdrasil.git
   cd Yggdrasil
   mkdir -p 3/3DPolyS_LE
   cp /path/to/3DPolyS-LE/Yggdrasil/3DPolyS_LE/build_tarballs.jl 3/3DPolyS_LE/
   ```

2. Build (all platforms, or a subset):

   ```bash
   julia --project=build_tarballs_3DPolyS_LE 3/3DPolyS_LE/build_tarballs.jl
   ```

   For a single platform (e.g. macOS ARM):

   ```bash
   julia 3/3DPolyS_LE/build_tarballs.jl aarch64-apple-darwin
   ```

   Debug build (interactive shell):

   ```bash
   julia 3/3DPolyS_LE/build_tarballs.jl --debug
   ```

### Option B: Build from this repo

Set `YGGDRASIL_DIR` to your Yggdrasil clone so the recipe can load `platforms/mpi.jl`:

```bash
export YGGDRASIL_DIR=/path/to/Yggdrasil
cd /path/to/3DPolyS-LE/Yggdrasil/3DPolyS_LE
julia build_tarballs.jl
```

## Output

Successful builds create tarballs under the current directory, e.g.:

- `3DPolyS_LE.v2026.1.0.x86_64-linux-gnu-libgfortran5-cxx11-mpi+mpich.tar.gz`
- `3DPolyS_LE.v2026.1.0.x86_64-apple-darwin14-libgfortran5-cxx11-mpi+mpich.tar.gz`
- … (and other platform/ABI variants)

Each tarball contains `bin/3dpolys_le`.

## Using the artifacts in Julia

### 1. Local Artifacts.toml (before the JLL is registered)

After building, you can host the tarballs (e.g. on GitHub Releases or a server) and reference them from this repo with an `Artifacts.toml` in the project.

Example layout (hashes and URLs must be filled from your build):

```toml
# Artifacts.toml (example – replace hashes and URLs with your build output)
[3dpolys_le]
arch = "x86_64"
os = "linux"
libc = "glibc"
git-tree-sha1 = "<sha1-from-build>"

[[3dpolys_le.download]]
url = "https://github.com/your-org/3DPolyS-LE/releases/download/v2025.3/3DPolyS_LE.v2026.1.0.x86_64-linux-gnu-libgfortran5-cxx11-mpi+mpich.tar.gz"
sha256 = "<sha256-from-build>"

# Repeat for other platforms: aarch64-linux-gnu, x86_64-apple-darwin, aarch64-apple-darwin, etc.
```

Then in Julia you can use the artifact to get the path to `3dpolys_le`:

```julia
using Artifacts
art = artifact"3dpolys_le"
path_3dpolys_le = joinpath(art, "bin", "3dpolys_le")
```

### 2. After the JLL is registered (General registry)

Once the recipe is merged into Yggdrasil and the JLL is registered (e.g. `3DPolyS_LE_jll`), add it as a dependency of **j3DPolySLE**:

In `Project.toml`:

```toml
[deps]
# ... existing deps ...
# 3DPolyS_LE_jll = "<uuid-from-registry>"   # Fortran binary for run
```

Then in code (e.g. in `Runner3dpolysLe.jl` or `JobRunner.jl`), use the JLL’s executable product when available:

```julia
# Prefer JLL binary if loaded
function get_3dpolys_le_path()
    if isdefined(Main, :ThreeDPolyS_LE_jll) || @isdefined(ThreeDPolyS_LE_jll)
        return ThreeDPolyS_LE_jll.three_dpolys_le()
    end
    # Fallback: system PATH or config
    return "3dpolys_le"
end
```

(Exact module name will match what the JLL’s `src/3DPolyS_LE_jll.jl` exports; it may be `three_dpolys_le` as the product name.)

## Adding artifacts for macOS and Linux

To ship artifacts for **macOS** and **Linux**:

1. **Build** the tarballs for the desired platforms (see above).
2. **Upload** the `.tar.gz` files to a stable URL (e.g. GitHub Releases).
3. **Create or update `Artifacts.toml`** in this repo with one entry per platform, each with:
   - `arch`, `os`, and optionally `libc` / `cxx11` / `libgfortran` to match the tarball.
   - `git-tree-sha1` (from `BinaryBuilder.runtime_tarball_hash(...)` or the build log).
   - `download.url` and `download.sha256` for the tarball.

Example minimal entries:

```toml
# Linux x86_64 (glibc, one ABI)
[3dpolys_le-x86_64-linux]
arch = "x86_64"
os = "linux"
libc = "glibc"
git-tree-sha1 = "..."
[[3dpolys_le-x86_64-linux.download]]
url = "https://..."
sha256 = "..."

# macOS x86_64
[3dpolys_le-x86_64-macos]
arch = "x86_64"
os = "macos"
git-tree-sha1 = "..."
[[3dpolys_le-x86_64-macos.download]]
url = "https://..."
sha256 = "..."

# macOS ARM (M1/M2)
[3dpolys_le-aarch64-macos]
arch = "aarch64"
os = "macos"
git-tree-sha1 = "..."
[[3dpolys_le-aarch64-macos.download]]
url = "https://..."
sha256 = "..."
```

You can add more entries for `aarch64-linux`, different `libgfortran` versions, etc., so that `artifact"3dpolys_le"` (or your chosen name) resolves on all supported Mac and Linux variants.

## Recipe details

- **Source:** GitLab repo `https://gitlab.com/togop/3DPolyS-LE.git`, ref `julia_port` (change `ref` in `build_tarballs.jl` for a tag or branch).
- **Dependencies:** CMake (host), HDF5_jll (with Fortran), MPICH_jll, CompilerSupportLibraries_jll.
- **Build:** CMake with `-DCMAKE_TOOLCHAIN_FILE`, `-DMPI_HOME`, `-DHDF5_ROOT`, OpenMP enabled, OpenACC disabled.
- **Output:** Single executable `bin/3dpolys_le`; no `install()` in CMake, so the script copies it from the project’s `py3dpolys_le/bin/` into `$prefix/bin/`.

## Submitting to Yggdrasil

To have the JLL built and registered automatically:

1. Fork [JuliaPackaging/Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil).
2. Add the recipe under a suitable folder, e.g. `3/3DPolyS_LE/build_tarballs.jl`.
3. Open a PR; CI will build for all platforms. After merge and registration, `3DPolyS_LE_jll` will be available in the General registry.

## Troubleshooting

- **MPI not found:** Ensure you are using the MPI-augmented platforms (the recipe includes `MPI.augment_platforms`). Building without MPI will fail because the project’s CMake requires MPI.
- **HDF5 Fortran not found:** The recipe depends on `HDF5_jll` (with Fortran). If your build is for an MPI platform, the matching HDF5 (MPI-enabled) artifact will be used.
- **`platforms/mpi.jl` not found:** Set `YGGDRASIL_DIR` to the root of your Yggdrasil clone when building from this repo.
