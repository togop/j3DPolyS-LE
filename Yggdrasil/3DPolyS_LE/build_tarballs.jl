# 3DPolyS-LE Yggdrasil / BinaryBuilder build recipe
# Builds the Fortran executable `3dpolys_le` for macOS and Linux with HDF5 and MPI.
#
# To build locally:
#   1. Clone Yggdrasil: git clone https://github.com/JuliaPackaging/Yggdrasil.git
#   2. Copy this file to Yggdrasil/3/3DPolyS_LE/build_tarballs.jl (or set YGGDRASIL_DIR)
#   3. From Yggdrasil root: julia 3/3DPolyS_LE/build_tarballs.jl --debug
#
# Or set YGGDRASIL_DIR to your Yggdrasil clone and run from this directory:
#   YGGDRASIL_DIR=/path/to/Yggdrasil julia build_tarballs.jl
using BinaryBuilder, Pkg
using Base.BinaryPlatforms

const YGGDRASIL_DIR = get(ENV, "YGGDRASIL_DIR", joinpath(@__DIR__, "..", ".."))
include(joinpath(YGGDRASIL_DIR, "platforms", "mpi.jl"))

name = "j3DPolySLE"
version = v"2026.1.0"

# Repository and ref for reproducible builds (use full commit hash)
repo = "https://gitlab.com/togop/3DPolyS-LE.git"
ref = "83cfb98fc7a4b23e056a81c50d71685fc5228169"

sources = [
    GitSource(repo, ref),
]

# Bash recipe: configure with CMake, build, install executable to prefix/bin
script = raw"""
cd ${WORKSPACE}/srcdir/3DPolyS-LE*

# Remove macOS resource fork files that can cause CMake errors
find /usr/share/cmake -name '._*' -delete 2>/dev/null || true

# Remove hardcoded compiler settings from CMakeLists.txt (lines 49-50)
# These override the cross-compiler toolchain
sed -i.bak '/set(CMAKE_Fortran_COMPILER "gfortran")/d' CMakeLists.txt
sed -i.bak '/set(CMAKE_C_COMPILER "gcc")/d' CMakeLists.txt

# Patch CMakeLists.txt to set Threads variables before HDF5 find_package
# Insert before line 57 (find_package(HDF5...))
sed -i.bak '56 a\
set(CMAKE_THREAD_LIBS_INIT "-lpthread")\
set(CMAKE_HAVE_THREADS_LIBRARY 1)\
set(CMAKE_USE_WIN32_THREADS_INIT 0)\
set(CMAKE_USE_PTHREADS_INIT 1)\
set(Threads_FOUND TRUE)\
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
' CMakeLists.txt

# Remove -static flag which doesn't work on macOS (line 153)
sed -i.bak 's/-static//g' CMakeLists.txt

# Fix duplicate symbol linker errors: build executable from sources only, do not link static libs.
# (sed/awk only to avoid Python site-module encoding issues in the container.)
sed -i.bak 's|\${SOURCE_FOLDER}/timers.f03 \${SOURCE_FOLDER}/lattice_data.f03|\${SOURCE_FOLDER}/timers.f03 \${SOURCE_FOLDER}/lib_conf.f90 \${SOURCE_FOLDER}/lattice_data.f03|' CMakeLists.txt
awk '
BEGIN { lib_block=0; exe_block=0 }
# Start of Libraries block (unique marker in CMakeLists)
/^# Libraries$/ { lib_block=1; next }
lib_block && /^add_executable\(\$\{TARGET\} \$\{SOURCE_FILES\}\)$/ {
  lib_block=0; exe_block=1
  print "# -------------------------------------------------------------------------------------------------------------------------"
  print "# Executable (all sources in one target to avoid duplicate symbol linker errors)"
  print "# -------------------------------------------------------------------------------------------------------------------------"
  print ""
  print "add_executable(\${TARGET} \${SOURCE_FILES})"
  print "target_link_libraries(\${TARGET} PRIVATE MPI::MPI_Fortran \${HDF5_Fortran_LIBRARIES} \${HDF5_Fortran_HL_LIBRARIES})"
  next
}
lib_block { next }
exe_block && /^##target_include_directories/ { next }
exe_block && /^add_dependencies\(\$\{TARGET\}/ { next }
exe_block && /^target_link_libraries\(\$\{TARGET\} PRIVATE \$\{RANDOMNUMBER\}/ { exe_block=0; next }
{ print }
' CMakeLists.txt > CMakeLists.txt.tmp && mv CMakeLists.txt.tmp CMakeLists.txt

# CMake expects HDF5 and MPI at prefix (provided by dependencies)
export HDF5_ROOT=${prefix}
export CMAKE_PREFIX_PATH=${prefix}

# Set compilers to full paths for HDF5 CMake config
export CC=${CC}
export FC=${FC}

# Help CMake find pthread for HDF5 in cross-compilation
export CMAKE_THREAD_LIBS_INIT="-lpthread"
export CMAKE_HAVE_THREADS_LIBRARY=1
export CMAKE_USE_WIN32_THREADS_INIT=0
export CMAKE_USE_PTHREADS_INIT=1
export THREADS_PREFER_PTHREAD_FLAG=ON

cmake -B build -S . \
    -DCMAKE_INSTALL_PREFIX=${prefix} \
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=${CC} \
    -DCMAKE_Fortran_COMPILER=${FC} \
    -DTHREADS_PREFER_PTHREAD_FLAG=ON \
    -DCMAKE_THREAD_LIBS_INIT="-lpthread" \
    -DCMAKE_HAVE_THREADS_LIBRARY=1 \
    -DCMAKE_USE_PTHREADS_INIT=1 \
    -DThreads_FOUND=TRUE \
    -DMPI_HOME=${prefix} \
    -DHDF5_ROOT=${prefix} \
    -DENABLE_OPENMP=ON \
    -DENABLE_OPENACC=OFF

cmake --build build --parallel ${nproc}

# Project does not define install(); executable is in source tree (EXECUTABLE_OUTPUT_PATH)
mkdir -p ${prefix}/bin
cp -f py3dpolys_le/bin/3dpolys_le ${prefix}/bin/3dpolys_le
chmod +x ${prefix}/bin/3dpolys_le

install_license LICENSE
"""

# MPI platform augmentation (same pattern as HDF5/MPICH in Yggdrasil)
augment_platform_block = """
    using Base.BinaryPlatforms
    $(MPI.augment)
    augment_platform!(platform::Platform) = augment_mpi!(platform)
"""

# Supported platforms: Linux and macOS (x86_64 and aarch64), with gfortran ABI variants
platforms = supported_platforms()
filter!(p -> !Sys.iswindows(p), platforms)
platforms = expand_gfortran_versions(platforms)

# Add MPI (MPICH) variants so we get HDF5 built with MPI and MPICH at prefix
platforms, platform_dependencies = MPI.augment_platforms(platforms; MPICH_compat="4.3.0 - 5")

products = [
    ExecutableProduct("3dpolys_le", :three_dpolys_le),
]

dependencies = [
    HostBuildDependency("CMake_jll"),
    Dependency("HDF5_jll"; compat="2.0"),
    Dependency("CompilerSupportLibraries_jll"),
]
append!(dependencies, platform_dependencies)

build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies;
    augment_platform_block,
    julia_compat="1.6",
    preferred_gcc_version=v"12",
)
