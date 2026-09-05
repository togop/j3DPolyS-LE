# --- select and verify a usable MPI
#
# The Fortran `mpi` module is a compiler-specific binary artefact, so an MPI is only
# usable if its modules were built by the compiler we build with. FindMPI takes the
# first wrapper on PATH, which on machines with the NVIDIA HPC SDK is HPC-X: its
# mpi.mod comes from nvfortran and gfortran rejects it. When building with GNU, skip
# wrappers inside the HPC SDK but otherwise keep PATH order, so a module-loaded MPI
# still wins. An explicit -DMPI_Fortran_COMPILER=... always takes precedence.
if(NOT MPI_Fortran_COMPILER AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  file(TO_CMAKE_PATH "$ENV{PATH}" _mpi_search_dirs)
  list(APPEND _mpi_search_dirs /usr/bin /usr/local/bin
    /usr/lib/x86_64-linux-gnu/openmpi/bin /usr/lib64/openmpi/bin)

  foreach(_dir IN LISTS _mpi_search_dirs)
    if(_dir MATCHES "hpc_sdk|nvhpc|nvidia")
      continue()
    endif()
    foreach(_name mpifort mpif90)
      if(EXISTS "${_dir}/${_name}")
        set(MPI_Fortran_COMPILER "${_dir}/${_name}"
          CACHE FILEPATH "MPI Fortran compiler wrapper")
        break()
      endif()
    endforeach()
    if(MPI_Fortran_COMPILER)
      break()
    endif()
  endforeach()
endif()

find_package(MPI COMPONENTS Fortran)

if(NOT MPI_Fortran_FOUND)
  message(FATAL_ERROR "No Fortran MPI found. Install an MPI built with "
    "${CMAKE_Fortran_COMPILER_ID} (on Debian/Ubuntu: apt-get install libopenmpi-dev) "
    "or point CMake at one with -DMPI_Fortran_COMPILER=<path to mpif90>.")
endif()

# Compile a probe against the selected MPI and keep its output, so a module mismatch
# reports the compiler error instead of failing later or silently skipping the target.
set(_mpi_probe_dir ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mpi_probe)
file(WRITE ${_mpi_probe_dir}/probe.f90 "use mpi\nend\n")
string(REPLACE ";" " " _mpi_probe_flags "${MPI_Fortran_COMPILE_OPTIONS}")

try_compile(hasMPI ${_mpi_probe_dir}/build
  SOURCES ${_mpi_probe_dir}/probe.f90
  CMAKE_FLAGS
    "-DINCLUDE_DIRECTORIES=${MPI_Fortran_INCLUDE_DIRS}"
    "-DCMAKE_Fortran_FLAGS=${_mpi_probe_flags}"
  LINK_LIBRARIES ${MPI_Fortran_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}
  OUTPUT_VARIABLE _mpi_probe_output)

if(NOT hasMPI)
  message(FATAL_ERROR "MPI library not usable with
          ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} (${CMAKE_Fortran_COMPILER})
          Wrapper: ${MPI_Fortran_COMPILER}
          Libs: ${MPI_Fortran_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}
          Include: ${MPI_Fortran_INCLUDE_DIRS}
          Opts: ${MPI_Fortran_COMPILE_OPTIONS} ${MPI_Fortran_LINK_FLAGS}

Compiler output of the `use mpi` probe:
${_mpi_probe_output}
Hint: `mpi.mod` must be built by the compiler used here
(${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}). Either install a matching MPI
  apt-get install libopenmpi-dev   # Debian/Ubuntu, gfortran-built OpenMPI
and/or select its wrapper
  cmake -Bcmake-build -S . -DMPI_Fortran_COMPILER=/usr/bin/mpif90
or build everything with the compiler that MPI was built with, e.g.
  cmake -Bcmake-build -S . -DCMAKE_Fortran_COMPILER=nvfortran
(the latter needs an HDF5 whose Fortran modules come from that compiler too).")
endif()

message(STATUS "MPI Fortran wrapper: ${MPI_Fortran_COMPILER}")

# --- end verify MPI
