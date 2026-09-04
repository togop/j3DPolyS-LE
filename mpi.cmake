# --- verify MPI actually works
find_package(MPI COMPONENTS Fortran)

if(NOT MPI_Fortran_FOUND)
  message(FATAL_ERROR "No Fortran MPI found. Install an MPI built with ${CMAKE_Fortran_COMPILER_ID} "
    "or point CMake at one with -DMPI_Fortran_COMPILER=<path to mpif90>.")
endif()

# The `mpi` module is a compiler-specific binary artefact: an MPI whose modules were
# built by another compiler (or another major version of the same compiler) is found by
# FindMPI but cannot be used. Compile a probe and keep the output so the reason is visible.
set(_mpi_probe_dir ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mpi_probe)
file(WRITE ${_mpi_probe_dir}/probe.f90 "use mpi\nend\n")

try_compile(hasMPI ${_mpi_probe_dir}/build
  SOURCES ${_mpi_probe_dir}/probe.f90
  CMAKE_FLAGS
    "-DINCLUDE_DIRECTORIES=${MPI_Fortran_INCLUDE_DIRS}"
    "-DCMAKE_Fortran_FLAGS=${MPI_Fortran_COMPILE_OPTIONS}"
  LINK_LIBRARIES ${MPI_Fortran_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}
  OUTPUT_VARIABLE _mpi_probe_output)

if(NOT hasMPI)
  message(FATAL_ERROR "MPI library not usable with
          ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} (${CMAKE_Fortran_COMPILER})
          Libs: ${MPI_Fortran_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}
          Include: ${MPI_Fortran_INCLUDE_DIRS}
          Opts: ${MPI_Fortran_COMPILE_OPTIONS} ${MPI_Fortran_LINK_FLAGS}

Compiler output of the `use mpi` probe:
${_mpi_probe_output}
Hint: this usually means the detected MPI was built by a different compiler than
${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}. Either select a matching MPI, e.g.
  cmake -Bcmake-build -S . -DMPI_Fortran_COMPILER=/usr/bin/mpif90
or build with the compiler that MPI was built with, e.g.
  cmake -Bcmake-build -S . -DCMAKE_Fortran_COMPILER=nvfortran")
endif()

# --- end verify MPI
