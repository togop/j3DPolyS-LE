# --- verify MPI actually works
find_package(MPI COMPONENTS Fortran)

if(NOT MPI_Fortran_FOUND)
  return()
endif()

set(CMAKE_REQUIRED_FLAGS ${MPI_Fortran_COMPILE_OPTIONS})
set(CMAKE_REQUIRED_INCLUDES ${MPI_Fortran_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_Fortran_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
include(CheckFortranSourceCompiles)

# Intel 2019 MPI doesn't yet have mpi_f08
check_fortran_source_compiles("use mpi; end" hasMPI SRC_EXT F90)

if(NOT hasMPI)
  message(STATUS "MPI library not functioning with
          ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}
          Libs: ${MPI_Fortran_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}
          Include: ${MPI_Fortran_INCLUDE_DIRS}
          Opts: ${MPI_Fortran_COMPILE_OPTIONS} ${MPI_Fortran_LINK_FLAGS}")
endif()

# --- end verify MPI
