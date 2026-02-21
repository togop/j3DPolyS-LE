cmake_minimum_required(VERSION 3.13)
cmake_policy(SET CMP0068 NEW)
project(py3dpolys_le VERSION 2025.3 LANGUAGES Fortran)
enable_language(Fortran C)

# Prevent C++ from taking over link language for Fortran static libraries (CMake 3.22+ with C enabled)
set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES FALSE)

# set default build paths
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PYTHONPATH})
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

include(mpi.cmake)
if(NOT hasMPI)
    return()
endif()

set(TARGET 3dpolys_le)
set(SOURCE_FOLDER sim_src)
# Unified sources compiled into executable only (no static libs from same sources = no duplicate symbol linker errors)
set(SOURCE_FILES ${SOURCE_FOLDER}/logging.f03 ${SOURCE_FOLDER}/randomnumber.f03 ${SOURCE_FOLDER}/repro_rng_mod.f03 ${SOURCE_FOLDER}/kinds.f03 ${SOURCE_FOLDER}/timers.f03 ${SOURCE_FOLDER}/lib_conf.f90 ${SOURCE_FOLDER}/lattice_data.f03 ${SOURCE_FOLDER}/polymer_model.f03 ${SOURCE_FOLDER}/gpu_replica_mod.f03 ${SOURCE_FOLDER}/gpu_flat_mod.f03 ${SOURCE_FOLDER}/analyse.f03 ${SOURCE_FOLDER}/3dpolys_le.f03)
set(PACKAGE py3dpolys_le/bin)

set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f03;f90)
set_source_files_properties(${SOURCE_FILES} PROPERTIES LANGUAGE Fortran)

set(CMAKE_CXX_STANDARD 14)

message( STATUS "PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")
message( STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

#SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/${PACKAGE})
#SET(LIBRARY_OUTPUT_PATH ${CMAKE_BUILD_FILES_DIRECTORY})
# set(CMAKE_BUILD_DIRECTORY ${CMAKE_BUILD_FILES_DIRECTORY})
# set(CMAKE_BUILD_FILES_DIRECTORY ${dir})

# generated export header will be placed here
include_directories(${PROJECT_BINARY_DIR})

message( STATUS "CMAKE_BUILD_DIRECTORY: ${CMAKE_BUILD_DIRECTORY}")
message( STATUS "CMAKE_BUILD_FILES_DIRECTORY: ${CMAKE_BUILD_FILES_DIRECTORY}")
message( STATUS "EXECUTABLE_OUTPUT_PATH: ${EXECUTABLE_OUTPUT_PATH}")
message( STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message( STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message( STATUS "Build type: " ${CMAKE_BUILD_TYPE})

# Default to gfortran/gcc only when not already using MPI wrappers (e.g. on HPC with FC=mpifort)
if(NOT CMAKE_Fortran_COMPILER MATCHES ".*(mpifort|mpif90|mpiifort).*")
  set(CMAKE_Fortran_COMPILER "gfortran")
  set(CMAKE_C_COMPILER "gcc")
endif()

# Optional OpenACC (GPU). OFF by default; set -DENABLE_OPENACC=ON to build with OpenACC when compiler supports it.
option(ENABLE_OPENACC "Build with OpenACC for GPU acceleration" OFF)
if(ENABLE_OPENACC)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(OPENACC_Fortran_FLAGS "-fopenacc")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    set(OPENACC_Fortran_FLAGS "-acc -Minfo=accel")
  else()
    message(WARNING "ENABLE_OPENACC=ON but OpenACC flags unknown for ${CMAKE_Fortran_COMPILER_ID}; trying -fopenacc")
    set(OPENACC_Fortran_FLAGS "-fopenacc")
  endif()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OPENACC_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${OPENACC_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OPENACC_Fortran_FLAGS}")
  # Link step needs -fopenacc so mkoffload/offload linking succeeds (GCC)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OPENACC_Fortran_FLAGS}")
  message(STATUS "OpenACC enabled: ${OPENACC_Fortran_FLAGS}")
else()
  message(STATUS "OpenACC disabled (use -DENABLE_OPENACC=ON to enable)")
endif()

message( STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER})
message( STATUS "C compiler: " ${CMAKE_C_COMPILER} )
message( STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH} )
message( STATUS "CMAKE_Fortran_MODULE_DIRECTORY: " ${CMAKE_Fortran_MODULE_DIRECTORY} )

find_package(OpenMP)

# Use HDF5 from environment when set by module (e.g. on Compute Canada: module load hdf5)
if(DEFINED ENV{HDF5_ROOT} AND NOT HDF5_ROOT)
  set(HDF5_ROOT "$ENV{HDF5_ROOT}")
endif()
if(DEFINED ENV{HDF5_DIR} AND NOT HDF5_DIR)
  set(HDF5_DIR "$ENV{HDF5_DIR}")
endif()
find_package(HDF5 REQUIRED COMPONENTS Fortran)

message( STATUS "HDF5_FOUND: " ${HDF5_FOUND} )
message( STATUS "HDF5_Fortran_LIBRARIES: " ${HDF5_Fortran_LIBRARIES})  # /usr/local/Cellar/hdf5/1.10.5_1/lib/
message( STATUS "HDF5_Fortran_HL_LIBRARIES: " ${HDF5_Fortran_HL_LIBRARIES} )
message( STATUS "HDF5_INCLUDE_DIRS: " ${HDF5_INCLUDE_DIRS} )
message( STATUS "HDF5_Fortran_INCLUDE_DIRS: " ${HDF5_Fortran_INCLUDE_DIRS} )
message( STATUS "HDF5_LIBRARIES: " ${HDF5_LIBRARIES} )
message( STATUS "HDF5_LIBRARY_DIRS: " ${HDF5_LIBRARY_DIRS} )

SET(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}")

# TODO find better way to set HDF5 lib folders

message(STATUS "-- CMAKE_SYSTEM_INFO_FILE: " ${CMAKE_SYSTEM_INFO_FILE})
message(STATUS "-- CMAKE_SYSTEM_NAME:      " ${CMAKE_SYSTEM_NAME})
message(STATUS "-- CMAKE_SYSTEM_PROCESSOR: " ${CMAKE_SYSTEM_PROCESSOR})
message(STATUS "-- CMAKE_SYSTEM:           " ${CMAKE_SYSTEM})
message(STATUS "-- OS_NAME:                " ${OS_NAME})
message(STATUS "-- OS_PLATFORM:            " ${OS_PLATFORM})

# HDF5 support

if(HDF5_FOUND)
    message( STATUS "Set link_directories: " ${HDF5_LIBRARY_DIRS})
    message( STATUS "Set include_directories: " ${HDF5_Fortran_INCLUDE_DIRS})
    include_directories(${HDF5_Fortran_INCLUDE_DIRS})
endif()

if(CMAKE_Fortran_COMPILER MATCHES "gfortran*")
    set(CMAKE_SKIP_BUILD_RPATH            FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH    TRUE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

    #Compiler gfortran Flags
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -fimplicit-none")   #  -lhdf5_fortran -fintrinsic-modules-path ${HDF5_INCLUDE_DIRS}") # /usr/local/include") #  TODO find metter way to include ${HDF5_INCLUDE_DIRS}")  #/usr/local/include")  # -lhdf5hl_fortran")          # -std=f2003 ${CMAKE_Fortran_FLAGS}
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g -O0 -Wall -fbacktrace -fbounds-check")  #   -fcheck=bounds
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3  -Wall -static")

    set(CMAKE_C_FLAGS "-c")
endif()

if (${CMAKE_BUILD_TYPE} MATCHES "Release")
    message( STATUS "Fortran compiler flags:" ${CMAKE_Fortran_FLAGS_RELEASE})
elseif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    message( STATUS "Fortran compiler flags:" ${CMAKE_Fortran_FLAGS_DEBUG})
endif()

# -------------------------------------------------------------------------------------------------------------------------
# Executable (all sources in one target to avoid duplicate symbol linker errors)
# -------------------------------------------------------------------------------------------------------------------------

add_executable(${TARGET} ${SOURCE_FILES})
set(LIBS_PRIVATE MPI::MPI_Fortran ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})
if(OpenMP_Fortran_FOUND)
  target_link_libraries(${TARGET} PRIVATE OpenMP::OpenMP_Fortran)
  message(STATUS "OpenMP enabled for CPU replica parallelization")
endif()
target_link_libraries(${TARGET} PRIVATE ${LIBS_PRIVATE})
if(ENABLE_OPENACC)
  target_compile_options(${TARGET} PRIVATE ${OPENACC_Fortran_FLAGS})
  target_link_options(${TARGET} PRIVATE ${OPENACC_Fortran_FLAGS})
  set_target_properties(${TARGET} PROPERTIES LINK_FLAGS "${OPENACC_Fortran_FLAGS}")
endif()
set_target_properties(${TARGET} PROPERTIES LINKER_LANGUAGE Fortran)

# data files
#configure_file(${EXECUTABLE_OUTPUT_PATH}/${TARGET} ${TARGET} COPYONLY)
#configure_file(${SOURCE_FOLDER}/__init__.py ${EXECUTABLE_OUTPUT_PATH}/__init__.py COPYONLY)
#configure_file(${PACKAGE}/run_3dpolys_le.sh ${EXECUTABLE_OUTPUT_PATH}/run_3dpolys_le.sh COPYONLY)
