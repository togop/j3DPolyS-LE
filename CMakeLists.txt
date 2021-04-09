cmake_minimum_required(VERSION 3.13)
#cmake_policy(VERSION 3.5)
cmake_policy(SET CMP0068 NEW)
project(3dpolys--le VERSION 2021.1 LANGUAGES Fortran)
enable_language(Fortran C)

#SET(CMAKE_C_COMPILER   /usr/local/bin/gcc-9)
#SET(CMAKE_CXX_COMPILER /usr/local/bin/g++-9)

# set default build paths
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

include(mpi.cmake)
if(NOT hasMPI)
    return()
endif()

set(TARGET 3dpolys-le)
set(PYTHONPATH sim-suite)
#set(PYTHON_FILES sim-suite/__init__.py sim-suite/hic_analysis.py sim-suite/hic_analysis_sim_vs_exp.py sim-suite/plot_hic.py sim-suite/plot_sim_stats_3d.py sim-suite/plot_configs_3d.py)
set(PYTHON_FILES sim-suite/__init__.py sim-suite/hic_analysis.py sim-suite/3dpolys-le_runner.py sim-suite/hic_converters.py sim-suite/job_runner.py)
set(SOURCE_FILES src/logging.f03 src/randomnumber.f03 src/kinds.f03 src/timers.f03 src/lattice_data.f03 src/polymer_model.f03 src/analyse.f03 src/3dpolys-le.f03)


set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f03;f90)
set_source_files_properties(${SOURCE_FILES} PROPERTIES LANGUAGE Fortran)

set(CMAKE_CXX_STANDARD 14)

message( STATUS "PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")
message( STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BUILD_FILES_DIRECTORY})
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

set(CMAKE_Fortran_COMPILER "gfortran")
set(CMAKE_C_COMPILER "gcc")

message( STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER})
message( STATUS "C compiler: " ${CMAKE_C_COMPILER} )
message( STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH} )
message( STATUS "CMAKE_Fortran_MODULE_DIRECTORY: " ${CMAKE_Fortran_MODULE_DIRECTORY} )

find_package(HDF5 REQUIRED COMPONENTS Fortran)

message( STATUS "HDF5_FOUND: " ${HDF5_FOUND} )
message( STATUS "HDF5_Fortran_LIBRARIES: " ${HDF5_Fortran_LIBRARIES})  # /usr/local/Cellar/hdf5/1.10.5_1/lib/
message( STATUS "HDF5_Fortran_HL_LIBRARIES: " ${HDF5_Fortran_HL_LIBRARIES} )
message( STATUS "HDF5_INCLUDE_DIRS: " ${HDF5_INCLUDE_DIRS} )
message( STATUS "HDF5_Fortran_INCLUDE_DIRS: " ${HDF5_Fortran_INCLUDE_DIRS} )
message( STATUS "HDF5_LIBRARIES: " ${HDF5_LIBRARIES} )
message( STATUS "HDF5_Fortran_LIBRARIES: " ${HDF5_Fortran_LIBRARIES} )
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
    set(CMAKE_Fortran_FLAGS "-O3 -fimplicit-none") #  -lhdf5_fortran -fintrinsic-modules-path ${HDF5_INCLUDE_DIRS}") # /usr/local/include") #  TODO find metter way to include ${HDF5_INCLUDE_DIRS}")  #/usr/local/include")  # -lhdf5hl_fortran")          # -std=f2003 ${CMAKE_Fortran_FLAGS}  -qopenmp
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -Wall -fbacktrace -fbounds-check")  #   -fcheck=bounds
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3  -Wall")

    set(CMAKE_C_FLAGS "-c")
endif()

if (${CMAKE_BUILD_TYPE} MATCHES "Release")
    message( STATUS "Fortran compiler flags:" ${CMAKE_Fortran_FLAGS_RELEASE})
elseif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    message( STATUS "Fortran compiler flags:" ${CMAKE_Fortran_FLAGS_DEBUG})
endif()

# -------------------------------------------------------------------------------------------------------------------------
# Libraries
# -------------------------------------------------------------------------------------------------------------------------

# fotrasn objs
SET(RANDOMNUMBER randomnumber)
ADD_LIBRARY(${RANDOMNUMBER} STATIC src/randomnumber.f03)
target_link_libraries(${RANDOMNUMBER})

SET(LOGGING logging)
ADD_LIBRARY(${LOGGING} STATIC src/logging.f03)
target_link_libraries(${LOGGING})

SET(KINDS kinds)
ADD_LIBRARY(${KINDS} STATIC src/kinds.f03)
target_link_libraries(${KINDS})

SET(TIMERS timers)
ADD_LIBRARY(${TIMERS} STATIC src/timers.f03)
target_link_libraries(${TIMERS} ${KINDS})

SET(LATTICE_DATA lattice_data)
ADD_LIBRARY(${LATTICE_DATA} STATIC src/lattice_data.f03)
target_link_libraries(${LATTICE_DATA})

SET(POLYMER_MODEL polymer_model)
ADD_LIBRARY(${POLYMER_MODEL} STATIC src/polymer_model.f03)
target_link_libraries(${POLYMER_MODEL} ${LOGGING} ${RANDOMNUMBER} ${TIMERS} ${KINDS})

SET(ANALYSE analyse)
ADD_LIBRARY(${ANALYSE} STATIC src/analyse.f03)
target_include_directories(${ANALYSE} PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(${ANALYSE} PRIVATE ${LATTICE_DATA} ${POLYMER_MODEL} ${RANDOMNUMBER} ${TIMERS} ${KINDS} ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})


# include_directories(src)

add_executable(${TARGET} ${SOURCE_FILES})
##target_include_directories(${TARGET} PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(${TARGET} ${RANDOMNUMBER} ${KINDS} ${TIMERS} ${LATTICE_DATA} ${POLYMER_MODEL} ${ANALYSE} MPI::MPI_Fortran) ## ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})
set_target_properties(${TARGET} PROPERTIES LINKER_LANGUAGE Fortran)

# data files
#configure_file(src/input_100.dat ${EXECUTABLE_OUTPUT_PATH}/input.dat COPYONLY)
#configure_file(src/data/no-boundary.csv ${EXECUTABLE_OUTPUT_PATH}/no-boundary.csv COPYONLY)
#configure_file(src/data/rex-sites.csv ${EXECUTABLE_OUTPUT_PATH}/rex-sites.csv COPYONLY)
#configure_file(src/data/mex-sites.csv ${EXECUTABLE_OUTPUT_PATH}/mex-sites.csv COPYONLY)
#configure_file(src/data/mex-norex-sites.csv ${EXECUTABLE_OUTPUT_PATH}/mex-norex-sites.csv COPYONLY)
#configure_file(src/data/strong_mex-sites.csv ${EXECUTABLE_OUTPUT_PATH}/strong_mex-sites.csv COPYONLY)
#configure_file(src/data/sim_chromosome_sizes.tsv ${EXECUTABLE_OUTPUT_PATH}/sim_chromosome_sizes.tsv COPYONLY)
#configure_file(src/data/chromosome_sizes.tsv ${EXECUTABLE_OUTPUT_PATH}/chromosome_sizes.tsv COPYONLY)
configure_file(sim-suite/mpirun_3dpolys-le.sh ${EXECUTABLE_OUTPUT_PATH}/mpirun_3dpolys-le.sh COPYONLY)
configure_file(sim-suite/run_3dpolys-le.sh ${EXECUTABLE_OUTPUT_PATH}/run_3dpolys-le.sh COPYONLY)
configure_file(sim-suite/run_3dpolys-le_analysis.sh ${EXECUTABLE_OUTPUT_PATH}/run_3dpolys-le_analysis.sh COPYONLY)
configure_file(sim-suite/run_3dpolys-le_runner.sh ${EXECUTABLE_OUTPUT_PATH}/run_3dpolys-le_runner.sh COPYONLY)
#configure_file(sim-suite/run_plot_chi2-min_time_hist.sh ${EXECUTABLE_OUTPUT_PATH}/run_plot_chi2-min_time_hist.sh COPYONLY)
#configure_file(sim-suite/run_plot_hic.sh ${EXECUTABLE_OUTPUT_PATH}/run_plot_hic.sh COPYONLY)
#configure_file(sim-suite/run_plot_hic3d.sh ${EXECUTABLE_OUTPUT_PATH}/run_plot_hic3d.sh COPYONLY)
configure_file(sim-suite/run_3dpolys-le_stats.sh ${EXECUTABLE_OUTPUT_PATH}/run_3dpolys-le_stats.sh COPYONLY)
configure_file(sim-suite/__init__.py ${EXECUTABLE_OUTPUT_PATH}/__init__.py COPYONLY)
configure_file(sim-suite/hic_analysis.py ${EXECUTABLE_OUTPUT_PATH}/hic_analysis.py COPYONLY)
configure_file(sim-suite/plot_hic.py ${EXECUTABLE_OUTPUT_PATH}/plot_hic.py COPYONLY)
#configure_file(sim-suite/plot_sim_stats_3d.py ${EXECUTABLE_OUTPUT_PATH}/plot_sim_stats_3d.py COPYONLY)
#configure_file(sim-suite/plot_configs_3d.py ${EXECUTABLE_OUTPUT_PATH}/plot_configs_3d.py COPYONLY)
#configure_file(sim-suite/plot_hic3d.py ${EXECUTABLE_OUTPUT_PATH}/plot_hic3d.py COPYONLY)
#configure_file(sim-suite/dcc_grid_sim.py ${EXECUTABLE_OUTPUT_PATH}/dcc_grid_sim.py COPYONLY)
#configure_file(sim-suite/dcc_grid_sim_perm_dir.py ${EXECUTABLE_OUTPUT_PATH}/dcc_grid_sim_perm_dir.py COPYONLY)
configure_file(sim-suite/3dpolys-le_stats.py ${EXECUTABLE_OUTPUT_PATH}/3dpolys-le_stats.py COPYONLY)
#configure_file(sim-suite/dcc_sim_stats_update.py ${EXECUTABLE_OUTPUT_PATH}/dcc_sim_stats_update.py COPYONLY)
configure_file(sim-suite/job_runner.py ${EXECUTABLE_OUTPUT_PATH}/job_runner.py COPYONLY)
configure_file(sim-suite/3dpolys-le_runner.py ${EXECUTABLE_OUTPUT_PATH}/3dpolys-le_runner.py COPYONLY)
#configure_file(sim-suite/plot_chi2-min_time_hist.py ${EXECUTABLE_OUTPUT_PATH}/plot_chi2-min_time_hist.py COPYONLY)
