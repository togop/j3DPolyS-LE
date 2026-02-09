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
set(SOURCE_FILES ${SOURCE_FOLDER}/logging.f03 ${SOURCE_FOLDER}/randomnumber.f03 ${SOURCE_FOLDER}/kinds.f03 ${SOURCE_FOLDER}/timers.f03 ${SOURCE_FOLDER}/lattice_data.f03 ${SOURCE_FOLDER}/polymer_model.f03 ${SOURCE_FOLDER}/analyse.f03 ${SOURCE_FOLDER}/3dpolys_le.f03)
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

# OpenMP CPU parallelization support (works on Linux, macOS, Windows)
option(ENABLE_OPENMP "Enable OpenMP CPU parallelization" OFF)

# Variable to track if we should link OpenMP (set after executable is created)
set(OPENMP_FORTRAN_FOUND FALSE)
set(OPENMP_USE_FLAGS FALSE)

if(ENABLE_OPENMP)
    # Find OpenMP - on Mac, OpenMP_C might not be found but OpenMP_Fortran can be
    # Try to find OpenMP Fortran component specifically
    find_package(OpenMP)
    # Check if Fortran component is found (even if C component is not)
    if(OpenMP_Fortran_FOUND)
        set(OPENMP_FORTRAN_FOUND TRUE)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
        message(STATUS "OpenMP enabled: Using ${OpenMP_Fortran_FLAGS}")
    else()
        # Try to set flags manually if FindOpenMP doesn't work
        set(OPENMP_USE_FLAGS TRUE)
        if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
            message(STATUS "OpenMP enabled: Using -fopenmp flag (manual)")
        elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopenmp")
            message(STATUS "OpenMP enabled: Using -qopenmp flag (manual)")
        else()
            message(WARNING "OpenMP requested but compiler support not detected")
        endif()
    endif()
    # Note: OpenMP_C not found is OK on Mac - we only need Fortran support
    if(OpenMP_Fortran_FOUND AND NOT OpenMP_C_FOUND)
        message(STATUS "OpenMP: Fortran support found (C support not required)")
    endif()
else()
    message(STATUS "OpenMP support disabled (use -DENABLE_OPENMP=ON to enable)")
endif()

# OpenACC GPU acceleration support
option(ENABLE_OPENACC "Enable OpenACC GPU acceleration" OFF)

if(ENABLE_OPENACC)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
        # NVIDIA HPC SDK / PGI Compiler
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc -gpu=cc70,cc75,cc80,cc86,cc89,cc90")
        set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Minfo=accel")
        message(STATUS "OpenACC enabled: Using NVIDIA HPC SDK compiler")
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        # GCC with OpenACC support
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenacc")
        set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fopt-info-optimized-omp")
        message(STATUS "OpenACC enabled: Using GCC with OpenACC support")
    else()
        message(WARNING "OpenACC requested but compiler may not support it. Compiler: ${CMAKE_Fortran_COMPILER_ID}")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenacc")
    endif()
else()
    message(STATUS "OpenACC support disabled (use -DENABLE_OPENACC=ON to enable)")
endif()

if(CMAKE_Fortran_COMPILER MATCHES "gfortran*")
    set(CMAKE_SKIP_BUILD_RPATH            FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH    TRUE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

    #Compiler gfortran Flags
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -fimplicit-none -fopenmp")   #  -lhdf5_fortran -fintrinsic-modules-path ${HDF5_INCLUDE_DIRS}") # /usr/local/include") #  TODO find metter way to include ${HDF5_INCLUDE_DIRS}")  #/usr/local/include")  # -lhdf5hl_fortran")          # -std=f2003 ${CMAKE_Fortran_FLAGS}
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
# Libraries
# -------------------------------------------------------------------------------------------------------------------------

# fotrasn objs
SET(RANDOMNUMBER randomnumber)
ADD_LIBRARY(${RANDOMNUMBER} STATIC ${SOURCE_FOLDER}/randomnumber.f03)
target_link_libraries(${RANDOMNUMBER})

SET(LOGGING logging)
ADD_LIBRARY(${LOGGING} STATIC ${SOURCE_FOLDER}/logging.f03)
target_link_libraries(${LOGGING})

SET(KINDS kinds)
ADD_LIBRARY(${KINDS} STATIC ${SOURCE_FOLDER}/kinds.f03)
target_link_libraries(${KINDS})

SET(TIMERS timers)
ADD_LIBRARY(${TIMERS} STATIC ${SOURCE_FOLDER}/timers.f03)
target_link_libraries(${TIMERS} ${KINDS})

SET(LIB_CONF lib_conf)
ADD_LIBRARY(${LIB_CONF} STATIC ${SOURCE_FOLDER}/lib_conf.f90)
target_link_libraries(${LIB_CONF} ${LOGGING} ${RANDOMNUMBER} ${TIMERS} ${KINDS})

SET(LATTICE_DATA lattice_data)
ADD_LIBRARY(${LATTICE_DATA} STATIC ${SOURCE_FOLDER}/lattice_data.f03)
target_link_libraries(${LATTICE_DATA})

SET(POLYMER_MODEL polymer_model)
ADD_LIBRARY(${POLYMER_MODEL} STATIC ${SOURCE_FOLDER}/polymer_model.f03)
target_link_libraries(${POLYMER_MODEL} ${LOGGING} ${RANDOMNUMBER} ${TIMERS} ${KINDS} ${LIB_CONF})

SET(POLYMER_MODEL_UNIFIED polymer_model_unified)
ADD_LIBRARY(${POLYMER_MODEL_UNIFIED} STATIC ${SOURCE_FOLDER}/polymer_model_unified.f03)
set_target_properties(${POLYMER_MODEL_UNIFIED} PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(${POLYMER_MODEL_UNIFIED} ${TIMERS} ${LATTICE_DATA} ${LOGGING} ${LIB_CONF})

SET(ANALYSE analyse)
ADD_LIBRARY(${ANALYSE} STATIC ${SOURCE_FOLDER}/analyse.f03)
target_include_directories(${ANALYSE} PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(${ANALYSE} PRIVATE ${LATTICE_DATA} ${POLYMER_MODEL} ${RANDOMNUMBER} ${TIMERS} ${KINDS} ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})

SET(ANALYSE_UNIFIED analyse_unified)
ADD_LIBRARY(${ANALYSE_UNIFIED} STATIC ${SOURCE_FOLDER}/analyse_unified.f03)
set_target_properties(${ANALYSE_UNIFIED} PROPERTIES LINKER_LANGUAGE Fortran)
target_include_directories(${ANALYSE_UNIFIED} PRIVATE ${HDF5_INCLUDE_DIRS})
add_dependencies(${ANALYSE_UNIFIED} ${ANALYSE})
target_link_libraries(${ANALYSE_UNIFIED} PRIVATE ${LATTICE_DATA} ${POLYMER_MODEL} ${ANALYSE} ${RANDOMNUMBER} ${TIMERS} ${KINDS} ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})


# include_directories(src)

add_executable(${TARGET} ${SOURCE_FILES})
##target_include_directories(${TARGET} PRIVATE ${HDF5_INCLUDE_DIRS})
# Ensure libraries are built before executable (for module files)
add_dependencies(${TARGET} ${POLYMER_MODEL_UNIFIED} ${POLYMER_MODEL} ${ANALYSE} ${ANALYSE_UNIFIED} ${LATTICE_DATA} ${LIB_CONF} ${TIMERS} ${LOGGING} ${RANDOMNUMBER} ${KINDS})
target_link_libraries(${TARGET} PRIVATE ${RANDOMNUMBER} ${KINDS} ${TIMERS} ${LIB_CONF} ${LATTICE_DATA} ${POLYMER_MODEL} ${POLYMER_MODEL_UNIFIED} ${ANALYSE} ${ANALYSE_UNIFIED} MPI::MPI_Fortran) ## ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})

# Link OpenMP if enabled (must be after executable is created)
if(ENABLE_OPENMP AND OPENMP_FORTRAN_FOUND)
    target_link_libraries(${TARGET} PRIVATE OpenMP::OpenMP_Fortran)
endif()

set_target_properties(${TARGET} PROPERTIES LINKER_LANGUAGE Fortran)

# data files
#configure_file(${EXECUTABLE_OUTPUT_PATH}/${TARGET} ${TARGET} COPYONLY)
#configure_file(${SOURCE_FOLDER}/__init__.py ${EXECUTABLE_OUTPUT_PATH}/__init__.py COPYONLY)
#configure_file(${PACKAGE}/run_3dpolys_le.sh ${EXECUTABLE_OUTPUT_PATH}/run_3dpolys_le.sh COPYONLY)
