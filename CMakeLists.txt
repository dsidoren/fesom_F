cmake_minimum_required(VERSION 3.16)
project(fesom3 Fortran)

# ----------------------------------------------------------------------
# FESOM3 — clean-architecture Fortran reimplementation of FESOM2.
# v1 anchor build (bit-identity oracle = FESOM2 v2.7.3): Intel + double
# precision. See cmake/fesom_flags.cmake and docs/LESSONS.md.
# ----------------------------------------------------------------------

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Default to the anchor (Release) build unless told otherwise.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Build type (Release = anchor; Debug = +checks, FP-identical)" FORCE)
endif()

# Precision options (mutually exclusive; DP is default = anchor).
option(USE_SINGLE_PRECISION "Build with single precision (WP=4)" OFF)
option(USE_HALF_PRECISION   "Build with half precision (WP=2, FP16; NVHPC only)" OFF)
if(USE_HALF_PRECISION AND USE_SINGLE_PRECISION)
    message(FATAL_ERROR "USE_HALF_PRECISION and USE_SINGLE_PRECISION are mutually exclusive.")
endif()
if(USE_HALF_PRECISION)
    message(STATUS "FESOM3 precision: HALF (WP=2)")
elseif(USE_SINGLE_PRECISION)
    message(STATUS "FESOM3 precision: SINGLE (WP=4)")
else()
    message(STATUS "FESOM3 precision: DOUBLE (WP=8, anchor)")
endif()

include(fesom_flags)

find_package(MPI REQUIRED COMPONENTS Fortran)

# Put .mod files in a private dir so the scanner resolves USE across subdirs.
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ----------------------------------------------------------------------
# Library: everything under src/ except the driver entry points.
# GLOB + CONFIGURE_DEPENDS so newly added module files are picked up on
# rebuild without re-running cmake by hand. CMake's Fortran dependency
# scanner orders the modules (params -> types -> infra -> ...) automatically.
# ----------------------------------------------------------------------
set(FESOM3_LIB_DIRS params types infra mesh oce ice forcing io step)
set(fesom3_lib_sources "")
foreach(_d ${FESOM3_LIB_DIRS})
    file(GLOB _s CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/${_d}/*.F90)
    list(APPEND fesom3_lib_sources ${_s})
endforeach()

add_library(fesom3 STATIC ${fesom3_lib_sources})
apply_fesom_compile_flags(fesom3)
target_link_libraries(fesom3 PUBLIC MPI::MPI_Fortran)
set_target_properties(fesom3 PROPERTIES LINKER_LANGUAGE Fortran)

# ----------------------------------------------------------------------
# Driver executables: one per file in src/drivers/.
# ----------------------------------------------------------------------
file(GLOB fesom3_drivers CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/drivers/*.F90)
foreach(_drv ${fesom3_drivers})
    get_filename_component(_name ${_drv} NAME_WE)
    add_executable(${_name} ${_drv})
    apply_fesom_compile_flags(${_name})
    target_link_libraries(${_name} PUBLIC fesom3)
    set_target_properties(${_name} PROPERTIES LINKER_LANGUAGE Fortran)
endforeach()

# ----------------------------------------------------------------------
# Tests: one executable per test/test_*.F90, registered with ctest.
# add_fesom_test(<name> <nranks>) builds + runs under mpirun.
# ----------------------------------------------------------------------
enable_testing()

# Oversubscribe-safe launcher (login-node runs of small rank counts).
set(FESOM3_MPIEXEC ${MPIEXEC_EXECUTABLE})
if(NOT FESOM3_MPIEXEC)
    set(FESOM3_MPIEXEC mpirun)
endif()
# The DKRZ UCX runtime settings (env/levante shells) target compute-node
# InfiniBand and fail on login nodes. Self-tests are tiny single-node runs, so
# default to the ob1/vader (shared-memory) path. Override for batch if needed.
set(FESOM3_MPI_LAUNCH_FLAGS "--mca;pml;ob1;--mca;btl;self,vader;--oversubscribe"
    CACHE STRING "Extra mpirun flags for ctest self-tests")

function(add_fesom_test test_name nranks)
    add_executable(${test_name} ${CMAKE_CURRENT_SOURCE_DIR}/test/${test_name}.F90)
    apply_fesom_compile_flags(${test_name})
    target_link_libraries(${test_name} PUBLIC fesom3)
    set_target_properties(${test_name} PROPERTIES
        LINKER_LANGUAGE Fortran
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    add_test(NAME ${test_name}_np${nranks}
             COMMAND ${FESOM3_MPIEXEC} ${FESOM3_MPI_LAUNCH_FLAGS} -np ${nranks}
                     $<TARGET_FILE:${test_name}>)
endfunction()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test/CMakeLists.txt)
    add_subdirectory(test)
endif()
