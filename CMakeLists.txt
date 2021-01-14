# Versions >= 3.9 have improved support for MPI: https://cliutils.gitlab.io/modern-cmake/chapters/packages/MPI.html
# The command add_compile_definitions, used in this file, was new in CMake version 3.12.
cmake_minimum_required(VERSION 3.12)

# Comment out the next line when speed is needed
set(CMAKE_BUILD_TYPE Debug)

# If on a Mac, use the Accelerate framework for BLAS/LAPACK.
# Otherwise, CMake may use a random BLAS/LAPACK it finds that is not as optimized.
# To check if we're on a Mac, evaluate 'uname -s'.
execute_process(COMMAND uname -s
                OUTPUT_VARIABLE UNAME_S
                OUTPUT_STRIP_TRAILING_WHITESPACE)
message("Result of uname -s: ${UNAME_S}")
if ("${UNAME_S}" STREQUAL "Darwin")
  message("Setting BLA_VENDOR to Apple")
  set(BLA_VENDOR Apple)
endif()

#set(CMAKE_MODULE_PATH cmake)
include(cmake/knownHosts.cmake)

project(booz_xform LANGUAGES CXX)

message("CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")
message("CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package (NetCDF REQUIRED)
message("NETCDF_INCLUDES=${NETCDF_INCLUDES}")
message("NETCDF_LIBRARIES=${NETCDF_LIBRARIES}")

message("BLA_VENDOR: ${BLA_VENDOR}")
find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)

# Tell "make" to print out the commands used for compiling and linking:
set(CMAKE_VERBOSE_MAKEFILE on)

include_directories(src)
include_directories("${NETCDF_INCLUDES}")

# Get a list of all sources:
file(GLOB SOURCES "src/*.cpp")
message("Sources: ${SOURCES}")
# Remove the 2 special files: the standalone executable and the python bindings:
list(FILTER SOURCES EXCLUDE REGEX ".*xbooz_xform.cpp$")
list(FILTER SOURCES EXCLUDE REGEX ".*bindings.cpp$")
message("Sources after filter: ${SOURCES}")

# Locate pybind11
execute_process(COMMAND python -m pybind11 --cmakedir
                OUTPUT_VARIABLE pybind11_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE)
message("Result of python -m pybind11 --cmakedir: ${pybind11_DIR}")
#set(pybind11_DIR "/Users/mattland/opt/miniconda3/envs/20201118-01/lib/python3.8/site-packages/pybind11/share/cmake/pybind11")
#message("New pybind11_DIR: ${pybind11_DIR}")

find_package(pybind11 REQUIRED PATHS pybind11_DIR)
pybind11_add_module(booz_xform ${SOURCES} src/bindings.cpp)

add_library(libbooz_xform ${SOURCES})
# The next line makes the library name libbooz_xform.a instead of liblibbooz_xform.a:
set_property(TARGET libbooz_xform PROPERTY OUTPUT_NAME booz_xform)

add_executable(xbooz_xform src/xbooz_xform.cpp)
#set_property(TARGET xbooz_xform PROPERTY CXX_STANDARD 11)

target_link_libraries(libbooz_xform PUBLIC ${BLAS_LIBRARIES} ${NETCDF_LIBRARIES})

# Driver and python module depend on the library:
target_link_libraries(xbooz_xform PUBLIC libbooz_xform)
target_link_libraries( booz_xform PUBLIC libbooz_xform)
