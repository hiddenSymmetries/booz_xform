# Versions >= 3.9 have improved support for MPI: https://cliutils.gitlab.io/modern-cmake/chapters/packages/MPI.html
# The command add_compile_definitions, used in this file, was new in CMake version 3.12.
cmake_minimum_required(VERSION 3.12)

# Comment out the next line when speed is needed
set(CMAKE_BUILD_TYPE Debug)
add_definitions(-O2)

#set(CMAKE_MODULE_PATH cmake)
include(cmake/knownHosts.cmake)

project(booz_xform LANGUAGES CXX)

message("CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")
message("CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package (NetCDF REQUIRED)
message("NETCDF_INCLUDES=${NETCDF_INCLUDES}")
message("NETCDF_LIBRARIES=${NETCDF_LIBRARIES}")

find_package(MPI REQUIRED)

# Tell "make" to print out the commands used for compiling and linking:
set(CMAKE_VERBOSE_MAKEFILE on)

set(SOURCE_DIR "src")
include_directories(${SOURCE_DIR})
include_directories(externalPackages/doctest)
include_directories(externalPackages/eigen)
include_directories("${NETCDF_INCLUDES}")

# Get a list of all sources:
file(GLOB SOURCES "${SOURCE_DIR}/*.cpp")
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
pybind11_add_module(booz_xform ${SOURCES} ${SOURCE_DIR}/bindings.cpp)

add_library(libbooz_xform ${SOURCES})
# The next line makes the library name libbooz_xform.a instead of liblibbooz_xform.a:
set_property(TARGET libbooz_xform PROPERTY OUTPUT_NAME booz_xform)

add_executable(xbooz_xform ${SOURCE_DIR}/xbooz_xform.cpp)
#set_property(TARGET xbooz_xform PROPERTY CXX_STANDARD 11)

target_link_libraries(libbooz_xform PUBLIC ${NETCDF_LIBRARIES})

# Driver and python module depend on the library:
target_link_libraries(xbooz_xform PUBLIC libbooz_xform)
target_link_libraries( booz_xform PUBLIC libbooz_xform)

# Set up unitTests executable
set(BOOZ_TESTS unitTests)
file(GLOB UNIT_TEST_SOURCES "tests/*.cpp")
message("Tests: ${UNIT_TEST_SOURCES}")
add_executable(${BOOZ_TESTS} ${UNIT_TEST_SOURCES})
target_link_libraries(${BOOZ_TESTS} PUBLIC libbooz_xform)
# Put the unitTests executable in the "tests" directory:
#set_property(TARGET ${BOOZ_TESTS} PROPERTY RUNTIME_OUTPUT_DIRECTORY tests)

# We need the c++11 standard both for doctest, and because the Intel compiler
# does not recognize std::begin/end otherwise.
set_property(TARGET libbooz_xform PROPERTY CXX_STANDARD 11)
set_property(TARGET ${BOOZ_TESTS} PROPERTY CXX_STANDARD 11)

# Do not build the unit tests by default, since they take some time to build:
#set_property(TARGET ${BOOZ_TESTS} PROPERTY EXCLUDE_FROM_ALL 1)


# Set up "make test" so it runs the scripts I want it to run.
# See https://cmake.org/cmake/help/latest/command/add_custom_target.html
#add_custom_target(test COMMAND ./runTests WORKING_DIRECTORY tests DEPENDS ${BOOZ_TESTS})
#add_custom_target(test COMMAND ./ WORKING_DIRECTORY tests DEPENDS ${BOOZ_TESTS})
