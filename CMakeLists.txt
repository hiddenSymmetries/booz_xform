# Versions >= 3.9 have improved support for MPI: https://cliutils.gitlab.io/modern-cmake/chapters/packages/MPI.html
# 3.18 is required because pybind11 >= 3's CMake package config needs >= 3.15,
# and because CMake 4 warns about compatibility requests below 3.10.
cmake_minimum_required(VERSION 3.18)

set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

# Default to an optimized build, but respect -DCMAKE_BUILD_TYPE=... on the
# command line.  setup.py passes Release, so pip-installed packages and the
# wheels on PyPI are optimized.  Use -DCMAKE_BUILD_TYPE=Debug for debugging;
# that also re-enables Eigen's bounds assertions.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# This next compile flag is needed to avoid a compile error on some systems. See
# https://pybind11.readthedocs.io/en/stable/faq.html#recursive-template-instantiation-exceeded-maximum-depth-of-256
add_definitions(-ftemplate-depth=1024)

#set(CMAKE_MODULE_PATH cmake)
include(cmake/knownHosts.cmake)

project(booz_xform LANGUAGES CXX)

message("CMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")
message("CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package (NetCDF REQUIRED)
message("NETCDF_INCLUDES=${NETCDF_INCLUDES}")
message("NETCDF_LIBRARIES=${NETCDF_LIBRARIES}")

# If in the future we want to use MPI, uncomment this next line:
#find_package(MPI REQUIRED)

# Apple's clang can compile OpenMP, but Apple does not ship the runtime
# library. Homebrew's libomp supplies it, and it is the only practical way to
# get an OpenMP build on macOS without switching compilers. The catch is that
# the formula is keg-only: it is deliberately NOT symlinked into the Homebrew
# prefix, so FindOpenMP finds neither libomp.dylib nor omp.h in the usual
# places, gives up, and silently configures a serial build. Point it at the keg
# instead. OpenMP_ROOT is honored by every find_ call inside FindOpenMP (policy
# CMP0074), which is what makes both the library and the header discoverable.
# An OpenMP_ROOT set by the user, or in the environment, always wins.
if(APPLE AND NOT DEFINED OpenMP_ROOT AND NOT DEFINED ENV{OpenMP_ROOT})
  find_program(BREW_EXECUTABLE brew)
  if(BREW_EXECUTABLE)
    execute_process(COMMAND ${BREW_EXECUTABLE} --prefix libomp
                    RESULT_VARIABLE brew_libomp_result
                    OUTPUT_VARIABLE brew_libomp_prefix
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET)
    # `brew --prefix <formula>` succeeds and prints a path even when the
    # formula is not installed, so check that the header is really there.
    if(brew_libomp_result EQUAL 0 AND EXISTS "${brew_libomp_prefix}/include/omp.h")
      set(OpenMP_ROOT "${brew_libomp_prefix}")
      message(STATUS "Using Homebrew's keg-only libomp for OpenMP: ${OpenMP_ROOT}")
    endif()
  endif()
endif()

find_package(OpenMP)

# OpenMP is optional, and a build without it is correct, just serial. That is
# easy to miss in a build log, and a wheel that quietly lost its parallelism is
# worth failing the build over, so the wheel builds (and CI) set this option.
option(BOOZ_XFORM_REQUIRE_OPENMP
       "Fail the build if OpenMP is not found, instead of building serially" OFF)
if(BOOZ_XFORM_REQUIRE_OPENMP AND NOT OpenMP_CXX_FOUND)
  message(FATAL_ERROR
    "BOOZ_XFORM_REQUIRE_OPENMP is ON but OpenMP was not found.\n"
    "On macOS, install Homebrew's libomp ('brew install libomp'), or set\n"
    "OpenMP_ROOT to a directory containing include/omp.h and lib/libomp.dylib.")
endif()

# Tell "make" to print out the commands used for compiling and linking:
set(CMAKE_VERBOSE_MAKEFILE on)

set(SOURCE_DIR "src/_booz_xform")
include_directories(${SOURCE_DIR})
include_directories(externalPackages/doctest)
include_directories(externalPackages/eigen)
include_directories("${NETCDF_INCLUDES}")
if(OpenMP_CXX_FOUND)
  message("OpenMP_CXX_INCLUDE_DIRS is ${OpenMP_CXX_INCLUDE_DIRS}")
  include_directories(${OpenMP_CXX_INCLUDE_DIRS})
  add_definitions(-DOPENMP)
endif()

# Get a list of all sources:
file(GLOB SOURCES "${SOURCE_DIR}/*.cpp")
message("Sources: ${SOURCES}")
# Remove the 2 special files: the standalone executable and the python bindings:
list(FILTER SOURCES EXCLUDE REGEX ".*xbooz_xform.cpp$")
list(FILTER SOURCES EXCLUDE REGEX ".*bindings.cpp$")
message("Sources after filter: ${SOURCES}")

# Locate pybind11.  setup.py passes -DPython_EXECUTABLE, so an install always
# asks the interpreter being installed into; a bare "cmake .." falls back to
# whatever python is on PATH.
if(NOT DEFINED Python_EXECUTABLE)
  set(Python_EXECUTABLE python)
endif()
execute_process(COMMAND ${Python_EXECUTABLE} -m pybind11 --cmakedir
                OUTPUT_VARIABLE pybind11_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE)
message("Result of ${Python_EXECUTABLE} -m pybind11 --cmakedir: ${pybind11_DIR}")

find_package(pybind11 REQUIRED PATHS ${pybind11_DIR})
pybind11_add_module(_booz_xform ${SOURCE_DIR}/bindings.cpp)

add_library(libbooz_xform ${SOURCES})
# The next line makes the library name libbooz_xform.a instead of liblibbooz_xform.a:
set_property(TARGET libbooz_xform PROPERTY OUTPUT_NAME booz_xform)
set_property(TARGET libbooz_xform PROPERTY POSITION_INDEPENDENT_CODE ON)

add_executable(xbooz_xform ${SOURCE_DIR}/xbooz_xform.cpp)
#set_property(TARGET xbooz_xform PROPERTY CXX_STANDARD 11)

target_link_libraries(libbooz_xform PUBLIC ${NETCDF_LIBRARIES})
if(OpenMP_CXX_FOUND)
  target_link_libraries(libbooz_xform PUBLIC OpenMP::OpenMP_CXX)
endif()
  
# Driver and python module depend on the library:
target_link_libraries(xbooz_xform PUBLIC libbooz_xform)
target_link_libraries(_booz_xform PUBLIC libbooz_xform)

# Pass the version number from pyproject.toml (via setup.py) to C++, following
# https://github.com/pybind/cmake_example/blob/master/CMakeLists.txt
# booz_xform::version is a namespace-scope const, so it has internal linkage and
# every target that uses it needs the definition: bindings.cpp exposes it as
# _booz_xform.__version__, driver.cpp prints it, and write_boozmn.cpp records it
# in the output file.  Left undefined for builds that do not go through
# setup.py, in which case the C++ code reports "development version".
if(BOOZ_XFORM_VERSION)
  foreach(target _booz_xform libbooz_xform xbooz_xform)
    target_compile_definitions(${target} PRIVATE VERSION_INFO=${BOOZ_XFORM_VERSION})
  endforeach()
endif()

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
