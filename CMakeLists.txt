cmake_minimum_required(VERSION 3.4)
project(example LANGUAGES CXX)

#include_directories(src)

# Get a list of all sources:
file(GLOB SOURCES "src/*.cpp")
message("Sources: ${SOURCES}")
# In the future, a driver could be removed this way:
#list(FILTER SOURCES EXCLUDE REGEX ".*xqsc.cpp$")
#message("Sources after filter: ${SOURCES}")

# Locate pybind11
execute_process(COMMAND python -m pybind11 --cmakedir
                OUTPUT_VARIABLE pybind11_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE)
message("Result of python -m pybind11 --cmakedir: ${pybind11_DIR}")
#set(pybind11_DIR "/Users/mattland/opt/miniconda3/envs/20201118-01/lib/python3.8/site-packages/pybind11/share/cmake/pybind11")
#message("New pybind11_DIR: ${pybind11_DIR}")

find_package(pybind11 REQUIRED PATHS pybind11_DIR)
pybind11_add_module(booz_xform ${SOURCES})
