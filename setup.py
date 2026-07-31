# -*- coding: utf-8 -*-

# All package metadata, including the version, lives in pyproject.toml.  This
# file exists only to drive the CMake build of the _booz_xform extension.
#
# It was adapted from the "official" pybind11 example at
# https://github.com/pybind/cmake_example
#
# See also https://www.benjack.io/2018/02/02/python-cpp-revisited.html

import os
import shlex
import shutil
import sys
import subprocess
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

# Convert distutils Windows platform specifiers to CMake -A arguments
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}


# CMakeCache.txt entries that hold absolute paths which can go stale between
# builds.  See cmake_cache_is_usable() below.
CACHED_PATH_ENTRIES = (
    "CMAKE_MAKE_PROGRAM",
    "CMAKE_CXX_COMPILER",
    "CMAKE_C_COMPILER",
    "pybind11_DIR",
    "Python_EXECUTABLE",
    "PYTHON_EXECUTABLE",
)


def cmake_cache_is_usable(build_temp, sourcedir):
    """Can the CMake cache left in build_temp by an earlier build be reused?

    The build directory persists between runs, but the tools it points at do
    not.  Every `pip install` builds in a fresh, randomly-named isolated
    environment, so the ninja executable and pybind11 CMake directory that CMake
    recorded during one run live in a directory pip deletes when that run ends.
    CMake does not notice on the next run: it reads CMAKE_MAKE_PROGRAM from the
    cache and dies with "no such file or directory", and nothing short of
    deleting the build directory by hand fixes it.  So the second and every
    later `pip install .` in a working tree used to fail.

    Rather than force a from-scratch configure every time, which would throw
    away incremental rebuilds for developers, check whether the recorded paths
    still exist, and whether the cache belongs to this source tree at all.
    """
    cache = os.path.join(build_temp, "CMakeCache.txt")
    if not os.path.isfile(cache):
        # Nothing cached, so nothing can be stale.
        return True

    with open(cache, errors="replace") as f:
        for line in f:
            # Cache entries look like NAME:TYPE=VALUE; comments start with //.
            name, separator, value = line.partition("=")
            if not separator or ":" not in name:
                continue
            key = name.split(":", 1)[0].strip()
            value = value.strip()
            if key == "CMAKE_HOME_DIRECTORY":
                # The tree was moved or renamed since the last build.
                if os.path.normpath(value) != os.path.normpath(sourcedir):
                    return False
            elif key in CACHED_PATH_ENTRIES:
                if os.path.isabs(value) and not os.path.exists(value):
                    return False
    return True


# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cfg = "Debug" if self.debug else "Release"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        # The version number comes from [project] version in pyproject.toml and
        # is handed to the C++ compiler as VERSION_INFO; see
        # src/_booz_xform/booz_xform.hpp.  Python code does not use it --
        # booz_xform.__version__ reads the installed distribution metadata --
        # but the C++ driver prints it and write_boozmn() records it in the
        # output file.
        #
        # Python_EXECUTABLE is what pybind11 >= 3 (FindPython) uses;
        # PYTHON_EXECUTABLE is the pybind11 2.x spelling.  Pass both so the
        # extension is always built against the interpreter that is installing
        # it, rather than whichever "python" happens to be first on PATH.
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={}".format(extdir),
            "-DPython_EXECUTABLE={}".format(sys.executable),
            "-DPYTHON_EXECUTABLE={}".format(sys.executable),
            "-DBOOZ_XFORM_VERSION={}".format(self.distribution.get_version()),
            "-DCMAKE_BUILD_TYPE={}".format(cfg),  # not used on MSVC, but no harm
        ]

        # Anything in the CMAKE_ARGS environment variable is passed straight
        # through, the same convention scikit-build uses.  The wheel builds use
        # it for -DBOOZ_XFORM_REQUIRE_OPENMP=ON; it is also the way to pass
        # -DNETCDF_DIR=... or any other CMake option to a `pip install`.
        cmake_args += shlex.split(os.environ.get("CMAKE_ARGS", ""))
        # Only the python extension module is needed for an install.  The
        # standalone xbooz_xform executable and the C++ unitTests executable are
        # separate CMake targets, built by hand or by CI; skipping them here
        # keeps wheel builds from compiling doctest five times per platform.
        build_args = ["--target", "_booz_xform"]

        if self.compiler.compiler_type != "msvc":
            # Using Ninja-build since it a) is available as a wheel and b)
            # multithreads automatically. MSVC would require all variables be
            # exported for Ninja to pick it up, which is a little tricky to do.
            # Users can override the generator with CMAKE_GENERATOR in CMake
            # 3.15+.
            if not cmake_generator:
                cmake_args += ["-GNinja"]

        else:

            # Single config generators are handled "normally"
            single_config = any(x in cmake_generator for x in {"NMake", "Ninja"})

            # CMake allows an arch-in-generator style for backward compatibility
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})

            # Specify the arch if using MSVC generator, but only if it doesn't
            # contain a backward-compatibility arch spec already in the
            # generator name.
            if not single_config and not contains_arch:
                cmake_args += ["-A", PLAT_TO_CMAKE[self.plat_name]]

            # Multi-config generators have a different way to specify configs
            if not single_config:
                cmake_args += [
                    "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)
                ]
                build_args += ["--config", cfg]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call, not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += ["-j{}".format(self.parallel)]

        if not cmake_cache_is_usable(self.build_temp, ext.sourcedir):
            print("Discarding the stale CMake cache in {} and configuring from "
                  "scratch.".format(self.build_temp))
            shutil.rmtree(self.build_temp, ignore_errors=True)

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )


setup(
    ext_modules=[CMakeExtension("booz_xform._booz_xform")],
    cmdclass={"build_ext": CMakeBuild},
)

# For guidance about packages involving both pure python and a pybind11 extension, see
# https://stackoverflow.com/questions/47599162/pybind11-how-to-package-c-and-python-code-into-a-single-package
