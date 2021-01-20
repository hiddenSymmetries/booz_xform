Getting started
===============

Overview
^^^^^^^^

``booz_xform`` is a package for computing Boozer coordinates in toroidal magnetohydrodynamic
equilibria, including both stellarators and tokamaks.
The package described here follows the same algorithm as the fortran 77 code of the same name
in `Stellopt <https://github.com/PrincetonUniversity/STELLOPT/tree/develop/BOOZ_XFORM>`_.
However the package here is written in C++, with python bindings.
The package here is also written so as to allow input from equilibrium codes other than VMEC.
It is also equipped with unit and regression tests and continuous integration.


Requirements
^^^^^^^^^^^^

``booz_xform`` requires a C++ compiler and python3. It uses `pybind11 package <https://pybind11.readthedocs.io/en/stable/>`_, which is installed automatically by the ``pip install`` step described in the next section.

``booz_xform`` also requires the NetCDF library. The C++ or fortran interfaces to NetCDF are
not required, only the standard C interface.

Finally, ``booz_xform`` requires MPI.


Installation
^^^^^^^^^^^^

Before attempting to install ``booz_xform``, make sure NetCDF and MPI are available on your system. On some HPC systems, this may require loading the relevant modules.

There are several ways you can install ``booz_xform``, depending on whether you are a user or developer.

1. Installation from PyPI
*************************

If you do not plan to edit the source code,
the recommended way to install ``booz_xform`` is to install
the latest release from `PyPI <https://pypi.org/project/booz_xform/>`_ using ``pip``::

    pip install -v booz_xform

This command will download and compile the code. At the start of the compilation step,
the ``cmake`` build system will search for the NetCDF and MPI header files and libraries.
Any of the environment variables ``NETCDF_DIR``, ``NETCDF_HOME``, or ``NETCDFDIR``
can be set to guide ``cmake`` in this search, and it will also look in standard locations such as ``/opt/local/include``.

The ``-v`` flag above (for verbose output) is optional, but it is useful since it allows you to see which compiler, NetCDF, and MPI libraries were used. This information can be found in the lines similar to the following in the output::

  ...
  -- The CXX compiler identification is AppleClang 11.0.0.11000033
  -- Detecting CXX compiler ABI info
  -- Detecting CXX compiler ABI info - done
  -- Check for working CXX compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ - skipped
  -- Detecting CXX compile features
  -- Detecting CXX compile features - done
  CMAKE_MODULE_PATH=
  CMAKE_CURRENT_SOURCE_DIR=/private/var/folders/_2/t14gsms50v1dmfz95bz32hk00000gn/T/pip-req-build-0fgxd2du
  Hello world from FindNetCDF
  -- Found NetCDF: /opt/local/lib/libnetcdf.dylib
  NETCDF_INCLUDES=/opt/local/include
  NETCDF_LIBRARIES=/opt/local/lib/libnetcdf.dylib
  -- Found MPI_CXX: /opt/local/lib/mpich-gcc8/libmpicxx.dylib (found version "3.1")
  ...

To change the compiler that is used to build the code, you can insert ``CXX=`` followed by the compiler name before ``pip``. For example, to use the Intel compiler ``mpiicpc``, use

.. code-block::

  CXX=mpiicpc pip install -v booz_xform
  
If the installation is successful, ``booz_xform`` will be added to your python environment. You should now be able to import the module from python::

  >>> import booz_xform

2. Installation from a local copy of the repository
***************************************************

If you prefer to see or edit the source code, you can first clone the repository using

.. code-block::

    git clone https://github.com/hiddenSymmetries/booz_xform.git

Then install the package to your local python environment with

.. code-block::

  cd booz_xform
  pip install -v .

The last line can be preceded by ``CXX=`` to select a specific compiler, as in the PyPI method above.

3. Building outside of the python package system
************************************************

If you are actively developing the code, you may wish to compile the source without
going through the python package system. In this case, you should have a local copy
of the repository, obtained with

.. code-block::

  git clone https://github.com/hiddenSymmetries/booz_xform.git

The code then can be built using the usual approach for a ``cmake`` project::

  cd booz_xform/build
  cmake ..
  make -j

In this case, the python extension library (usually ending in ``.so``), the standalone executable ``xbooz_xform``,
and the library ``libbooz_xform.a`` will all be created in the ``build`` directory.
