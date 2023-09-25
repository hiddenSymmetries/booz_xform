Getting started
===============


Requirements
^^^^^^^^^^^^

``booz_xform`` requires a C++ compiler and python3.  The python
interface uses arrays from the ``numpy`` package, and building
``booz_xform`` requires the `pybind11 package
<https://pybind11.readthedocs.io/en/stable/>`_.  These packages are
installed automatically by the ``pip install`` step described in the
next section.

``booz_xform`` also requires the NetCDF library. The C++ or fortran
interfaces to NetCDF are not required, only the standard C interface.

OpenMP is an optional dependency. If found, OpenMP is used to
parallelize the calculation over magnetic surfaces.  MPI is not used.

There are some additional python packages that are optional dependencies.
The :doc:`plotting` routines require the `matplotlib
<https://matplotlib.org/>`_ package.  The python unit tests require
``matplotlib`` as well as ``scipy``.  These packages are not needed for
carrying out the coordinate transformation, so they are not
automatically installed by ``pip``, and you should install them
separately if you want to use these features.

Installation
^^^^^^^^^^^^

Before attempting to install ``booz_xform``, make sure NetCDF 
is available on your system. On some HPC systems, this may require
loading the relevant module.

There are several ways you can install ``booz_xform``, depending on
whether you are a user or developer.

1. Installation from PyPI
*************************

If you do not plan to edit the source code, the recommended way to
install ``booz_xform`` is to install the latest release from `PyPI
<https://pypi.org/project/booz_xform/>`_ using ``pip``::

    pip install -v booz_xform

This command will download and compile the code. At the start of the
compilation step, the ``cmake`` build system will search for the
NetCDF and MPI header files and libraries.  Any of the environment
variables ``NETCDF_DIR``, ``NETCDF_HOME``, or ``NETCDFDIR`` can be set
to guide ``cmake`` in this search, and it will also look in standard
locations such as ``/opt/local/include``.

The ``-v`` flag above (for verbose output) is optional, but it is
useful since it allows you to see which compiler, NetCDF, and MPI
libraries were used. This information can be found in the lines
similar to the following in the output::

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
  ...

To change the compiler that is used to build the code, you can insert
``CXX=`` followed by the compiler name before ``pip``. For example, to
use the Intel compiler ``icpc``, use

.. code-block::

  CXX=icpc pip install -v booz_xform
  
If the installation is successful, ``booz_xform`` will be added to
your python environment. You should now be able to import the module
from python::

  >>> import booz_xform

On some systems, you may not have permission to install packages to
the default location. In this case, add the ``--user`` flag to ``pip``
so the package can be installed for your user only::

    pip install -v --user booz_xform

  
2. Installation from a local copy of the repository
***************************************************

If you prefer to see or edit the source code, you can first clone the
repository using

.. code-block::

    git clone https://github.com/hiddenSymmetries/booz_xform.git

Then install the package to your local python environment with

.. code-block::

  cd booz_xform
  pip install -v -e .

The ``-e`` flag is not mandatory but it can be helpful during
development. This flag makes the installation editable, in that edits
to the pure python source in ``src/booz_xform`` are immediately
reflected in the package you import into python. Without this flag,
you would need to re-install the package for changes to become active.

The ``pip install`` line can be preceded by ``CXX=`` to select a
specific compiler, as in the PyPI method above.

Again, if you encounter a permissions error trying to install packages
to the default location, add the ``--user`` flag::

    pip install -v -e --user .


3. Installation without pip from a local copy of the repository
***************************************************************

This option is similar to option 2, and is well suited for code
development. This option gives a somewhat faster installation than
option 2, but you must first manually ensure all the dependencies in
``pyproject.toml`` are installed.

From a cloned copy of the repository, run

.. code-block::

  python setup.py develop

This command can be preceded by ``CXX=`` to select a
specific compiler, as in the other methods above.


4. Building outside of the python package system
************************************************

If you are actively developing the code, you may wish to compile the
C++ source without going through the python package system. In this
case, you should have a local copy of the repository, obtained with

.. code-block::

  git clone https://github.com/hiddenSymmetries/booz_xform.git

You must also have the ``pybind11`` python package installed, as well
as ``cmake``.  The code then can be built using the usual approach for
a ``cmake`` project::

  cd booz_xform/build
  cmake ..
  make -j

In this case, the python extension library ``_booz_xform`` (with a
filename usually ending in ``.so``), the standalone executable
``xbooz_xform``, and the library ``libbooz_xform.a`` will all be
created in the ``build`` directory. Note that in this approach, no
python package is installed.  You can import only the ``Booz_xform``
class with ``import _booz_xform``, which loads the compiled extension
without importing the pure python functions.
