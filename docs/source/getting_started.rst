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

Then, the recommended way to install ``booz_xform`` is to install it from `PyPI <https://pypi.org/project/booz_xform/>`_ using ``pip``::

    pip install booz_xform

This command will download and compile the code. At the start of the compilation step,
the ``cmake`` build system will search for the NetCDF and MPI header files and libraries.
Any of the environment variables ``NETCDF_DIR``, ``NETCDF_HOME``, or ``NETCDFDIR``
can be set to guide ``cmake`` in this search, and it will also look in standard locations such as ``/opt/local/include``.

If you prefer to see or edit the source code, you can first clone the repository using::

    git clone https://github.com/landreman/booz_xform.git

Then install the package to your local python environment with::

  cd booz_xform
  pip install -e .

The ``-e`` in the last command is optional.

