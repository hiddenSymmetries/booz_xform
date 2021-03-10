=========================
booz\_xform documentation
=========================

``booz_xform`` is a package for computing Boozer coordinates in toroidal magnetohydrodynamic
equilibria, including both stellarators and tokamaks.
The package described here follows the same algorithm as the fortran 77 code of the same name
in `Stellopt <https://github.com/PrincetonUniversity/STELLOPT/tree/develop/BOOZ_XFORM>`_.
However the package here is written in C++, with python bindings.
The package here is also written so as to allow input from equilibrium codes other than VMEC,
it is parallelized using OpenMP, and it includes functions for plotting output.
It is also equipped with unit and regression tests and continuous integration.

.. toctree::
   :maxdepth: 3
   :caption: Contents

   getting_started
   theory
   usage
   plotting
   api
   developer
   source
