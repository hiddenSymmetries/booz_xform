API Reference
=============


C++
^^^

All definitions for the ``booz_xform`` code reside in the namelist ``booz_xform::``,
which will be omitted for brevity here.

1D and 2D arrays are handled using the `Eigen <http://eigen.tuxfamily.org/>`_ library.

All floating point values have the type ``boozfloat``, which is merely
a ``typedef`` for ``double``.  All 1D arrays of integers (such as
``xm``) have the type ``IntVector``, which is merely a ``typedef`` for
``Eigen::ArrayXi``.  All 1D arrays of floating point values have the
type ``Vector``, which is merely a ``typedef`` for
``Eigen::ArrayXd``. 2D arrays of floating point values have the type
``Matrix``, which is a ``typedef`` for ``Eigen::ArrayXXd``.  If
desired, the entire code could be switched to single precision by
changing these ``typedef`` statements, which are defined in
``vector_matrix.hpp``.

The main functionality of the ``booz_xform`` code is contained in a single class,
``Booz_xform``:

.. doxygenclass:: booz_xform::Booz_xform
   :members:

..
   :undoc-members:
      
Python
^^^^^^

