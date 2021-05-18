API Reference
=============

Python
^^^^^^

The ``booz_xform`` module provides a class ``Booz_xform``, described
in detail here.  The same python module also includes several
functions for plotting, described on the :doc:`plotting` page.

The class ``Booz_xform`` provides a small number of functions to load
input data, drive the calculation, and save output data.  In addition,
all the input and output data are available as properties of the
class. The properties are scalars, 1D numpy arrays, and 2D numpy
arrays.

Properties of the class are labeled below as being an input or
output. The input quantities can be either read or written to in
python, and should be set before calling
:func:`~booz_xform.Booz_xform.run()`. (Many of these properties can be
set by calling :func:`~booz_xform.Booz_xform.read_wout()`).  The
output quantities are read-only, and they are populated when
:func:`~booz_xform.Booz_xform.run()` is called.

.. autoclass:: booz_xform.Booz_xform
   :members:
   :undoc-members:
   :member-order: groupwise

.. note:: While 1D and 2D array parameters can appear on the left-hand
          side of assignment operations, individual array elements or
          array slices cannot. So for instance, if ``b`` is an
          instance of :class:`booz_xform.Booz_Xform`, ``b.rmnc = [[1,
          0.1],[1.1, 0.1]]`` will behave as expected, but
          ``b.rmnc[0,0] = 2.0`` will not. This is due to the copying
          of data when interfacing python to C++.


C++
^^^

All definitions for the ``booz_xform`` code reside in the namespace ``booz_xform::``,
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

Public variables below are labeled below as being an input or
output. The input quantities should be set before calling run(). (Many
of these input quantities can be set by calling read_wout()).  The
output quantities are populated when run() is called.

The main functionality of the ``booz_xform`` code is contained in a single class,
``Booz_xform``:

.. doxygenclass:: booz_xform::Booz_xform
   :members:

..
   :undoc-members:
      
