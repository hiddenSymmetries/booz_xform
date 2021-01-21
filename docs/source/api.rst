API Reference
=============


C++
^^^

All definitions for the ``booz_xform`` code reside in the namelist ``booz_xform::``,
which will be omitted for brevity here.

All floating point values have the type ``boozfloat``, which is merely a ``typedef`` for ``double``. (The entire code could be switched to single precision by changing this ``typedef`` if desired.)
All 1D arrays of integers (such as ``xm``) have the type ``IntVector``,
which is merely a ``typedef`` for ``std::valarray<int>``.
All 1D arrays of floating point values have the type ``Vector``,
which is merely a ``typedef`` for ``std::valarray<boozfloat>``. 2D arrays of
floating point values are stored with the class ``Matrix``. This
class does not have many features of matrices, such as matrix-vector multiplication,
just the small number of operations needed for the Boozer-coordinate transformation,
such as these:
      
.. doxygenclass:: booz_xform::Matrix
   :members:

The main functionality of the ``booz_xform`` code is contained in a single class,
``Booz_xform``:

.. doxygenclass:: booz_xform::Booz_xform
   :members:

..
   :undoc-members:
      
Python
^^^^^^

