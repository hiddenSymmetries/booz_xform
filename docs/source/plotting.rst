Plotting
========

The ``booz_xform`` python module includes several routines for plotting.
These functions take a :class:`~booz_xform.Booz_xform` instance as an argument.
This object can be one which was used to drive the coordinate transformation.
Or, this :class:`~booz_xform.Booz_xform` instance could be one in which results from an earlier
transformation were loaded using the :meth:`~booz_xform.Booz_xform.read_boozmn` function.

The available plotting routines follow.

.. automodule:: booz_xform
   :members: surfplot, symplot


