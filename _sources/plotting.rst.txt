Plotting
========

The ``booz_xform`` python module includes several routines for plotting.
These functions take a :class:`~booz_xform.Booz_xform` instance as an argument.
This object can be one which was used to drive the coordinate transformation.
Or, this :class:`~booz_xform.Booz_xform` instance could be one in which results from an earlier
transformation were loaded using the :meth:`~booz_xform.Booz_xform.read_boozmn` function.

The plotting routine require the python `matplotlib
<https://matplotlib.org/>`_ package. This package must be installed
manually since it is not required for the core functionality of
``booz_xform`` and so is not installed by ``pip``.

All plotting routines use matplotlib's current axis. When using these
routines in a script, you typically need to call ``plt.show()`` at the
end to actually display the figure.

A `gallery of plots that can be generated
<https://github.com/hiddenSymmetries/booz_xform/blob/main/examples/plots_demo.ipynb>`_
can be found in ``plots_demo`` notebook in the ``examples`` directory.
The full API of the available plotting routines follow.

.. automodule:: booz_xform
   :members: surfplot, symplot, modeplot


