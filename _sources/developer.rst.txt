Developer notes
===============

Testing
^^^^^^^

The ``booz_xform`` package includes both C++ and python tests,
including unit and regression tests, and continuous integration.

C++ tests
*********

Unit tests of the C++ code are handled using the `doctest
<https://github.com/onqtam/doctest>`_ header-only library. A copy of
the library in ``externalPackages/doctest`` is used. Source code for
the C++ tests is located in the ``tests`` directory. Whenever
``booz_xform`` is compiled, a ``unitTests`` executable is compiled
containing the C++ tests. To access this executable, you should
compile the code manually (i.e. outside of ``pip``) using ``cd build;
cmake ..; make -j``, which will create the ``unitTests`` executable in
the ``build`` directory.

Python tests
************

The python tests require the packages ``scipy`` and ``matplotlib``,
which the core part of ``booz_xform`` does not require.

Python tests are based on the standard ``unittest`` python module.
Source code for the python tests is located in the ``tests`` directory.
The python tests will use the installed version of the ``booz_xform`` python package,
not necessarily a shared library compiled manually in the ``build`` directory.
One way to run the python tests is to call

.. code-block::

   python -m unittest

from the repository home directory or from the ``tests``
directory. You can also execute individual ``*.py`` test files
directly.

The python regression tests make use of files in the ``tests/test_files`` directory.


Continuous integration
**********************

The C++ and python tests are automatically run after every commit to
the repository.  This automation is handled by GitHub Actions, and
controlled by the script ``.github/workflows/ci.yml``.
To view the results of the continuous integration runs, you can click on the "Actions"
link from the `GitHub repository page <https://github.com/hiddenSymmetries/booz_xform>`_,
or you can directly visit `<https://github.com/hiddenSymmetries/booz_xform/actions>`_.
