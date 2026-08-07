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

Python tests are based on the standard ``unittest`` python module, and
can also be run with ``pytest``. Source code for the python tests is located in
the ``tests`` directory. The python tests will use the installed version
of the ``booz_xform`` python package, not necessarily a shared library
compiled manually in the ``build`` directory. To run the python tests,
call

.. code-block::

   pytest

from the repository home directory. You can also run a single file with
``pytest tests/test_regression.py``, or execute individual ``*.py`` test
files directly. ``pip install ".[test]"`` installs ``pytest`` along with
the other packages the tests need.

The python regression tests make use of files in the ``tests/test_files`` directory.


Continuous integration
**********************

The C++ and python tests are automatically run after every commit to
the repository.  This automation is handled by GitHub Actions, and
controlled by the script ``.github/workflows/ci.yml``.
To view the results of the continuous integration runs, you can click on the "Actions"
link from the `GitHub repository page <https://github.com/hiddenSymmetries/booz_xform>`_,
or you can directly visit `<https://github.com/hiddenSymmetries/booz_xform/actions>`_.


Versions and releases
^^^^^^^^^^^^^^^^^^^^^

The version number is declared in exactly one place, ``[project]
version`` in ``pyproject.toml``. Everything else -- the wheel and sdist
metadata that PyPI displays, ``booz_xform.__version__``, the version
compiled into the C++ code and recorded in ``boozmn_*.nc`` output files,
the version shown in this documentation, and the git tag -- is derived
from it. ``tests/test_version.py`` and the ``check-version`` job in
``.github/workflows/release.yml`` both fail if those derivations ever
disagree.

Releases are cut by tagging with `tbump
<https://github.com/your-tools/python-tbump>`_ and publishing a GitHub
release, which triggers ``.github/workflows/release.yml``. That workflow
builds binary wheels for Linux and macOS with `cibuildwheel
<https://cibuildwheel.pypa.io/>`_, smoke-tests each one, and uploads
them together with the source distribution to PyPI. The full procedure
is documented for maintainers in ``dev/RELEASING.md`` in the repository.
