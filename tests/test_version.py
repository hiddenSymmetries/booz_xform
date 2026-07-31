#!/usr/bin/env python3

"""Guards on the version wiring described in dev/RELEASING.md.

The version literal lives in exactly one place, [project] version in
pyproject.toml.  Everything else derives from it, and these tests fail if any
of those derivations drift apart.
"""

import unittest
from importlib.metadata import version

import booz_xform as bx


class VersionTest(unittest.TestCase):
    def test_version_matches_installed_distribution_metadata(self):
        """bx.__version__ must agree with what pip and PyPI see."""
        self.assertEqual(bx.__version__, version("booz_xform"))

    def test_version_is_pep440(self):
        """A PEP 440 parse is the same check PyPI applies at upload time."""
        try:
            from packaging.version import Version
        except ImportError:
            self.skipTest("packaging is not installed")
        Version(bx.__version__)

    def test_cpp_version_matches_python_version(self):
        """The version compiled into the C++ code must match the metadata.

        setup.py passes [project] version to CMake as BOOZ_XFORM_VERSION; the
        C++ side reports it via booz_xform::version, which is what
        write_boozmn() records in the output file.  A build that did not go
        through setup.py leaves it at the placeholder, which is not an error
        here.
        """
        cpp_version = bx._booz_xform.__version__
        if cpp_version in ("", "dev", "development version"):
            self.skipTest(
                "extension was built without BOOZ_XFORM_VERSION "
                "(not via setup.py)"
            )
        self.assertEqual(cpp_version, bx.__version__)


if __name__ == "__main__":
    unittest.main()
