"""booz_xform: transformation to Boozer coordinates."""

from importlib.metadata import PackageNotFoundError, version as _version

# The single source of truth for the version is [project] version in
# pyproject.toml.  Reading it back out of the installed distribution metadata
# means booz_xform.__version__ can never disagree with what pip and PyPI
# report.  See dev/RELEASING.md.
try:
    __version__ = _version("booz_xform")
except PackageNotFoundError:
    # Imported from a source tree with no installed distribution (e.g. an
    # in-place `setup.py build_ext` build).  Fall back to the version the C++
    # extension was compiled with, which is "development version" unless that
    # build went through setup.py.
    from ._booz_xform import __version__

from ._booz_xform import Booz_xform, omp_max_threads
from .plots import surfplot, symplot, modeplot, wireplot
