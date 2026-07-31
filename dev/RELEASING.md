# Releasing booz_xform

Maintainer notes. Users do not need this file.

## How the version is wired

There is exactly **one** version literal in the repository:

```toml
# pyproject.toml
[project]
version = "0.1.0"
```

Everything else derives from it:

| Consumer | How it gets the version |
|---|---|
| wheel / sdist metadata, PyPI | setuptools reads `[project] version` directly |
| `booz_xform.__version__` | `importlib.metadata.version("booz_xform")`, i.e. read back out of the installed distribution |
| `booz_xform::version` in C++, i.e. `_booz_xform.__version__` | `setup.py` passes it to CMake as `-DBOOZ_XFORM_VERSION`, which becomes the `VERSION_INFO` preprocessor macro |
| the `version` attribute in `boozmn_*.nc` output files | `booz_xform::version` |
| what `xbooz_xform` prints on startup | `booz_xform.__version__` |
| `docs/source/conf.py` | `importlib.metadata.version("booz_xform")` |
| the git tag | written by `tbump`, then re-checked by CI at release time |

`[tool.tbump.version] current` in `pyproject.toml` duplicates the literal, but
`tbump` itself keeps the two in step and errors out if they disagree.

Two independent guards keep this honest:

- `tests/test_version.py` fails if `booz_xform.__version__`, the installed
  distribution metadata, or the version compiled into the C++ code disagree.
- The `check-version` job in `.github/workflows/release.yml` refuses to build or
  publish anything if the release tag does not match `pyproject.toml`.

## Cutting a release

The `release.yml` workflow is used both for test releases to TestPyPI and for
real releases to PyPI. Triggered manually (`workflow_dispatch`), it uploads to
TestPyPI; triggered by publishing a GitHub Release, it uploads to PyPI.
Publishing uses **PyPI Trusted Publishing (OIDC)**, so there is no API token in
GitHub secrets and nothing to rotate.

### 1. Rehearse against TestPyPI

From the GitHub Actions tab, run the **Release** workflow via *Run workflow*
(`workflow_dispatch`) with *Upload to TestPyPI* checked. This builds the full
wheel matrix plus the sdist, smoke-tests every wheel, and uploads to TestPyPI —
without creating a tag or touching real PyPI.

Then install the result on a clean machine and confirm it works:

```bash
python -m venv /tmp/bx && /tmp/bx/bin/pip install -U pip
/tmp/bx/bin/pip install --index-url https://test.pypi.org/simple/ \
    --extra-index-url https://pypi.org/simple/ booz_xform
/tmp/bx/bin/python -c "import booz_xform as bx; print(bx.__version__); print(bx.omp_max_threads())"
```

The `--extra-index-url` is required: TestPyPI does not mirror numpy.

If you want to test the release machinery itself, use a release-candidate
version (`tbump 0.1.0rc1`) so a botched rehearsal does not burn the next clean
version number — **a version can never be re-uploaded to PyPI, even after
deletion.**

### 2. Bump and tag

Version bumps are done with `tbump`, which is pip-installed as part of
booz_xform's `[dev]` optional dependencies:

```bash
git checkout main && git pull      # tbump pushes to the current branch's upstream
tbump 0.2.0
```

`tbump` shows exactly what it plans to do and asks for confirmation. It rewrites
the version in `pyproject.toml`, runs the test suite as a pre-commit hook,
commits as `Release 0.2.0`, creates the annotated tag `v0.2.0`, and pushes the
branch and tag together.

Run it from a branch that has an upstream. On a local-only branch, `tbump`
renders an incomplete `git push` command and the push fails.

### 3. Publish the GitHub Release

Go to *Releases → Draft a new release*, choose the existing `v0.2.0` tag, use
*Generate release notes*, edit as needed, and **Publish**.

Publishing triggers the release workflow. Once all the wheels are built
(~30 minutes), GitHub asks you to approve the final upload to PyPI (if you set
up the `pypi` environment with a required reviewer, as recommended below). After
you approve, the wheels and sdist go to PyPI.

### 4. Verify

```bash
pip download --no-deps --only-binary :all: booz_xform==0.2.0
```

Check that an installed wheel is parallel and self-contained. `omp_max_threads()`
returns 1 both on a single-core machine and when the wheel was built without
OpenMP, so on macOS also confirm that `libomp` was vendored:

```bash
python -c "import booz_xform as bx; print(bx.__version__, bx.omp_max_threads())"
python -c "import booz_xform, os; print(os.path.dirname(booz_xform.__file__))"
ls $(python -c "import booz_xform, os; print(os.path.dirname(booz_xform.__file__))")/.dylibs
```

## Failure modes and what to do

**The `check-version` job fails.** The tag and `pyproject.toml` disagree. Delete
the GitHub release *and* the tag (`git push --delete origin v0.2.0`), fix the
version with `tbump`, and re-release. Nothing was uploaded, so the version
number is still usable.

**A wheel build fails after the release is published.** Nothing has been
uploaded (the publish job runs only after the whole matrix succeeds), but the
GitHub release exists. Fix the problem, delete the release and tag, and re-cut.
Since nothing reached PyPI, you may reuse the version number.

**A wheel builds but its test step fails to install scipy/matplotlib.** The
wheels are tested with the `[test]` extra, and a brand-new CPython release
sometimes has no scipy or matplotlib wheel yet for one architecture. Either wait
for those wheels, or add a `test-skip` entry to `[tool.cibuildwheel]` in
`pyproject.toml` for that one build identifier — the wheel is still built and
published, just not smoke-tested.

**The upload partially succeeded.** PyPI uploads are effectively atomic per
file, but if you end up with some files uploaded and some not, you cannot
overwrite them — bump to the next patch version and release again.

## What is deliberately *not* built

- **Windows wheels.** Would need a MSVC-compatible NetCDF (and HDF5) build. Not
  attempted; there is no Windows CI either. Windows users fall back to the sdist,
  or better, to conda.
- **musllinux (Alpine) wheels.** Would need a separate apk toolchain for netcdf
  and hdf5.
- **PyPy wheels.** The bindings are built against the CPython C API via pybind11.
- **Free-threaded (`t`) builds.** Not exercised; the C++ side uses OpenMP and has
  no GIL discipline beyond pybind11's defaults.

The sdist covers all of these for anyone with a C++ compiler and NetCDF.

## Platform notes for the wheel builds

- **Linux** uses the `manylinux_2_28` (AlmaLinux 8) images. NetCDF is not in the
  base repos there, so `before-all` installs `epel-release` and then
  `netcdf-devel` via `dnf`. `auditwheel` vendors `libnetcdf`, `libhdf5`,
  `libgomp` and their dependencies into the wheel. arm64 builds run on native
  `ubuntu-24.04-arm` runners — no QEMU emulation.
- **macOS** uses Homebrew's `netcdf` and `libomp`, and `delocate` vendors
  `libnetcdf`, `libhdf5` and `libomp` into the wheel. `libomp` is what makes
  `find_package(OpenMP)` succeed with Apple clang, so the macOS wheels are
  OpenMP-parallel; without it they would be single-threaded. CMake does not
  search Homebrew's prefix by default on Apple Silicon, so `NETCDF_DIR` is set
  from `brew --prefix` in `[tool.cibuildwheel.macos] environment`.

  **The libomp trap.** The `libomp` formula is *keg-only*: `brew install libomp`
  does not symlink it into `/opt/homebrew`. CMake's `FindOpenMP` looks for
  `libomp.dylib` and `omp.h` on the default search paths, finds neither, and —
  because OpenMP is optional — configures a perfectly working *serial* build.
  The only trace is one easily-missed `Could NOT find OpenMP` line in a
  thousand-line log. This is what happened to the 0.1.0rc2 macOS wheels.
  `CMakeLists.txt` now runs `brew --prefix libomp` and sets `OpenMP_ROOT` from
  it, which makes both the library and the header discoverable (`OpenMP_ROOT` is
  honored by every `find_` call inside `FindOpenMP`, policy CMP0074).

  Note that a developer machine is a *bad* place to test this: if `libomp` was
  ever `brew link --force`-ed, or another formula left an `omp.h` in
  `/opt/homebrew/include`, `find_package(OpenMP)` succeeds locally for reasons
  that do not exist on a clean CI runner. To reproduce the CI environment,
  configure with the prefix hidden:

  ```bash
  cmake .. -DCMAKE_IGNORE_PATH="/opt/homebrew/lib;/opt/homebrew/include" \
      -DNETCDF_INCLUDES=/opt/homebrew/include \
      -DNETCDF_LIBRARIES=/opt/homebrew/lib/libnetcdf.dylib
  ```

  As a backstop, the wheel builds and `ci.yml` pass
  `CMAKE_ARGS=-DBOOZ_XFORM_REQUIRE_OPENMP=ON`, which turns "OpenMP not found"
  from a silent downgrade into a failed build.

  `MACOSX_DEPLOYMENT_TARGET` is set per-runner in `release.yml` (14.0 for arm64
  on `macos-14`, 15.0 for x86_64 on `macos-15-intel`). **It must be at least the
  SDK version the Homebrew bottles were built against**, or delocate will refuse
  the wheel. This is why the wheels do not support older macOS; lowering the
  floor would mean building NetCDF against an older SDK.
- **OpenMP caveat:** the wheels vendor `libgomp` (Linux) or `libomp` (macOS). If
  a user's NumPy/SciPy brings a different OpenMP runtime into the same process,
  oversubscription or (rarely) a crash is possible. This is the standard
  situation for OpenMP-using wheels.
- **Build time.** Only the `_booz_xform` CMake target is built during an install
  (see `setup.py`); the standalone `xbooz_xform` executable and the C++
  `unitTests` executable are not, which keeps the matrix from recompiling doctest
  on every platform. Those targets are still built by hand and by `ci.yml`.

## The xbooz_xform command

`pip install booz_xform` installs `xbooz_xform` as a console script backed by
`booz_xform/_cli.py`, a Python port of `booz_xform::driver()` with the same input
file format and output. The C++ executable of the same name cannot be shipped in
a wheel — `auditwheel`/`delocate` rewrite library paths relative to the package
directory, which an executable installed into `bin/` does not end up next to. A
source checkout can still build the C++ one:

```bash
mkdir build && cd build && cmake .. && make xbooz_xform
```

If `driver.cpp` changes, `_cli.py` should be updated to match, and vice versa.
`tests/test_cli.py` covers the Python one.

---

## Steps that must be done manually, once

These require account access and cannot be scripted from the repository.

- [ ] **Create the PyPI trusted publisher.** On PyPI → *Your projects* →
      `booz_xform` → *Publishing* → *Add a new publisher*:
      - Owner: `hiddenSymmetries`
      - Repository name: `booz_xform`
      - Workflow name: `release.yml`
      - Environment name: `pypi`

- [ ] **Create the TestPyPI trusted publisher.** Same values, on
      test.pypi.org, but environment name `testpypi`. If `booz_xform` does not
      exist on TestPyPI yet, use *Publishing → Add a pending publisher*.

- [ ] **Create the two GitHub environments.** Repo *Settings → Environments*:
      `pypi` and `testpypi`. The names must match the workflow exactly.
      Recommended: on the `pypi` environment, add yourself as a *required
      reviewer*. That gives you a final human confirmation between "wheels built
      successfully" and "uploaded to PyPI forever."

- [ ] **Delete the old PyPI API tokens** and the `pypi_api_token` /
      `test_pypi_api_token` repository secrets, which the removed `pypi.yml`
      workflow used. Trusted publishing replaces them.

- [ ] **First release under this system:** rehearse with `0.1.0rc1` end to end
      before cutting `0.1.0`.
