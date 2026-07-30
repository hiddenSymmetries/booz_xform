"""Command-line driver for booz_xform, installed as the ``xbooz_xform`` script.

This is a Python port of ``booz_xform::driver()`` in
``src/_booz_xform/driver.cpp``: same input file format, same output file, same
console output.  It exists because CMake's standalone ``xbooz_xform`` executable
cannot be shipped inside a binary wheel -- ``auditwheel``/``delocate`` rewrite
library paths relative to the package directory, which the executable does not
end up next to once installed.  Building on the extension module instead means
the command works identically for wheel installs, source installs, and source
checkouts.

Source checkouts still get the C++ executable from ``cmake --build . --target
xbooz_xform``; the two are interchangeable.
"""

import sys

from . import Booz_xform, __version__

USAGE = """Usage:  xbooz_xform <inputfile>

Where <inputfile> is a text file with 2 or 3 lines, formatted as follows:

<mboz> <nboz>
<extension>
<compute_surfs>

Here,

<mboz> is the maximum poloidal Fourier mode number to include in the Boozer representation.

<nboz> is the maximum toroidal Fourier mode number to include in the Boozer representation.

<extension> is the VMEC extension from the VMEC input file to process.
For instance, to process the file wout_li383_1.4m.nc, <extension> should be li383_1.4m
Note that the wout file to load must be in the current working directory.

<compute_surfs> is a list of integers giving the surfaces of the VMEC input file to
transform. Note that the transformation is only performed on half-grid surfaces.
The first half-grid surface has index 0. The outermost available surface has index NS-2,
where NS is the VMEC input parameter. You can omit <compute_surfs> from the input file,
in which case the transformation will be performed on all half-grid surfaces."""


def _parse_input_file(filename):
    """Read an in_booz-style input file.

    Returns (mboz, nboz, extension, compute_surfs).  Whitespace and line breaks
    are handled the same way as the C++ ``>>`` extraction operators, i.e. all
    tokens are read in order and newlines are not significant.
    """
    with open(filename) as f:
        tokens = f.read().split()

    if len(tokens) < 3:
        raise RuntimeError(
            "Unable to read mboz, nboz, and extension from " + filename
        )

    try:
        mboz, nboz = int(tokens[0]), int(tokens[1])
    except ValueError:
        raise RuntimeError("Unable to read mboz and nboz")
    print("Read mboz = {}, nboz = {}".format(mboz, nboz))

    extension = tokens[2]
    print("Read extension =", extension)

    # As in driver.cpp, reading stops at the first token that is not an integer
    # rather than being treated as an error.
    compute_surfs = []
    for token in tokens[3:]:
        try:
            value = int(token)
        except ValueError:
            print("Unable to read more.")
            break
        print("Read value", value)
        compute_surfs.append(value)

    return mboz, nboz, extension, compute_surfs


def main(argv=None):
    """Entry point for the ``xbooz_xform`` command."""
    if argv is None:
        argv = sys.argv[1:]

    if len(argv) != 1:
        print(USAGE)
        return 0

    print("This is xbooz_xform " + __version__)

    mboz, nboz, extension, compute_surfs = _parse_input_file(argv[0])

    b = Booz_xform()
    b.read_wout("wout_" + extension + ".nc")

    b.mboz = max(b.mboz, mboz)
    b.nboz = max(b.nboz, nboz)

    # If no surfaces are specified, we do not need to modify the default
    # compute_surfs set up by read_wout().
    if compute_surfs:
        b.compute_surfs = sorted(compute_surfs)
    else:
        print("No compute_surfs specified, so including all half-grid surfaces.")

    print("About to run transformation with compute_surfs =",
          " ".join(str(s) for s in b.compute_surfs))

    b.run()
    b.write_boozmn("boozmn_" + extension + ".nc")

    print("Good bye.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
