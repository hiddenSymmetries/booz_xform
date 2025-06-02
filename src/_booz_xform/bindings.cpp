#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "booz_xform.hpp"

namespace py = pybind11;

using namespace booz_xform;
using namespace pybind11::literals; // For documenting function arguments

PYBIND11_MODULE(_booz_xform, m) {
  m.doc() = "Transformation to Boozer coordinates";

  py::class_<Booz_xform>(m, "Booz_xform", R"(
A class holding the input and output information for a
transformation to Boozer coordinates)")

    .def(py::init())

    .def("read_wout", &Booz_xform::read_wout, R"(
Read input information from a VMEC ``wout*.nc`` file.

:param filename: The full name of the file to load.
:param flux: If true, the poloidal flux is read. (Defaults to false)
)",
	 "filename"_a,
     "flux"_a=false)

    .def("init_from_vmec", &Booz_xform::init_from_vmec, R"(
Handle conversion of the radial grids used in vmec to to the radial
grid used by booz_xform. This involves truncation of the leading 0 for
quantities on vmec's half grid, and interpolation for quantities on
vmec's full grid. This function is called by the higher level function
:func:`~booz_xform.Booz_xform.read_wout()`, so usually users do not
need to call this function directly.

:param ns: The number of radial vmec surfaces.
:param iotas: Iota on vmec's half grid.
:param rmnc0: Vmec's original rmnc array, on the full grid.
:param rmns0: Vmec's original rmns array, on the full grid.
  For stellarator-symmetric configurations this array is ignored and need not
  be specified.
:param zmnc0: Vmec's original zmnc array, on the full grid.
  For stellarator-symmetric configurations this array is ignored and need not
  be specified.
:param zmns0: Vmec's original zmns array, on the full grid.
:param lmnc0: Vmec's original lmnc array, on the half grid.
  For stellarator-symmetric configurations this array is ignored and need not
  be specified.
:param lmns0: Vmec's original lmns array, on the half grid.
:param bmnc0: Vmec's original bmnc array, on the half grid.
:param bmns0: Vmec's original bmns array, on the half grid.
  For stellarator-symmetric configurations this array is ignored and need not
  be specified.
:param bsubumnc0: Vmec's original bsubumnc array, on the half grid.
:param bsubumns0: Vmec's original bsubumns array, on the half grid.
  For stellarator-symmetric configurations this array is ignored and need not
  be specified.
:param bsubvmnc0: Vmec's original bsubvmnc array, on the half grid.
:param bsubvmns0: Vmec's original bsubvmns array, on the half grid.
  For stellarator-symmetric configurations this array is ignored and need not
  be specified.
:param phip: Phip on vmec's full grid. (Defaults to empty Vector)
:param chi: chi on vmec's full grid. (Defaults to empty Vector)
:param pres: pressure on vmec's full grid. (Defaults to empty Vector)
:param phi: phi on vmec's full grid. (Defaults to empty Vector)
)",
	 "ns"_a,
	 "iotas"_a,
	 "rmnc0"_a,
	 "rmns0"_a,
	 "zmnc0"_a,
	 "zmns0"_a,
	 "lmnc0"_a,
	 "lmns0"_a,
	 "bmnc0"_a,
	 "bmns0"_a,
	 "bsubumnc0"_a,
	 "bsubumns0"_a,
	 "bsubvmnc0"_a,
	 "bsubvmns0"_a,
     "phips"_a=defaultInitPtr,
     "chi"_a=defaultInitPtr,
     "pres"_a=defaultInitPtr,
     "phi"_a=defaultInitPtr)

    .def("run", &Booz_xform::run, R"(
Run the transformation to Boozer coordinates on all the selected
flux surfaces. The input arrays should be initialized before this step.)")

    .def("write_boozmn", &Booz_xform::write_boozmn, R"(
Write the results of the transformation to a classic ``boozmn_*.nc``
NetCDF file. This function should only be called after
:func:`booz_xform.Booz_xform.run`.

:param filename: The full name of the file to save.)",
	 "filename"_a)

    .def("read_boozmn", &Booz_xform::read_boozmn, R"(
Read in the results of an earlier transformation from a classic
``boozmn_*.nc`` NetCDF file.

:param filename: The full name of the file to load.)",
	 "filename"_a)

    // End of functions. Now come properties that are inputs.

    .def_readwrite("verbose", &Booz_xform::verbose, R"(
(int, input) Set this to 0 for no output to stdout, 1 for some
output, 2 for lots of output)")

    .def_readwrite("asym", &Booz_xform::asym, R"(
(bool, input) True if the configuration is not
stellarator-symmetric, false if the configuration is
stellarator-symmetric.)")

    .def_readwrite("nfp", &Booz_xform::nfp, R"(
(int, input) Number of field periods, i.e. the discrete toroidal
rotation symmetry.)")

    .def_readwrite("mpol", &Booz_xform::mpol, R"(
(int, input) Maximum poloidal mode number for the input arrays rmnc,
rmns, zmnc, zmns, lmnc, and lmns.)")

    .def_readwrite("ntor", &Booz_xform::ntor, R"(
(int, input) Maximum toroidal mode number (divided by nfp) for the
input arrays rmnc, rmns, zmnc, zmns, lmnc, and lmns.)")

    .def_readwrite("mnmax", &Booz_xform::mnmax, R"(
(int, input) Number of Fourier modes for the input arrays rmnc,
rmns, zmnc, zmns, lmnc, and lmns.)")

    .def_readwrite("mpol_nyq", &Booz_xform::mpol_nyq, R"(
(int, input) Maximum poloidal mode number for the input arrays bmnc,
bmns, bsubumnc, bsubumns, bsubvmnc, and bsubvmns.)")

    .def_readwrite("ntor_nyq", &Booz_xform::ntor_nyq, R"(
(int, input) Maximum toroidal mode number (divided by nfp) for the
input arrays bmnc, bmns, bsubumnc, bsubumns, bsubvmnc, and bsubvmns.)")

    .def_readwrite("mnmax_nyq", &Booz_xform::mnmax_nyq, R"(
(int, input) Total number of Fourier modes for the input arrays
bmnc, bmns, bsubumnc, bsubumns, bsubvmnc, and bsubvmns.)")

    .def_readwrite("xm", &Booz_xform::xm, R"(
(1D integer array of length mnmax, input) The poloidal Fourier mode
numbers for the input arrays rmnc, rmns, zmnc, zmns, lmnc, and lmns.)")

    .def_readwrite("xn", &Booz_xform::xn, R"(
(1D integer array of length mnmax, input) The toroidal Fourier mode
numbers for the input arrays rmnc, rmns, zmnc, zmns, lmnc, and
lmns. The values should all be integer multiples of nfp.)")

    .def_readwrite("xm_nyq", &Booz_xform::xm_nyq, R"(
(1D integer array of length mnmax_nyq, input) The poloidal Fourier
mode numbers for the input arrays bmnc, bmns, bsubumnc, bsubumns,
bsubvmnc, and bsubvmns.)")

    .def_readwrite("xn_nyq", &Booz_xform::xn_nyq, R"(
(1D integer array of length mnmax_nyq, input) The toroidal Fourier
mode numbers for the input arrays bmnc, bmns, bsubumnc, bsubumns,
bsubvmnc, and bsubvmns. The values should all be integer multiples of
nfp.)")

    .def_readwrite("ns_in", &Booz_xform::ns_in, R"(
(int, input) Number of flux surfaces on which the input data is
supplied. The transformation to Boozer coordinates is not necessarily
run on all of these surfaces, only the ones indicated by
compute_surfs.)")

    .def_readwrite("s_in", &Booz_xform::s_in, R"(
(1D float array of length ns_in, input) The values of normalized
toroidal flux s for the input data. Here, s is the toroidal flux
normalized to the value at the plasma boundary. This information is
not needed for the coordinate transformation itself, but is useful for
plotting output.)")

    .def_readwrite("iota", &Booz_xform::iota, R"(
(1D float array of length ns_in, input) Rotational transform on the
input radial surfaces.)")

    .def_readwrite("phi", &Booz_xform::phi, R"(
(1D float array of length ns_in, input) Toroidal flux normalized by 2*pi,
evaluated on all the magnetic surfaces for which input data was provided.)")

    .def_readwrite("chi", &Booz_xform::chi, R"(
(1D float array of length ns_in, input) Poloidal flux normalized by 2*pi,
evaluated on all the magnetic surfaces for which input data was provided.)")

    .def_readwrite("phip", &Booz_xform::phip, R"(
(1D float array of length ns_in, input) Toroidal flux normalized by 2*pi,
evaluated on all the magnetic surfaces for which input data was provided.)")

    .def_readwrite("pres", &Booz_xform::pres, R"(
(1D float array of length ns_in, input) Pressure on full vmec grid.)")

    .def_readwrite("rmnc", &Booz_xform::rmnc, R"(
(2D float array of size mnmax x ns_in, input) cos(m * theta_0 - n *
zeta_0) Fourier modes of the major radius R.)")

    .def_readwrite("rmns", &Booz_xform::rmns, R"(
(2D float array of size mnmax x ns_in, input) sin(m * theta_0 - n *
zeta_0) Fourier modes of the major radius R of the flux surfaces. For
stellarator-symmetric configurations, this array is not used and need
not be specified.)")

    .def_readwrite("zmnc", &Booz_xform::zmnc, R"(
(2D float array of size mnmax x ns_in, input) cos(m * theta_0 - n *
zeta_0) Fourier modes of the Cartesian coordinate Z of the flux
surfaces. For stellarator-symmetric configurations, this array is not
used and need not be specified.)")

    .def_readwrite("zmns", &Booz_xform::zmns, R"(
(2D float array of size mnmax x ns_in, input) sin(m * theta_0 - n *
zeta_0) Fourier modes of the Cartesian coordinate Z of the flux
surfaces.)")

    .def_readwrite("lmnc", &Booz_xform::lmnc, R"(
(2D float array of size mnmax_nyq x ns_in, input) cos(m * theta_0 - n
* zeta_0) Fourier modes of lambda = theta^* - theta_0, the difference
between the original poloidal angle theta_0 and the straight field
line poloidal angle theta^*. For stellarator-symmetric configurations,
this array is not used and need not be specified.)")

    .def_readwrite("lmns", &Booz_xform::lmns, R"(
(2D float array of size mnmax_nyq x ns_in, input) sin(m * theta_0 - n
* zeta_0) Fourier modes of lambda = theta^* - theta_0, the difference
between the original poloidal angle theta_0 and the straight field
line poloidal angle theta^*.)")

    .def_readwrite("bmnc", &Booz_xform::bmnc, R"(
(2D float array of size mnmax_nyq x ns_in, input) cos(m * theta_0 - n *
zeta_0) Fourier modes of the magnetic field strength B.)")

    .def_readwrite("bmns", &Booz_xform::bmns, R"(
(2D float array of size mnmax_nyq x ns_in, input) sin(m * theta_0 - n
* zeta_0) Fourier modes of the magnetic field strength B.  For
stellarator-symmetric configurations, this array is not used and need
not be specified.)")

    .def_readwrite("bsubumnc", &Booz_xform::bsubumnc, R"(
(2D float array of size mnmax_nyq x ns_in, input) cos(m * theta_0 -
n * zeta_0) Fourier modes of B dot (d r / d theta_0) where r is the
position vector.)")

    .def_readwrite("bsubumns", &Booz_xform::bsubumns, R"(
(2D float array of size mnmax_nyq x ns_in, input) sin(m * theta_0 -
n * zeta_0) Fourier modes of B dot (d r / d theta_0) where r is the
position vector.  For stellarator-symmetric configurations, this array
is not used and need not be specified.)")

    .def_readwrite("bsubvmnc", &Booz_xform::bsubvmnc, R"(
(2D float array of size mnmax_nyq x ns_in, input) cos(m * theta_0 -
n * zeta_0) Fourier modes of B dot (d r / d zeta_0) where r is the
position vector.)")

    .def_readwrite("bsubvmns", &Booz_xform::bsubvmns, R"(
(2D float array of size mnmax_nyq x ns_in, input) sin(m * theta_0 -
n * zeta_0) Fourier modes of B dot (d r / d zeta_0) where r is the
position vector.  For stellarator-symmetric configurations, this array
is not used and need not be specified.)")

    .def_readwrite("mboz", &Booz_xform::mboz, R"(
(int, input) Maximum poloidal mode number for representing output
quantities in Boozer coordinates.)")

    .def_readwrite("nboz", &Booz_xform::nboz, R"(
(int, input) Maximum toroidal mode number (divided by nfp) for
representing output quantities in Boozer coordinates. For example, if
nboz=2 and nfp=10, the toroidal modes used will be n=-20, -10, 0, 10,
20.)")

    .def_readwrite("compute_surfs", &Booz_xform::compute_surfs, R"(
(1D integer array, input) Indices of ns_in-sized radial arrays,
specifying the flux surfaces for which the transformation to Boozer
coordinates will be performed. All values should be >= 0 and < ns_in.
The array compute_surfs is similar to the array jlist in the earlier
fortran booz_xform program, with compute_surfs = jlist - 2.)")

    .def_readwrite("aspect", &Booz_xform::aspect, R"(
(float, input) The aspect ratio of the configuration. This value is
not used for anything by booz_xform, and does not need to be set. It
is provided as a means to pass this value from the input equilibrium
to booz_xform output files.)")

    .def_readwrite("toroidal_flux", &Booz_xform::toroidal_flux, R"(
(float, input) The boundary toroidal flux of the configuration. This
value is not used for anything by booz_xform, and does not need to be
set. It is provided as a means to pass this value from the input
equilibrium to booz_xform output files.)")

    // End of inputs. Now come the outputs.

    .def_readonly("ns_b", &Booz_xform::ns_b, R"(
(int, output) Number of flux surfaces on which output data are
available.)")

    .def_readonly("s_b", &Booz_xform::s_b, R"(
(1D float array of length ns_b, output) Values of normalized
toroidal flux s defining the magnetic surfaces for the output data.)")

    .def_readonly("mnboz", &Booz_xform::mnboz, R"(
(int, output) Total number of Fourier modes for output data.)")

    .def_readonly("xm_b", &Booz_xform::xm_b, R"(
(1D int array of length mnboz, output) Poloidal mode numbers for the
output data.)")

    .def_readonly("xn_b", &Booz_xform::xn_b, R"(
(1D int array of length mnboz, output) Toroidal mode numbers for the
output data. These values are all integer multiples of nfp.)")

    .def_readonly("bmnc_b", &Booz_xform::bmnc_b, R"(
(2D float array of size mnboz x ns_b, output) cos(m * theta_B - n *
zeta_B) Fourier modes of the magnetic field strength in Boozer
coordinates.)")

    .def_readonly("bmns_b", &Booz_xform::bmns_b, R"(
(2D float array of size mnboz x ns_b, output) sin(m * theta_B - n *
zeta_B) Fourier modes of the magnetic field strength in Boozer
coordinates. If the configuration is stellarator-symmetric, this
quantity is zero so the array will have size 0 x 0.)")

    .def_readonly("gmnc_b", &Booz_xform::gmnc_b, R"(
(2D float array of size mnboz x ns_b, output) cos(m * theta_B - n *
zeta_B) Fourier modes (with respect to Boozer coordinates) of the
Jacobian of (psi, theta_B, zeta_B) coordinates.)")

    .def_readonly("gmns_b", &Booz_xform::gmns_b, R"(
(2D float array of size mnboz x ns_b, output) sin(m * theta_B - n *
zeta_B) Fourier modes (with respect to Boozer coordinates) of the
Jacobian of (psi, theta_B, zeta_B) coordinates. If the configuration is
stellarator-symmetric, this quantity is zero so this array will have
size 0 x 0.)")

    .def_readonly("rmnc_b", &Booz_xform::rmnc_b, R"(
(2D float array of size mnboz x ns_b, output) cos(m * theta_B - n *
zeta_B) Fourier modes (with respect to Boozer coordinates) of the
major radius R of the flux surfaces.)")

    .def_readonly("rmns_b", &Booz_xform::rmns_b, R"(
(2D float array of size mnboz x ns_b, output) sin(m * theta_B - n *
zeta_B) Fourier modes (with respect to Boozer coordinates) of the
major radius R of the flux surfaces. If the configuration is
stellarator-symmetric, this quantity is zero so this array will have
size 0 x 0.)")

    .def_readonly("zmnc_b", &Booz_xform::zmnc_b, R"(
(2D float array of size mnboz x ns_b, output) cos(m * theta_B - n *
zeta_B) Fourier modes (with respect to Boozer coordinates) of the
Cartesian coordinate Z of the flux surfaces. If the configuration is
stellarator-symmetric, this quantity is zero so array will have size 0
x 0.)")

    .def_readonly("zmns_b", &Booz_xform::zmns_b, R"(
(2D float array of size mnboz x ns_b, output) sin(m * theta_B - n *
zeta_B) Fourier modes (with respect to Boozer coordinates) of the
Cartesian coordinate Z of the flux surfaces.)")

    .def_readonly("numnc_b", &Booz_xform::numnc_b, R"(
(2D float array of size mnboz x ns_b, output) cos(m * theta_B - n *
zeta_B) Fourier modes (with respect to Boozer coordinates) of the
toroidal angle difference nu = zeta_B - zeta_0. If the configuration
is stellarator-symmetric, this quantity is zero so this array will
have size 0 x 0.)")

    .def_readonly("numns_b", &Booz_xform::numns_b, R"(
(2D float array of size mnboz x ns_b, output) sin(m * theta_B - n *
zeta_B) Fourier modes (with respect to Boozer coordinates) of the
toroidal angle difference nu = zeta_B - zeta_0.)")

    .def_readonly("Boozer_G", &Booz_xform::Boozer_G, R"(
(1D float array of length ns_b, output) Coefficient of grad zeta_B
in the covariant representation of the magnetic field vector B in
Boozer coordinates, evaluated on the magnetic surfaces used for output
quantities.)")

    .def_readonly("Boozer_G_all", &Booz_xform::Boozer_G_all, R"(
(1D float array of length ns_in, output) Coefficient of grad zeta_B
in the covariant representation of the magnetic field vector B in
Boozer coordinates, evaluated on all the magnetic surfaces for which
input data was provided.)")

    .def_readonly("Boozer_I", &Booz_xform::Boozer_I, R"(
(1D float array of length ns_b, output) Coefficient of grad theta_B
in the covariant representation of the magnetic field vector B in
Boozer coordinates, evaluated on the magnetic surfaces used for output
quantities.)")

    .def_readonly("Boozer_I_all", &Booz_xform::Boozer_I_all, R"(
(1D float array of length ns_in, output) Coefficient of grad theta_B
in the covariant representation of the magnetic field vector B in
Boozer coordinates, evaluated on all the magnetic surfaces for which
input data was provided.)")

    ;

  // Trick for passing version number from setup.py, from
  // https://github.com/pybind/cmake_example/blob/master/src/main.cpp
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}

// https://github.com/pybind/pybind11/issues/2271

// https://github.com/pybind/pybind11/issues/1042

// Possible ways to wrap a custom Matrix class, without using Eigen:
// https://github.com/tdegeus/pybind11_examples/tree/master/09_numpy_cpp-custom-matrix
// https://stackoverflow.com/questions/42645228/cast-numpy-array-to-from-custom-c-matrix-class-using-pybind11
