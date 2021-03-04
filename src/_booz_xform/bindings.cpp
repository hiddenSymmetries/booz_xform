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
)",
	 "filename"_a)
    
    .def("init_from_vmec", &Booz_xform::init_from_vmec, R"(
Handle radial interpolation of vmec results to the half grid.
)")
    
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
    
    .def_readwrite("verbose", &Booz_xform::verbose, R"(
Set this to 0 for no output to stdout, 1 for some output, 2 for lots
of output)")
    
    .def_readwrite("asym", &Booz_xform::asym) 
    .def_readwrite("nfp", &Booz_xform::nfp)
    .def_readwrite("s_in", &Booz_xform::s_in)
    .def_readwrite("mpol", &Booz_xform::mpol, R"(
Maximum poloidal mode number for the input arrays rmnc, rmns, zmnc, and zmns)")
    
    .def_readwrite("ntor", &Booz_xform::ntor)
    .def_readwrite("mnmax", &Booz_xform::mnmax)
    .def_readwrite("mpol_nyq", &Booz_xform::mpol_nyq)
    .def_readwrite("ntor_nyq", &Booz_xform::ntor_nyq)
    .def_readwrite("mnmax_nyq", &Booz_xform::mnmax_nyq)
    .def_readwrite("xm", &Booz_xform::xm)
    .def_readwrite("xn", &Booz_xform::xn)
    .def_readwrite("xm_nyq", &Booz_xform::xm_nyq)
    .def_readwrite("xn_nyq", &Booz_xform::xn_nyq)
    .def_readwrite("ns_in", &Booz_xform::ns_in)
    .def_readwrite("iota", &Booz_xform::iota)
    .def_readwrite("rmnc", &Booz_xform::rmnc)
    .def_readwrite("rmns", &Booz_xform::rmns)
    .def_readwrite("zmnc", &Booz_xform::zmnc)
    .def_readwrite("zmns", &Booz_xform::zmns)
    .def_readwrite("lmnc", &Booz_xform::lmnc)
    .def_readwrite("lmns", &Booz_xform::lmns)
    .def_readwrite("bmnc", &Booz_xform::bmnc)
    .def_readwrite("bmns", &Booz_xform::bmns)
    .def_readwrite("bsubumnc", &Booz_xform::bsubumnc)
    .def_readwrite("bsubumns", &Booz_xform::bsubumns)
    .def_readwrite("bsubvmnc", &Booz_xform::bsubvmnc)
    .def_readwrite("bsubvmns", &Booz_xform::bsubvmns)
    .def_readwrite("mboz", &Booz_xform::mboz)
    .def_readwrite("nboz", &Booz_xform::nboz)
    .def_readwrite("compute_surfs", &Booz_xform::compute_surfs)
    // End of inputs. Now come the outputs.
    .def_readonly("ns_b", &Booz_xform::ns_b)
    .def_readonly("s_b", &Booz_xform::s_b)
    .def_readonly("mnboz", &Booz_xform::mnboz)
    .def_readonly("xm_b", &Booz_xform::xm_b)
    .def_readonly("xn_b", &Booz_xform::xn_b)
    .def_readonly("bmnc_b", &Booz_xform::bmnc_b)
    .def_readonly("bmns_b", &Booz_xform::bmns_b)
    .def_readonly("gmnc_b", &Booz_xform::gmnc_b)
    .def_readonly("gmns_b", &Booz_xform::gmns_b)
    .def_readonly("rmnc_b", &Booz_xform::rmnc_b)
    .def_readonly("rmns_b", &Booz_xform::rmns_b)
    .def_readonly("zmnc_b", &Booz_xform::zmnc_b)
    .def_readonly("zmns_b", &Booz_xform::zmns_b)
    .def_readonly("pmnc_b", &Booz_xform::pmnc_b)
    .def_readonly("pmns_b", &Booz_xform::pmns_b);
    
}

// https://github.com/pybind/pybind11/issues/2271

// https://github.com/pybind/pybind11/issues/1042

// Possible ways to wrap a custom Matrix class, without using Eigen:
// https://github.com/tdegeus/pybind11_examples/tree/master/09_numpy_cpp-custom-matrix
// https://stackoverflow.com/questions/42645228/cast-numpy-array-to-from-custom-c-matrix-class-using-pybind11
