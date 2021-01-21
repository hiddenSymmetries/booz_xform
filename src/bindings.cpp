#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "booz_xform.hpp"

namespace py = pybind11;

using namespace booz_xform;

int add(int i, int j) {
  Vector v;

  v.resize(5, 3.14);
  //std::cout << "Hello world! ";
  std::cout << v << std::endl;

  Booz_xform booz;
  booz.read_boozmn("/Users/mattland/pyBooz_xform/boozmn_n3are_R7.75B5.7.nc");
  
  return i + j;
}

PYBIND11_MODULE(booz_xform, m) {
  m.doc() = "Transformation to Boozer coordinates";
  m.def("add", &add, "A little function that adds two numbers.");
  /*
  py::class_<Vector>(m, "Vector", py::buffer_protocol())
    .def_buffer([](Vector &vec) -> py::buffer_info {
      return py::buffer_info(
			     &vec[0], // Pointer to buffer
			     sizeof(boozfloat), // Size of one scalar
			     py::format_descriptor<boozfloat>::format(), // Python struct-style format descriptor
			     1, // Number of dimensions
			     { vec.size() }, // Buffer dimensions
			     { sizeof(boozfloat)} // Strides (in bytes) for each index
			     );
    });
  */
  py::class_<Matrix>(m, "Matrix", py::buffer_protocol())
    .def_buffer([](Matrix &mat) -> py::buffer_info {
      // This ".def_buffer" section handles conversion from C++ to python
      return py::buffer_info(
			     &mat[0], // Pointer to buffer
			     sizeof(boozfloat), // Size of one scalar
			     py::format_descriptor<boozfloat>::format(), // Python struct-style format descriptor
			     2, // Number of dimensions
			     { mat.nrows(), mat.ncols() }, // Buffer dimensions
			     { sizeof(boozfloat), // Strides (in bytes) for each index
			       sizeof(boozfloat) * mat.nrows()}
			     );
    })
    .def(py::init([](py::buffer b) {
      // This ".def" section handles conversion from python to C++.
      // This section is based on the pybind11 page
      // https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html
      // specifically the example following "To create a C++ function that can take a Python buffer object as an argument"

      // Request a buffer descriptor from Python
      py::buffer_info info = b.request();

      // Some sanity checks ...
      if (info.format != py::format_descriptor<boozfloat>::format())
	throw std::runtime_error("Incompatible format: expected a double array!");

      if (info.ndim != 2)
	throw std::runtime_error("Incompatible buffer dimension!");

      Matrix mat(info.shape[0], info.shape[1]);
      boozfloat * data = static_cast<boozfloat *>(info.ptr);
      // Strides are in bytes. For convenience, scale by sizeof(boozfloat)
      info.strides[0] /= sizeof(boozfloat);
      info.strides[1] /= sizeof(boozfloat);
      // Copy data. This approach is slow but it should handle both
      // row-major and col-major numpy arrays with any strides.
      for (int j = 0; j < info.shape[0]; j++) {
	for (int k = 0; k < info.shape[1]; k++) {
	  mat(j, k) = data[j * info.strides[0] + k * info.strides[1]];
	}
      }

      return mat;
    }));
  
  py::class_<Booz_xform>(m, "Booz_xform")
    .def(py::init())
    .def("read_boozmn", &Booz_xform::read_boozmn)
    .def("read_wout", &Booz_xform::read_wout)
    .def("run", &Booz_xform::run)
    .def("write_boozmn", &Booz_xform::write_boozmn)
    .def_readwrite("mboz", &Booz_xform::mboz)
    .def_readwrite("nboz", &Booz_xform::nboz)
    .def_readwrite("xmb", &Booz_xform::xmb)
    .def_readwrite("xnb", &Booz_xform::xnb)
    .def_readwrite("mpol", &Booz_xform::mpol)
    .def_readwrite("ntor", &Booz_xform::ntor)
    .def_readwrite("ns", &Booz_xform::ns)
    .def_readwrite("iotas", &Booz_xform::iotas)
    .def_readwrite("rmnc", &Booz_xform::rmnc)
    .def_readwrite("verbose", &Booz_xform::verbose)
    .def_readwrite("jlist", &Booz_xform::jlist)
    .def_readwrite("bmnc_b", &Booz_xform::bmnc_b)
    .def_readwrite("bmns_b", &Booz_xform::bmns_b)
    .def_readwrite("gmnc_b", &Booz_xform::gmnc_b)
    .def_readwrite("gmns_b", &Booz_xform::gmns_b)
    .def_readwrite("rmnc_b", &Booz_xform::rmnc_b)
    .def_readwrite("rmns_b", &Booz_xform::rmns_b)
    .def_readwrite("zmnc_b", &Booz_xform::zmnc_b)
    .def_readwrite("zmns_b", &Booz_xform::zmns_b)
    .def_readwrite("pmnc_b", &Booz_xform::pmnc_b)
    .def_readwrite("pmns_b", &Booz_xform::pmns_b);
    
}

// https://github.com/pybind/pybind11/issues/2271

// https://github.com/pybind/pybind11/issues/1042
