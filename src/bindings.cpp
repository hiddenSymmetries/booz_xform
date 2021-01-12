#include <iostream>
#include <pybind11/pybind11.h>
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
  m.doc() = "My first pybind11 module";
  m.def("add", &add, "A little function that adds two numbers.");
  py::class_<Booz_xform>(m, "Booz_xform")
    .def(py::init())
    .def("read_boozmn", &Booz_xform::read_boozmn)
    .def_readwrite("mboz", &Booz_xform::mboz)
    .def_readwrite("nboz", &Booz_xform::nboz)
    .def_readwrite("verbose", &Booz_xform::verbose)
    .def("testfunc1", &Booz_xform::testfunc1)
    .def("testfunc2", &Booz_xform::testfunc2);
    
}
