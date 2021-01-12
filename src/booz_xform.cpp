#include <iostream>
#include <pybind11/pybind11.h>
#include "booz_xform.hpp"

// namespace py = pybind11;

using namespace booz_xform;

Booz_xform::Booz_xform() {
  std::cout << "Hello world from constructor" << std::endl;
}

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
}
