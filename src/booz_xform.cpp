#include <pybind11/pybind11.h>

// namespace py = pybind11;

int add(int i, int j) {
  return i + j;
}

PYBIND11_MODULE(booz_xform, m) {
  m.doc() = "My first pybind11 module";
  m.def("add", &add, "A little function that adds two numbers.");
}
