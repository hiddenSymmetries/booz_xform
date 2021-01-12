#include <iostream>
#include "booz_xform.hpp"

using namespace booz_xform;

void Booz_xform::defaults() {
  completed = false;
  verbose = 1;
  mboz = 16;
  nboz = 16;
}

Booz_xform::Booz_xform() {
  std::cout << "Hello world from constructor" << std::endl;
  defaults();
}

void Booz_xform::testfunc1() {
  std::cout << "Hello world from testfunc1" << std::endl;
}

void Booz_xform::testfunc2(int j) {
  std::cout << "Hello world from testfunc2. j=" << j << std::endl;
}
