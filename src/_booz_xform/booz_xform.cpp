#include <iostream>
#include "booz_xform.hpp"

using namespace booz_xform;

void Booz_xform::defaults() {
  completed = false;
  verbose = 1;
  mboz = 16;
  nboz = 16;
  aspect = 0;
}

Booz_xform::Booz_xform() {
  // if (verbose > 0) std::cout << "Booz_xform object created." << std::endl;
  defaults();
}
