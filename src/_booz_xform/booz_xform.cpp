#include <iostream>
#include "booz_xform.hpp"

using namespace booz_xform;

void Booz_xform::defaults() {
  completed = false;
  ns_in = 0;
  ns_b = 0;
  mnboz = 0;
  n_theta_zeta = 0;
  verbose = 1;
  mboz = 16;
  nboz = 16;
  aspect = 0;
}

Booz_xform::Booz_xform() {
  // if (verbose > 0) std::cout << "Booz_xform object created." << std::endl;
  defaults();
}
