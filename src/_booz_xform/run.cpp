#include <iostream>
#include "booz_xform.hpp"

using namespace booz_xform;

void Booz_xform::run() {
  if (verbose > 0) {
    std::cout << "compute_surfs (0-based indices): ";
    for (int j = 0; j < compute_surfs.size(); j++) std::cout << " " << compute_surfs[j];
    std::cout << std::endl;
  }

  init();

  if (verbose > 0) {
    std::cout << "                   |        outboard (theta=0)      |      inboard (theta=pi)      |" << std::endl;
    std::cout << "thread js_b js zeta| |B|input  |B|Boozer    Error   | |B|input  |B|Boozer    Error |" << std::endl;
    std::cout << "------------------------------------------------------------------------------------" << std::endl;
  }

  // In the next line would could try to add something like "schedule(static, 1)" for more order.
#pragma omp parallel for
  for (int j = 0; j < compute_surfs.size(); j++) {
    surface_solve(j);
  }
}
