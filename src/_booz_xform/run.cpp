#include <iostream>
#ifdef OPENMP
#include <omp.h>
#endif
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
    std::cout << "             OUTBOARD (theta=0)        SURFACE         INBOARD (theta=pi)" << std::endl;
    std::cout << "------------------------------------------------------------------------------" << std::endl;
    std::cout << "zeta    |B|input  |B|Boozer    Error            |B|input  |B|Boozer    Error"
	      << std::endl << std::endl;
  }

  #pragma omp parallel for
  for (int j = 0; j < compute_surfs.size(); j++) {
#ifdef OPENMP
    std::cout << "Thread " << omp_get_thread_num() << " is handling index " << j << std::endl;
#endif
    surface_solve(j);
  }
}
