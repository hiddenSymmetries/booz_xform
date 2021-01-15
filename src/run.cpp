#include <iostream>
#include "booz_xform.hpp"

using namespace booz_xform;

void Booz_xform::run() {
  std::cout << "jlist: ";
  for (int j = 0; j < jlist.size(); j++) std::cout << " " << jlist[j];
  std::cout << std::endl;

  init();

  if (verbose > 0) {
    std::cout << "             OUTBOARD (u=0)              JS          INBOARD (u=pi)" << std::endl;
    std::cout << "-----------------------------------------------------------------------------" << std::endl;
    std::cout << "  v     |B|vmec    |B|booz    Error             |B|vmec    |B|booz    Error"
	      << std::endl << std::endl;
  }
  
  for (int j = 0; j < jlist.size(); j++)
    surface_solve(j);
}
