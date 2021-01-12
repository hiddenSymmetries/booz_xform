#include <iostream>
#include "booz_xform.hpp"

using namespace booz_xform;

void Booz_xform::run() {
  std::cout << "jlist: ";
  for (int j = 0; j < jlist.size(); j++) std::cout << " " << jlist[j];
  std::cout << std::endl;
}
