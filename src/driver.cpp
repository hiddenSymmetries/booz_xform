#include <iostream>
#include "booz_xform.hpp"

//using namespace booz_xform;

int booz_xform::driver(int argc, char* argv[]) {
  std::cout << "This is xbooz_xform." << std::endl;

  booz_xform::Vector v;

  v.resize(5, 3.14);
  std::cout << v << std::endl;

  booz_xform::Booz_xform booz;
  //booz.read_boozmn("/Users/mattland/pyBooz_xform/boozmn_n3are_R7.75B5.7.nc");

  // booz.read_wout("../tests/test_files/wout_li383_1.4m.nc");
    
  std::cout << "Good bye." << std::endl;

  return 0;
}
