#ifndef BOOZ_XFORM_H
#define BOOZ_XFORM_H

#include <string>
#include "vector_matrix.hpp"

namespace booz_xform {  

  const boozfloat pi = 3.141592653589793;
  const boozfloat mu0 = (4.0e-7) * pi;
  
  int driver(int, char**);
  
  class Booz_xform {
  private:
    void defaults();
    
  public:
    int verbose;
    int mboz, nboz;
    
    Booz_xform();
    void read_boozmn(std::string);
    void testfunc1();
    void testfunc2(int);
  };
  
}

#endif

