#ifndef BOOZ_XFORM_H
#define BOOZ_XFORM_H

#include <string>
#include <vector>
#include "vector_matrix.hpp"

namespace booz_xform {  

  const boozfloat pi = 3.141592653589793;
  const boozfloat mu0 = (4.0e-7) * pi;
  
  int driver(int, char**);
  
  class Booz_xform {
  private:
    bool completed;
    
    void defaults();
    
  public:
    int verbose;
    int mboz, nboz;
    bool asym;
    int mpol, ntor, mnmax, mnmax_nyq, ns, nfp;
    Matrix rmnc, rmns, zmnc, zmns, lmnc, lmns, bmnc, bmns;
    Matrix bsubumnc, bsubumns, bsubvmnc, bsubvmns;
    Vector iotas, xm, xn, xm_nyq, xn_nyq;
    std::vector<int> jlist;
    
    Booz_xform();
    void read_boozmn(std::string);
    void read_wout(std::string);
    void run();
    void testfunc1();
    void testfunc2(int);
  };
  
}

#endif

