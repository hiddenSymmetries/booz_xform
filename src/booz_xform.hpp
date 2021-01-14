#ifndef BOOZ_XFORM_H
#define BOOZ_XFORM_H

#include <string>
#include <vector>
#include "vector_matrix.hpp"

namespace booz_xform {  

  const boozfloat pi = 3.141592653589793;
  const boozfloat twopi = 2.0 * pi;
  const boozfloat mu0 = (4.0e-7) * pi;
  
  int driver(int, char**);
  
  class Booz_xform {
  private:
    int nu2_b; //!< Not sure
    int nu3_b; //!< Not sure
    int nu, nv; //!< Number of real-space grid points in the VMEC poloidal and toroidal angles
    boozfloat d_theta, d_zeta; //!< Spacing between grid points in the VMEC poloidal and toroidal angles
    int n_theta_zeta; //!< Number of grid points in the combined theta-zeta grid
    boozfloat hs; //!< Spacing in s between adjacent flux surfaces
    Vector fourier_factor; //!< Normalization factor used in Fourier transforms
    Vector sfull; //!< Not sure
    bool completed;
    boozfloat ntorsum[2]; //!< The number of VMEC modes with m == 0 and m <= 1.
    Vector theta_grid; //!< Theta values of the (theta, zeta) grid, reshaped from 2D -> 1D.
    Vector zeta_grid; //!< Theta values of the (theta, zeta) grid, reshaped from 2D -> 1D.
    Matrix cosm; //!< Stores cos(m*theta) for xm vs theta_grid
    Matrix cosn; //!< Stores cos(n*theta) for xn vs zeta_grid
    Matrix sinm; //!< Stores sin(m*theta) for xm vs theta_grid
    Matrix sinn; //!< Stores sin(n*theta) for xn vs zeta_grid
    Matrix cosm_nyq; //!< Stores cos(m*theta) for xm_nyq vs theta_grid
    Matrix cosn_nyq; //!< Stores cos(n*theta) for xn_nyq vs zeta_grid
    Matrix sinm_nyq; //!< Stores sin(m*theta) for xm_nyq vs theta_grid
    Matrix sinn_nyq; //!< Stores sin(n*theta) for xn_nyq vs zeta_grid
    Vector w, d_w_d_theta, d_w_d_zeta, bmod;
    
    void defaults();
    void init_trig(Matrix&, Matrix&, Matrix&, Matrix&, int, int);
    
  public:
    int verbose;
    int mboz, nboz;
    bool asym;
    int mpol, ntor, mnmax, mpol_nyq, ntor_nyq, mnmax_nyq, ns, nfp, mnboz;
    Matrix rmnc, rmns, zmnc, zmns, lmnc, lmns, bmnc, bmns;
    Matrix bsubumnc, bsubumns, bsubvmnc, bsubvmns;
    Vector iotas, xm, xn, xm_nyq, xn_nyq;
    std::vector<int> jlist;
    Vector xmb, xnb;
    int ns_b; //!< Number of surface on which the transformation is calculated
    Matrix pmns, pmnc; //!< Difference between Boozer vs Vmec toroidal angle as function of Vmec angles
    Matrix bmnc_b, rmnc_b, zmns_b, pmns_b, gmnc_b;
    Matrix bmns_b, rmns_b, zmnc_b, pmnc_b, gmns_b;
    Vector Boozer_I, Boozer_G; //!< The covariant components of B in Boozer coordinates.
    
    Booz_xform();
    void read_boozmn(std::string);
    void read_wout(std::string);
    void run();
    void init();
    void surface_solve(int);
    void testfunc1();
    void testfunc2(int);
  };
  
}

/* Translating between fortran and C++ variable names:

   nu_b, nv_b -> nu, nv
   nunv -> n_theta_zeta
   dth, dzt -> d_theta, d_zeta
   scl -> fourier_factor
   thgrd -> theta_grid
   ztgrd -> zeta_grid
   jsize -> ns_b
   wp -> w
   wt -> d_w_d_theta
   wz -> d_w_d_zeta
 */
#endif

