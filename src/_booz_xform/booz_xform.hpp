#ifndef BOOZ_XFORM_H
#define BOOZ_XFORM_H

#include <string>
#include <valarray>
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
    bool completed;
    Vector theta_grid; //!< Theta values of the (theta, zeta) grid, reshaped from 2D -> 1D.
    Vector zeta_grid; //!< Theta values of the (theta, zeta) grid, reshaped from 2D -> 1D.
    Matrix cosm; //!< Stores cos(m*theta) for xm
    Matrix cosn; //!< Stores cos(n*zeta) for xn
    Matrix sinm; //!< Stores sin(m*theta) for xm
    Matrix sinn; //!< Stores sin(n*zeta) for xn
    Matrix cosm_nyq; //!< Stores cos(m*theta) for xm_nyq
    Matrix cosn_nyq; //!< Stores cos(n*zeta) for xn_nyq
    Matrix sinm_nyq; //!< Stores sin(m*theta) for xm_nyq
    Matrix sinn_nyq; //!< Stores sin(n*zeta) for xn_nyq
    
    void defaults();
    void check_accuracy(int, int, Vector&, Vector&, Vector&, Matrix&, Matrix&, Matrix&, Matrix&);
    
  public:
    int verbose;
    int mboz, nboz;
    bool asym; //!< false if the configuration is stellarator-symmetric, true otherwise.

    //! The number of flux surfaces on which input data is stored.
    /**
     * The transformation to Boozer coordinates is not necessarily run on all
     * of these surfaces, only the ones indicated by compute_surfs.
     */
    int ns_in;
    
    int mpol, ntor, mnmax, mpol_nyq, ntor_nyq, mnmax_nyq, nfp, mnboz;
    Matrix rmnc, rmns, zmnc, zmns, lmnc, lmns, bmnc, bmns;
    Matrix bsubumnc, bsubumns, bsubvmnc, bsubvmns;
    Vector iota;
    IntVector xm, xn, xm_nyq, xn_nyq;

    /** List of the 0-based indices of the surfaces on which to perform the transformation.
     */
    IntVector compute_surfs;

    /** Values of normalized toroidal flux for which the input data is
	stored. These numbers are used only for plotting.
     */
    Vector s_in;
    
    /** Values of normalized toroidal flux for which the output data is
	stored. These numbers are used only for plotting.
     */
    Vector s_b;
    
    IntVector xm_b, xn_b;
    int ns_b; //!< Number of surfaces on which the transformation is calculated
    Matrix bmnc_b, rmnc_b, zmns_b, pmns_b, gmnc_b;
    Matrix bmns_b, rmns_b, zmnc_b, pmnc_b, gmns_b;
    Vector Boozer_G; //!< The toroidal covariant component of B in Boozer coordinates.
    Vector Boozer_I; //!< The poloidal covariant component of B in Boozer coordinates.

    //! Constructor
    /**
     * Create a Booz_xform object with no data.
     */
    Booz_xform();

    //! Read input data from a VMEC wout_*.nc file.
    /**
     * This method also calls init_from_vmec().
     *
     * @param[in] filename The name of the VMEC wout file to load
     */
    void read_wout(std::string filename);

    //! Handle radial dimension of vmec arrays.
    /**
     * This method also handles radial interpolation of the full-grid
     * quantities rmnc, rmns, zmnc, and zmns onto the half-grid
     * points.  It also discards the first radial index for half-grid
     * quantities, which is populated with zeros in vmec.
     */
    void init_from_vmec(int ns,
			Vector& iotas,
			Matrix& rmnc,
			Matrix& rmns,
			Matrix& zmnc,
			Matrix& zmns,
			Matrix& lmnc,
			Matrix& lmns,
			Matrix& bmnc,
			Matrix& bmns,
			Matrix& bsubumnc,
			Matrix& bsubumns,
			Matrix& bsubvmnc,
			Matrix& bsubvmns);
    
    //! Carry out the transformation calculation
    /**
     * This method is the main computationally intensive step, and it
     * should be called after the equilibrium data is loaded.
     */
    void run();

    //! Write results to a boozmn_*.nc output file
    /**
     * This method writes a classic output file, and it should be called
     * after run() has completed.
     *
     * @param[in] filename The full name of the boozmn_*.nc file to write.
     */
    void write_boozmn(std::string filename);
    
    void init();
    void surface_solve(int);
    void read_boozmn(std::string);
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
   lam -> lambda
   lt -> d_lambda_d_theta
   lz -> d_lambda_d_zeta
   xjac -> d_Boozer_d_vmec
 */
#endif

