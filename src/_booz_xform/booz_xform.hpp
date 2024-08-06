#ifndef BOOZ_XFORM_H
#define BOOZ_XFORM_H

#include <string>
#include <valarray>
#include "vector_matrix.hpp"

namespace booz_xform {

  const boozfloat pi = 3.141592653589793;
  const boozfloat twopi = 2.0 * pi;
  const boozfloat mu0 = (4.0e-7) * pi;

  static Vector defaultInit = Vector::Zero(0);
  static Vector& defaultInitPtr = defaultInit;

  int driver(int, char**);

  // Trick for passing version number from setup.py to C++, from
  // https://github.com/pybind/cmake_example/blob/master/src/main.cpp
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
#ifdef VERSION_INFO
  const std::string version = MACRO_STRINGIFY(VERSION_INFO);
#else
  const std::string version = "development version";
#endif

  class Booz_xform {
  private:
    int nu2_b; //!< Not sure
    int nu3_b; //!< Not sure
    int ntheta; //!< Number of real-space grid points in the original poloidal angle
    int nzeta; //!< Number of real-space grid points in the original toroidal angle
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

    /* Public variables are labeled below as being an input or
       output. The input quantities should be set before calling
       run(). (Many of these input quantities can be set by calling
       read_wout()).  The output quantities are populated when run()
       is called.
    */

    /** (input) Set this to 0 for no output to stdout, 1 for some
	output, 2 for lots of output.
    */
    int verbose;

    /** (input) false if the configuration is stellarator-symmetric,
	true otherwise.
     */
    bool asym;

    /** (input) Number of field periods, i.e. the discrete toroidal
	rotation symmetry.
     */
    int nfp;

    /** (input) Maximum poloidal mode number for the input arrays
	rmnc, rmns, zmnc, zmns, lmnc, and lmns.
    */
    int mpol;

    /** (input) Maximum toroidal mode number (divided by nfp) for the
	input arrays rmnc, rmns, zmnc, zmns, lmnc, and lmns.
    */
    int ntor;

    /** (input) mnmax Number of Fourier modes for the input arrays
	rmnc, rmns, zmnc, zmns, lmnc, and lmns.
     */
    int mnmax;

    /** (input) Maximum poloidal mode number for the input arrays
	bmnc, bmns, bsubumnc, bsubumns, bsubvmnc, and bsubvmns.
     */
    int mpol_nyq;

    /** (input) Maximum toroidal mode number (divided by nfp) for the
	input arrays bmnc, bmns, bsubumnc, bsubumns, bsubvmnc, and
	bsubvmns.
     */
    int ntor_nyq;

    /** (input) Total number of Fourier modes for the input arrays
	bmnc, bmns, bsubumnc, bsubumns, bsubvmnc, and bsubvmns.
     */
    int mnmax_nyq;

    /** (size mnmax, input) The poloidal Fourier mode numbers for the
	input arrays rmnc, rmns, zmnc, zmns, lmnc, and lmns.
    */
    IntVector xm;

    /** (size mnmax, input) The toroidal Fourier mode numbers for the
	input arrays rmnc, rmns, zmnc, zmns, lmnc, and lmns. The
	values should all be integer multiples of nfp.
    */
    IntVector xn;

    /** (size mnmax_nyq, input) The poloidal Fourier mode numbers for
	the input arrays bmnc, bmns, bsubumnc, bsubumns, bsubvmnc, and
	bsubvmns.
    */
    IntVector xm_nyq;

    /** (size mnmax_nyq, input) The toroidal Fourier mode numbers for
	the input arrays bmnc, bmns, bsubumnc, bsubumns, bsubvmnc, and
	bsubvmns. The values should all be integer multiples of nfp.
    */
    IntVector xn_nyq;

    /** (input) The number of flux surfaces on which input data is
	supplied.  The transformation to Boozer coordinates is not
	necessarily run on all of these surfaces, only the ones
	indicated by compute_surfs.
     */
    int ns_in;

    /** (size ns_in, input) Values of normalized toroidal flux for
	which the input data is stored. These numbers are used only
	for plotting.
     */
    Vector s_in;

    /** (size ns_in, input) Rotational transform on the input radial
	surfaces.
    */
    Vector iota;

    /** (size mnmax x ns_in, input) cos(m * theta_0 - n * zeta_0)
	Fourier modes of the major radius R.
     */
    Matrix rmnc;

    /** (size mnmax x ns_in, input) sin(m * theta_0 - n * zeta_0)
	Fourier modes of the major radius R. For stellarator-symmetric
	configurations, this array is not used and need not be
	specified.
     */
    Matrix rmns;

    /** (size mnmax x ns_in, input) cos(m * theta_0 - n * zeta_0)
	Fourier modes of the Cartesian coordinate Z of the flux
	surfaces. For stellarator-symmetric configurations, this array
	is not used and need not be specified.
     */
    Matrix zmnc;

    /** (size mnmax x ns_in, input) sin(m * theta_0 - n * zeta_0)
	Fourier modes of the Cartesian coordinate Z of the flux
	surfaces.
     */
    Matrix zmns;

    /** (size mnmax x ns_in, input) cos(m * theta_0 - n zeta_0)
	Fourier modes of lambda = theta^* - theta_0, the difference
	between the original poloidal angle theta_0 and the straight
	field line poloidal angle theta^*. For stellarator-symmetric
	configurations, this array is not used and need not be
	specified.
    */
    Matrix lmnc;

    /** (size mnmax x ns_in, input) sin(m * theta_0 - n zeta_0)
	Fourier modes of lambda = theta^* - theta_0, the difference
	between the original poloidal angle theta_0 and the straight
	field line poloidal angle theta^*.
    */
    Matrix lmns;

    /** (size mnmax_nyq x ns_in, input) cos(m * theta_0 - n * zeta_0)
	Fourier modes of the magnetic field strength B.
     */
    Matrix bmnc;

    /** (size mnmax_nyq x ns_in, input) sin(m * theta_0 - n * zeta_0)
	Fourier modes of the magnetic field strength B. For
	stellarator-symmetric configurations, this array is not used
	and need not be specified.
     */
    Matrix bmns;

    /** (size mnmax_nyq x ns_in, input) cos(m * theta_0 - n * zeta_0)
	Fourier modes of B dot (d r / d theta_0) where r is the
	position vector.
    */
    Matrix bsubumnc;

    /** (size mnmax_nyq x ns_in, input) sin(m * theta_0 - n * zeta_0)
	Fourier modes of B dot (d r / d theta_0) where r is the
	position vector.  For stellarator-symmetric configurations,
	this array is not used and need not be specified.
    */
    Matrix bsubumns;

    /** (size mnmax_nyq x ns_in, input) cos(m * theta_0 - n * zeta_0)
	Fourier modes of B dot (d r / d zeta_0) where r is the
	position vector.
    */
    Matrix bsubvmnc;

    /** (size mnmax_nyq x ns_in, input) sin(m * theta_0 - n * zeta_0)
	Fourier modes of B dot (d r / d zeta_0) where r is the
	position vector.  For stellarator-symmetric configurations,
	this array is not used and need not be specified.
    */
    Matrix bsubvmns;

    /** (input) Maximum poloidal mode number for representing output
	quantities in Boozer coordinates.
    */
    int mboz;

    /** (input) Maximum toroidal mode number (divided by nfp) for
	representing output quantities in Boozer coordinates. For
	example, if nboz=2 and nfp=10, the toroidal modes used will be
	n=-20, -10, 0, 10, 20.
     */
    int nboz;

    /** (input) Indices of ns_in-sized radial arrays, specifying the
	flux surfaces for which the transformation to Boozer
	coordinates will be performed. All values should be >= 0 and <
	ns_in.  The array compute_surfs is similar to the array jlist
	in the earlier fortran booz_xform program, with compute_surfs
	= jlist - 2.
     */
    IntVector compute_surfs;

    /** (input) The aspect ratio of the configuration. This value is
	not used for anything by booz_xform, and does not need to be
	set. It is provided as a means to pass this value from the
	input equilibrium to booz_xform output files.
     */
    boozfloat aspect;

    /** (input) The boundary toroidal flux of the configuration (not
	divided by (2*pi)). This value is not used for anything by
	booz_xform, and does not need to be set. It is provided as a
	means to pass this value from the input equilibrium to
	booz_xform output files.
     */
    boozfloat toroidal_flux;

    // End of the inputs. Now come the outputs.

    /** (output) Number of surfaces on which the transformation is
	calculated.
     */
    int ns_b;

    /** (size ns_b, output) Values of normalized toroidal flux for
	which the output data is stored. These numbers are used only
	for plotting.
     */
    Vector s_b;

    /** (output) Total number of Fourier modes for output data.
     */
    int mnboz;

    /** (size mnboz, output) Poloidal mode numbers used for the
	Fourier representation of output quantities, i.e. functions of
	the Boozer angles.
    */
    IntVector xm_b;

    /** (size mnboz, output) Toroidal mode numbers used for the
	Fourier representation of output quantities, i.e. functions of
	the Boozer angles.
    */
    IntVector xn_b;

    /** (size mnboz x ns_b, output) cos(m * theta_B - n * zeta_B)
	Fourier modes of the magnetic field strength in Boozer
	coordinates.
    */
    Matrix bmnc_b;

    /** (size mnboz x ns_b, output) sin(m * theta_B - n * zeta_B)
	Fourier modes of the magnetic field strength in Boozer
	coordinates. If the configuration is stellarator-symmetric,
	this quantity is zero so the array will have size 0 x 0.
    */
    Matrix bmns_b;

    /** (size mnboz x ns_b, output) cos(m * theta_B - n * zeta_B)
	Fourier modes (with respect to Boozer coordinates) of the
	Jacobian of (psi, theta_B, zeta_B) coordinates.
    */
    Matrix gmnc_b;

    /** (size mnboz x ns_b, output) sin(m * theta_B - n * zeta_B)
	Fourier modes (with respect to Boozer coordinates) of the
	Jacobian of (psi, theta_B, zeta_B) coordinates. If the
	configuration is stellarator-symmetric, this quantity is zero
	so this array will have size 0 x 0.
    */
    Matrix gmns_b;

    /** (size mnboz x ns_b, output) cos(m * theta_B - n * zeta_B)
	Fourier modes (with respect to Boozer coordinates) of the
	major radius R of the flux surfaces.
    */
    Matrix rmnc_b;

    /** (size mnboz x ns_b, output) sin(m * theta_B - n * zeta_B)
	Fourier modes (with respect to Boozer coordinates) of the
	major radius R of the flux surfaces. If the configuration is
	stellarator-symmetric, this quantity is zero so this array
	will have size 0 x 0.
    */
    Matrix rmns_b;

    /** (size mnboz x ns_b, output) cos(m * theta_B - n * zeta_B)
	Fourier modes (with respect to Boozer coordinates) of the
	Cartesian coordinate Z of the flux surfaces. If the
	configuration is stellarator-symmetric, this quantity is zero
	so array will have size 0 x 0.
    */
    Matrix zmnc_b;

    /** (size mnboz x ns_b, output) sin(m * theta_B - n * zeta_B)
	Fourier modes (with respect to Boozer coordinates) of the
	Cartesian coordinate Z of the flux surfaces.
    */
    Matrix zmns_b;

    /** (size mnboz x ns_b, output) cos(m * theta_B - n * zeta_B)
	Fourier modes (with respect to Boozer coordinates) of the
	toroidal angle difference nu = zeta_B - zeta_0. If the
	configuration is stellarator-symmetric, this quantity is zero
	so this array will have size 0 x 0.
    */
    Matrix numnc_b;

    /** (size mnboz x ns_b, output) sin(m * theta_B - n * zeta_B)
	Fourier modes (with respect to Boozer coordinates) of the
	toroidal angle difference nu = zeta_B - zeta_0.
    */
    Matrix numns_b;

    /** (size ns_b, output) Coefficient of grad zeta_B in the
	covariant representation of the magnetic field vector B in
	Boozer coordinates, evaluated on the magnetic surfaces used
	for output quantities.
    */
    Vector Boozer_G;

    /** (size ns_in, output) Coefficient of grad zeta_B in the
	covariant representation of the magnetic field vector B in
	Boozer coordinates, evaluated on all the magnetic surfaces for
	which input data was provided.
    */
    Vector Boozer_G_all;

    /** (size ns_b, output) Coefficient of grad theta_B in the
	covariant representation of the magnetic field vector B in
	Boozer coordinates, evaluated on the magnetic surfaces used
	for output quantities.
    */
    Vector Boozer_I;

    /** (size ns_in, output) Coefficient of grad theta_B in the
	covariant representation of the magnetic field vector B in
	Boozer coordinates, evaluated on all the magnetic surfaces for
	which input data was provided.
    */
    Vector Boozer_I_all;

    /** (size ns_in + 1, output) Uniformly spaced grid going from 0 to the boundary
    toroidal flux (not divided by (2*pi)), evaluated on full vmec grid.
    */
    Vector phi;

    /** (size ns_in + 1, output) The derivative of the toroidal flux (not divided by (2*pi))
    with respect to s, evaluated on full vmec grid.
    */
    Vector phip;

    /** (size ns_in + 1, output) Uniformly spaced grid going from 0 to the boundary
    poloidal flux (not divided by (2*pi)), evaluated on full vmec grid.
    */
    Vector chi;

    /** Pressure evaluated on the full vmec grid **/
    Vector pres;

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
    void read_wout(std::string filename, bool flux=false);

    //! Handle radial dimension of vmec arrays.
    /**
     * This method also handles radial interpolation of the full-grid
     * quantities rmnc, rmns, zmnc, and zmns onto the half-grid
     * points.  It also discards the first radial index for half-grid
     * quantities, which is populated with zeros in vmec.
     *
     * The input parameters to this function are the variables of the
     * same name in a VMEC wout file.
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
			Matrix& bsubvmns,
      Vector& phips=defaultInitPtr,
      Vector& chi=defaultInitPtr,
      Vector& pres=defaultInitPtr,
      Vector& phi=defaultInitPtr);

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

    //! Read previously calculated results from a boozmn_*.nc output file
    /**
     * @param[in] filename The full name of the boozmn_*.nc file to read.
     */
    void read_boozmn(std::string filename);

    void init();

    void surface_solve(int js_b);

  };

}

/* Translating between fortran and C++ variable names:

   nu_b, nv_b -> ntheta, nzeta
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
