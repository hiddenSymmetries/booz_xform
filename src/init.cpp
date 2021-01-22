#include <iostream>
#include <iomanip>
#include <cassert>
#include "booz_xform.hpp"
#include "init_trig.hpp"

using namespace booz_xform;

/** Initialize variables that are common to all surfaces that share
 *  the same (mpol, ntor) resolution.
 *
 *  In the original fortran code, this work was done in read_wout_booz.f,
 *  setup_booz.f, boozer_coords.f, and foranl.f.
 */
void Booz_xform::init() {
  // First comes some initialization from read_wout_booz.f:
  
  // Set resolution for the real-space grids:
  nu = 2 * (2 * mboz + 1);
  nv = 2 * (2 * nboz + 1);
  if (nboz == 0) nv = 1;
  nu2_b = nu / 2 + 1; // Note integer division, although nu is always even so there is never rounding.
  
  mnboz = (2 * nboz + 1) * (mboz - 1) + nboz + 1;
  xmb.setZero(mnboz);
  xnb.setZero(mnboz);
  if (verbose > 0) {
    std::cout << "Initializing with mboz=" << mboz << ", nboz=" << nboz << std::endl;
    std::cout << "nu = " << nu << ", nv = " << nv << std::endl;
  }
  
  // Done with the bits from read_wout_booz.f.
  // Now comes steps from setup_booz.f.
  
  int nmin, j = 0;
  for (int m = 0; m <= mboz - 1; m++) {
    nmin = -nboz;
    if (m == 0) nmin = 0;
    for (int n = nmin; n <= nboz; n++) {
      xmb[j] = m;
      xnb[j] = n * nfp;
      j++;
    }
  }
  assert (j == mnboz);

  // A section about "fac" and "scl" has been moved to the end of
  // surface_solve(). Look for "fourier_factor" there.

  /*
  ntorsum[0] = 0.0;
  ntorsum[1] = 1.0;
  for (j = 0; j < mnmax; j++) {
    if (xm[j] == 0) ntorsum[0]++;
    if (xm[j] <= 1) ntorsum[1]++;
  }
  */

  hs = 1.0 / (ns - 1.0);
  /*
  sfull.setZero(ns);
  for (j = 2; j <= ns; j++) {
    sfull[j-1] = sqrt(hs * (j - 1));
  }
  */

  // Done with the steps from setup.booz.f.
  // Now comes some work from boozer_coords.f.

  if (asym) {
    nu3_b = nu;
  } else {
    nu3_b = nu2_b;
  }
  n_theta_zeta = nu3_b * nv;

  // Done with the steps from boozer_coords.f.
  // Now come the steps from foranl.f.
  // Note that in foranl, "nu" means nu3_b.
  
  if (asym) {
    d_theta = twopi / nu3_b;
    // But nu3_b = nu in this case, so d_theta = twopi / nu;
  } else {
    d_theta = twopi / (2 * (nu3_b - 1));
    // But nu3_b = nu2_b = nu / 2 + 1 in this case, so
    // d_theta = twopi / (2 * (nu / 2)) = twopi / nu.
  }

  d_zeta = twopi / (nfp * nv);

  // Initialize grids of (theta, zeta) for integration:
  theta_grid.setZero(nu3_b * nv);
  zeta_grid.setZero(nu3_b * nv);
  int index = 0;
  for (int j_theta = 0; j_theta < nu3_b; j_theta++) {
    for (int j_zeta = 0; j_zeta < nv; j_zeta++) {
      theta_grid[index] = j_theta * d_theta;
      zeta_grid[index] = j_zeta * d_zeta;
      index++;
    }
  }

  cosm.setZero(n_theta_zeta, mpol);
  sinm.setZero(n_theta_zeta, mpol);
  cosn.setZero(n_theta_zeta, ntor + 1);
  sinn.setZero(n_theta_zeta, ntor + 1);
  cosm_nyq.setZero(n_theta_zeta, mpol_nyq + 1);
  sinm_nyq.setZero(n_theta_zeta, mpol_nyq + 1);
  cosn_nyq.setZero(n_theta_zeta, ntor_nyq + 1);
  sinn_nyq.setZero(n_theta_zeta, ntor_nyq + 1);
  cosm_b.setZero(n_theta_zeta, mboz + 1);
  sinm_b.setZero(n_theta_zeta, mboz + 1);
  cosn_b.setZero(n_theta_zeta, nboz + 1);
  sinn_b.setZero(n_theta_zeta, nboz + 1);
  
  init_trig(theta_grid, zeta_grid,
	    cosm, sinm, cosn, sinn,
	    mpol - 1, ntor, nfp);
  
  init_trig(theta_grid, zeta_grid,
	    cosm_nyq, sinm_nyq, cosn_nyq, sinn_nyq,
	    mpol_nyq, ntor_nyq, nfp);

  /*
  std::cout << std::setprecision(15) << "cosm:" << std::endl << cosm << std::endl;
  std::cout << std::setprecision(15) << "cosn:" << std::endl << cosn << std::endl;
  std::cout << std::setprecision(15) << "sinm:" << std::endl << sinm << std::endl;
  std::cout << std::setprecision(15) << "sinn:" << std::endl << sinn << std::endl;
  */
  ns_b = jlist.size();
  Boozer_I.setZero(ns_b);
  Boozer_G.setZero(ns_b);
  
  wmns.setZero(mnmax_nyq, ns_b); // Note mnmax_nyq instead of mnboz for this one.
  bmnc_b.setZero(mnboz, ns_b);
  rmnc_b.setZero(mnboz, ns_b);
  zmns_b.setZero(mnboz, ns_b);
  pmns_b.setZero(mnboz, ns_b);
  gmnc_b.setZero(mnboz, ns_b);
  if (asym) {
    wmnc.setZero(mnmax_nyq, ns_b); // Note mnmax_nyq instead of mnboz for this one.
    bmns_b.setZero(mnboz, ns_b);
    rmns_b.setZero(mnboz, ns_b);
    zmnc_b.setZero(mnboz, ns_b);
    pmnc_b.setZero(mnboz, ns_b);
    gmns_b.setZero(mnboz, ns_b);
  }
  
  r.setZero(n_theta_zeta);
  z.setZero(n_theta_zeta);
  lambda.setZero(n_theta_zeta);
  d_lambda_d_theta.setZero(n_theta_zeta);
  d_lambda_d_zeta.setZero(n_theta_zeta);
  w.setZero(n_theta_zeta);
  d_w_d_theta.setZero(n_theta_zeta);
  d_w_d_zeta.setZero(n_theta_zeta);
  bmod.setZero(n_theta_zeta);
  p.setZero(n_theta_zeta);
  d_p_d_theta.setZero(n_theta_zeta);
  d_p_d_zeta.setZero(n_theta_zeta);
  zeta_Boozer_grid.setZero(n_theta_zeta);
  theta_Boozer_grid.setZero(n_theta_zeta);
  d_Boozer_d_vmec.setZero(n_theta_zeta);
  boozer_jacobian.setZero(n_theta_zeta);
}
