#include <iostream>
#include <cassert>
#include "booz_xform.hpp"

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
  xmb.resize(mnboz, 0.0);
  xnb.resize(mnboz, 0.0);

  // Done with the bits from read_wout_booz.f.
  // Now comes steps from setup_booz.f.
  
  int nmin, j = 0;
  for (int m = 0; m <= mboz - 1; m++) {
    nmin = -nboz;
    if (m == 0) nmin = 0;
    for (int n = nmin; n <= nboz; n++) {
      xmb[j] = m;
      xnb[j] = n;
      j++;
    }
  }
  assert (j == mnboz);

  boozfloat fac;
  if (asym) {
    fac = 2.0 / (nu * nv);
  } else {
    fac = 2.0 / ((nu2_b - 1) * nv);
    // Equivalently, fac = 4.0 / (nu * nv)
  }

  fourier_factor.resize(mnboz, fac);
  fourier_factor[0] = fac / 2.0;

  ntorsum[0] = 0.0;
  ntorsum[1] = 1.0;
  for (j = 0; j < mnmax; j++) {
    if (xm[j] == 0) ntorsum[0]++;
    if (xm[j] <= 1) ntorsum[1]++;
  }

  hs = 1.0 / (ns - 1.0);
  /*
  sfull.resize(ns, 0.0);
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
  theta_grid.resize(nu3_b * nv, 0.0);
  zeta_grid.resize(nu3_b * nv, 0.0);
  int index = 0;
  for (int j_theta = 0; j_theta < nu3_b; j_theta++) {
    for (int j_zeta = 0; j_zeta < nv; j_zeta++) {
      theta_grid[index] = j_theta * d_theta;
      zeta_grid[index] = j_zeta * d_zeta;
      index++;
    }
  }

  cosm.resize(n_theta_zeta, mpol, 0.0);
  sinm.resize(n_theta_zeta, mpol, 0.0);
  cosn.resize(n_theta_zeta, ntor + 1, 0.0);
  sinn.resize(n_theta_zeta, ntor + 1, 0.0);
  cosm_nyq.resize(n_theta_zeta, mpol_nyq + 1, 0.0);
  sinm_nyq.resize(n_theta_zeta, mpol_nyq + 1, 0.0);
  cosn_nyq.resize(n_theta_zeta, ntor_nyq + 1, 0.0);
  sinn_nyq.resize(n_theta_zeta, ntor_nyq + 1, 0.0);
  
  init_trig(cosm, sinm, cosn, sinn, mpol - 1, ntor);
  init_trig(cosm_nyq, sinm_nyq, cosn_nyq, sinn_nyq, mpol_nyq, ntor_nyq);

  ns_b = jlist.size();
  Boozer_I.resize(ns_b, 0.0);
  Boozer_G.resize(ns_b, 0.0);
  
  pmns.resize(mnmax_nyq, ns_b, 0.0); // Note mnmax_nyq instead of mnboz for this one.
  bmnc_b.resize(mnboz, ns_b, 0.0);
  rmnc_b.resize(mnboz, ns_b, 0.0);
  zmns_b.resize(mnboz, ns_b, 0.0);
  pmns_b.resize(mnboz, ns_b, 0.0);
  gmnc_b.resize(mnboz, ns_b, 0.0);
  if (asym) {
    pmnc.resize(mnmax_nyq, ns_b, 0.0); // Note mnmax_nyq instead of mnboz for this one.
    bmns_b.resize(mnboz, ns_b, 0.0);
    rmns_b.resize(mnboz, ns_b, 0.0);
    zmnc_b.resize(mnboz, ns_b, 0.0);
    pmnc_b.resize(mnboz, ns_b, 0.0);
    gmns_b.resize(mnboz, ns_b, 0.0);
  }
  
  r.resize(n_theta_zeta, 0.0);
  z.resize(n_theta_zeta, 0.0);
  lambda.resize(n_theta_zeta, 0.0);
  d_lambda_d_theta.resize(n_theta_zeta, 0.0);
  d_lambda_d_zeta.resize(n_theta_zeta, 0.0);
  w.resize(n_theta_zeta, 0.0);
  d_w_d_theta.resize(n_theta_zeta, 0.0);
  d_w_d_zeta.resize(n_theta_zeta, 0.0);
  bmod.resize(n_theta_zeta, 0.0);
}

///////////////////////////////////////////////////

/** Evaluate sin(m*theta), cos(m*theta), sin(n*zeta), and cos(n*zeta) for a range
 *  of (m, n, theta, zeta) values.
 *
 *  This function is called trigfunc in the fortran version.
 *
 *  This function assumes that all of the matrices have been allocated
 *  to the appropriate size beforehand.
 */
void Booz_xform::init_trig(Matrix& cosmx, Matrix& sinmx, Matrix& cosnx, Matrix& sinnx,
			   int mmax, int nmax) {
  int j, m, n;

  // Initialize modes with (m,n) = 0 or 1:
  for (j = 0; j < n_theta_zeta; j++) {
    cosmx(j, 0) = 1.0;
    cosnx(j, 0) = 1.0;
    cosmx(j, 1) = cos(theta_grid[j]);
    sinmx(j, 1) = sin(theta_grid[j]);
  }
  if (nmax >= 1) {
    for (j = 0; j < n_theta_zeta; j++) {
      cosnx(j, 1) = cos(nfp * zeta_grid[j]);
      sinnx(j, 1) = sin(nfp * zeta_grid[j]);
    }
  }

  // Evaluate the rest of the m values:
  for (m = 2; m <= mmax; m++) {
    for (j = 0; j < n_theta_zeta; j++) {
      cosmx(j, m) = cosmx(j, m - 1) * cosmx(j, 1) - sinmx(j, m - 1) * sinmx(j, 1);
      sinmx(j, m) = sinmx(j, m - 1) * cosmx(j, 1) + cosmx(j, m - 1) * sinmx(j, 1);
    }
  }
  
  // Evaluate the rest of the n values:
  for (n = 2; n <= nmax; n++) {
    for (j = 0; j < n_theta_zeta; j++) {
      cosnx(j, n) = cosnx(j, n - 1) * cosnx(j, 1) - sinnx(j, n - 1) * sinnx(j, 1);
      sinnx(j, n) = sinnx(j, n - 1) * cosnx(j, 1) + cosnx(j, n - 1) * sinnx(j, 1);
    }
  }
}
