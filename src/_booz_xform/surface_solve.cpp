#include <iostream>
#include <fstream>
#include <iomanip>
#include "booz_xform.hpp"
#include "init_trig.hpp"

using namespace booz_xform;

/** Compute the real-space quantities on one flux surface that enter the
 *  Fourier integrals of eq (11): R, Z, lambda, nu, |B|, the Boozer angles
 *  on the VMEC grid, the Jacobian of eq (12), and the trigonometric
 *  tables in Boozer angles. This is shared by surface_solve() and by the
 *  derivative routines in derivatives.cpp.
 *
 * @param[in] js_b The index into compute_surfs indicating the surface.
 * @param[out] st The real-space quantities.
 */
void Booz_xform::surface_state(int js_b, SurfaceState& st) {
  int js = compute_surfs[js_b];
  st.js = js;
  int jmn, m, n, abs_n, j;

  st.cosm_b.setZero(n_theta_zeta, mboz + 1);
  st.sinm_b.setZero(n_theta_zeta, mboz + 1);
  st.cosn_b.setZero(n_theta_zeta, nboz + 1);
  st.sinn_b.setZero(n_theta_zeta, nboz + 1);

  st.wmns.setZero(mnmax_nyq); // Note mnmax_nyq instead of mnboz for this one.
  st.wmnc.setZero(mnmax_nyq);
  st.r.setZero(n_theta_zeta);
  st.z.setZero(n_theta_zeta);
  st.lambda.setZero(n_theta_zeta);
  st.d_lambda_d_theta.setZero(n_theta_zeta);
  st.d_lambda_d_zeta.setZero(n_theta_zeta);
  st.w.setZero(n_theta_zeta);
  st.d_w_d_theta.setZero(n_theta_zeta);
  st.d_w_d_zeta.setZero(n_theta_zeta);
  st.bmod.setZero(n_theta_zeta);

  // This next bit corresponds to transpmn.f in the fortran version.

  // Compute the part of p that is independent of lambda (the
  // covariant source terms in eq (10)). The net p in real space is
  //
  // p(FINAL) = { SUM(m,n)[wmns(m,n)*SIN(arg) + wmnc(m,n)*COS(arg)] - I * lambda } / D
  //
  // WHERE arg = m * theta - n * zeta and D = G + iota * I
  st.I = 0;
  st.G = 0;
  for (jmn = 0; jmn < mnmax_nyq; jmn++) {
    m = xm_nyq[jmn];
    n = xn_nyq[jmn];

    if (m != 0) {
      st.wmns(jmn) = bsubumnc(jmn, js) / m;
      if (asym) st.wmnc(jmn) = -bsubumns(jmn, js) / m;
    } else if (n != 0) {
      st.wmns(jmn) = -bsubvmnc(jmn, js) / n;
      if (asym) st.wmnc(jmn) = bsubvmns(jmn, js) / n;
    } else {
      // So m = n = 0.
      st.I = bsubumnc(jmn, js);
      st.G = bsubvmnc(jmn, js);
      // Without loss of generality, these modes of p can be set to 0.
    }
  }

  // End of the part corresponding to transpmn.f in the fortran version.

  // Beginning of the part corresponding to vcoords.f::vcoords_rz() in
  // the fortran version.  In this part, we transform R, Z, and lambda
  // from Fourier space to the real-space (theta, zeta) grid. At the
  // same time we also evaluate the derivatives of lambda with respect
  // to theta and zeta. R, Z, and lambda all use non-Nyquist
  // frequencies.  Lambda is already on the half grid. However, in the
  // original vmec file R and Z were on the full grid, so
  // interpolation to the half grid was needed. This interpolation was
  // done back when we loaded the VMEC file, in read_wout().

  int sign;
  boozfloat tcos, tsin;
  for (jmn = 0; jmn < mnmax; jmn++) {
    m = xm[jmn];
    n = xn[jmn];
    abs_n = std::abs(n) / nfp;
    if (n < 0) {
      sign = -1;
    } else {
      sign = 1;
    }

    for (j = 0; j < n_theta_zeta; j++) {
      tcos = cosm(j, m) * cosn(j, abs_n) + sinm(j, m) * sinn(j, abs_n) * sign;
      tsin = sinm(j, m) * cosn(j, abs_n) - cosm(j, m) * sinn(j, abs_n) * sign;
      st.r[j]                +=  tcos * rmnc(jmn, js);
      st.z[j]                +=  tsin * zmns(jmn, js);
      st.lambda[j]           +=  tsin * lmns(jmn, js);
      st.d_lambda_d_theta[j] +=  tcos * lmns(jmn, js) * m;
      st.d_lambda_d_zeta[j]  += -tcos * lmns(jmn, js) * n;
      if (!asym) continue;
      st.r[j]                +=  tsin * rmns(jmn, js);
      st.z[j]                +=  tcos * zmnc(jmn, js);
      st.lambda[j]           +=  tcos * lmnc(jmn, js);
      st.d_lambda_d_theta[j] += -tsin * lmnc(jmn, js) * m;
      st.d_lambda_d_zeta[j]  +=  tsin * lmnc(jmn, js) * n;
    }
  }
  // End of the part corresponding to vcoords_rz() in the fortran version.

  // Beginning of the part corresponding to vcoords_w().  The goal of
  // this next section is to transform several functions from Fourier
  // space to the real-space (theta, zeta) grid. One of these function
  // is "w", the right-hand side of eq (10), and its derivatives with
  // respect to theta and zeta. The other function is |B|. These
  // functions are all Nyquist-frequency quantities that are already
  // known on the half grid. This section is just like the previous
  // one, but with the difference that we use the Nyquist values for
  // (m,n).
  for (jmn = 0; jmn < mnmax_nyq; jmn++) {
    m = xm_nyq[jmn];
    n = xn_nyq[jmn];
    abs_n = std::abs(n / nfp);
    if (n < 0) {
      sign = -1;
    } else {
      sign = 1;
    }

    for (j = 0; j < n_theta_zeta; j++) {
      tcos = cosm_nyq(j, m) * cosn_nyq(j, abs_n) + sinm_nyq(j, m) * sinn_nyq(j, abs_n) * sign;
      tsin = sinm_nyq(j, m) * cosn_nyq(j, abs_n) - cosm_nyq(j, m) * sinn_nyq(j, abs_n) * sign;
      st.w[j]           +=  tsin * st.wmns(jmn);
      st.d_w_d_theta[j] +=  tcos * st.wmns(jmn) * m;
      st.d_w_d_zeta[j]  += -tcos * st.wmns(jmn) * n;
      st.bmod[j]        +=  tcos * bmnc(jmn, js);
      if (!asym) continue;
      st.w[j]           +=  tcos * st.wmnc(jmn);
      st.d_w_d_theta[j] += -tsin * st.wmnc(jmn) * m;
      st.d_w_d_zeta[j]  +=  tsin * st.wmnc(jmn) * n;
      st.bmod[j]        +=  tsin * bmns(jmn, js);
    }
  }
  // End of the part corresponding to vcooords_w() in fortran.

  // The next part is located in harfun.f in the fortran version.
  st.iota = iota[js];
  st.GI = st.G + st.iota * st.I;
  boozfloat one_over_GI = 1.0 / st.GI;

  // Get total nu from eq (10):
  st.nu = one_over_GI * (st.w - st.I * st.lambda);

  // Get Boozer poloidal and toroidal angles from eq (3):
  st.theta_Boozer_grid = theta_grid + st.lambda + st.iota * st.nu;
  st.zeta_Boozer_grid = zeta_grid + st.nu;

  // Derivatives of Eq (10) with respect to theta or zeta:
  st.d_nu_d_zeta  = one_over_GI * (st.d_w_d_zeta  - st.I * st.d_lambda_d_zeta);
  st.d_nu_d_theta = one_over_GI * (st.d_w_d_theta - st.I * st.d_lambda_d_theta);

  // Eq (12):
  st.d_Boozer_d_vmec = (1.0 + st.d_lambda_d_theta) * (1.0 + st.d_nu_d_zeta)
    + (st.iota - st.d_lambda_d_zeta) * st.d_nu_d_theta;

  // Consider removing this test in the next line for speed:
  //if (d_Boozer_d_vmec.max() * d_Boozer_d_vmec.min() <= 0)
  //  throw std::runtime_error("d_Boozer_d_vmec crosses through 0!");

  // End of the part located in harfun.f in the fortran version.

  // Now come lines from boozer.f
  // The goal now is to evaluate eq (11).

  init_trig(st.theta_Boozer_grid, st.zeta_Boozer_grid,
	    st.cosm_b, st.sinm_b, st.cosn_b, st.sinn_b, mboz, nboz, nfp);

  if (!asym) {
    // Only integrate in theta half-way around
    int i = nzeta * (nu2_b - 1) + 1;
    int imax = i - 1 + nzeta;
    for (int jmn = 0; jmn <= mboz; jmn++) {
      for (j = 0; j < nzeta; j++) {
	st.cosm_b(j, jmn) = 0.5 * st.cosm_b(j, jmn);   // theta = 0
	st.sinm_b(j, jmn) = 0.5 * st.sinm_b(j, jmn);   // Should be 0
      }
      for (j = i - 1; j < imax; j++) { // i-1 because C is 0-based.
	st.cosm_b(j, jmn) = 0.5 * st.cosm_b(j, jmn); // theta = pi
	st.sinm_b(j, jmn) = 0.5 * st.sinm_b(j, jmn); // Should be 0
      }
    }
  }

  // Jacobian of (psi, theta_Boozer, zeta_Boozer), where 2 * pi * psi is toroidal flux:
  st.boozer_jacobian = st.GI / (st.bmod * st.bmod);
}

/** Normalization of the Fourier integral of eq (11) for Boozer mode jmn.
 *
 *  This quantity is called "fac" and "scl" in the fortran version.
 */
boozfloat Booz_xform::fourier_factor(int jmn) {
  boozfloat fourier_factor0;
  if (asym) {
    fourier_factor0 = 2.0 / (ntheta * nzeta);
  } else {
    fourier_factor0 = 2.0 / ((nu2_b - 1) * nzeta);
    // Equivalently, fac = 4.0 / (ntheta * nzeta)
  }
  // There is a factor of 1/2 for the m=0 element.
  if (jmn == 0) fourier_factor0 *= 0.5;
  return fourier_factor0;
}

/** Execute the transformation to Boozer coordinates for a single flux surface.
 *
 *  This subroutine corresponds to boozer_coords.f in the fortran version.
 *
 * @param[in] js_b The index into compute_surfs indicating the surface
 * on which the Boozer transformation will be calculated.
 */
void Booz_xform::surface_solve(int js_b) {
  int js = compute_surfs[js_b];
  // js is the 0-based radial index into the input data arrays
  // js_b is the 0-based index among the subset of surfaces on which the transformation is done.
  if (verbose > 1) std::cout << "Solving for js_b=" << js_b << ", js=" << js << std::endl;

  // The work arrays in st need to be independent on each thread, so
  // threads do not over-write each other's arrays. There is some time
  // spent allocating the memory here, on each surface, rather than
  // allocating just once in init(). The time spent is small compared
  // to the loops that come later.
  SurfaceState st;
  surface_state(js_b, st);
  Boozer_I[js_b] = st.I;
  Boozer_G[js_b] = st.G;

  int jmn, m, n, abs_n, j, sign;
  boozfloat tcos, tsin, fac;
  for (jmn = 0; jmn < mnboz; jmn++) {
    m = xm_b[jmn];
    n = xn_b[jmn];
    abs_n = std::abs(n / nfp);
    if (n < 0) {
      sign = -1;
    } else {
      sign = 1;
    }

    fac = fourier_factor(jmn);

    for (j = 0; j < n_theta_zeta; j++) {
      tcos = (st.cosm_b(j, m) * st.cosn_b(j, abs_n)
	      + st.sinm_b(j, m) * st.sinn_b(j, abs_n) * sign)
	* st.d_Boozer_d_vmec[j] * fac;

      tsin = (st.sinm_b(j, m) * st.cosn_b(j, abs_n)
	      - st.cosm_b(j, m) * st.sinn_b(j, abs_n) * sign)
	* st.d_Boozer_d_vmec[j] * fac;

      bmnc_b(jmn, js_b) +=  tcos * st.bmod[j];
      rmnc_b(jmn, js_b) +=  tcos * st.r[j];
      zmns_b(jmn, js_b) +=  tsin * st.z[j];
      numns_b(jmn, js_b) += tsin * st.nu[j];
      gmnc_b(jmn, js_b) +=  tcos * st.boozer_jacobian[j];
      if (!asym) continue;
      bmns_b(jmn, js_b) +=  tsin * st.bmod[j];
      rmns_b(jmn, js_b) +=  tsin * st.r[j];
      zmnc_b(jmn, js_b) +=  tcos * st.z[j];
      numnc_b(jmn, js_b) += tcos * st.nu[j];
      gmns_b(jmn, js_b) +=  tsin * st.boozer_jacobian[j];
    }
  }

  // End of the part from boozer.f

  if (verbose > 0) check_accuracy(js, js_b, st.bmod,
				  st.theta_Boozer_grid, st.zeta_Boozer_grid,
				  st.cosm_b, st.cosn_b, st.sinm_b, st.sinn_b);
}
