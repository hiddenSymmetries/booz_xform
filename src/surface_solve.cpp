#include "booz_xform.hpp"

using namespace booz_xform;

/** Execute the transformation to Boozer coordinates for a single flux surface.
 *
 *  This subroutine corresponds to boozer_coords.f in the fortran version.
 */
void Booz_xform::surface_solve(int js_b) {
  int js = jlist[js_b];
  // js is the radial index into the original VMEC arrays
  // js_b is the index among the subset of surfaces on which the transformation is done.

  int jmn, m, n, abs_n, j;
  
  // This next bit corresponds to transpmn.f in the fortran version.
  
  // Compute the part of p that is independent of lambda (the
  // covariant source terms in eq (10)). The net p in real space is
  //
  // p(FINAL) = { SUM(m,n)[pmns(m,n)*SIN(arg) + pmnc(m,n)*COS(arg)] - I * lambda } / D
  //
  // WHERE arg = m * theta - n * zeta and D = G + iota * I
  for (jmn = 0; jmn < mnmax_nyq; jmn++) {
    m = xm_nyq[jmn];
    n = xn_nyq[jmn];
    
    if (m != 0) {
      pmns(jmn, js_b) = bsubumnc(jmn, js) / m;
      if (asym) pmnc(jmn, js_b) = -bsubumns(jmn, js) / m;
    } else if (n != 0) {
      pmns(jmn, js_b) = -bsubvmnc(jmn, js) / n;
      if (asym) pmnc(jmn, js_b) = bsubvmns(jmn, js) / n;
    } else {
      // So m = n = 0.
      Boozer_I[js_b] = bsubvmnc(jmn, js);
      Boozer_G[js_b] = bsubumnc(jmn, js);
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
  r = 0;
  z = 0;
  lambda = 0;
  d_lambda_d_theta = 0;
  d_lambda_d_zeta = 0;
  for (jmn = 0; jmn < mnmax; jmn++) {
    m = xm[jmn];
    n = xn[jmn] / nfp;
    abs_n = std::abs(n);
    if (n < 0) {
      sign = -1;
    } else {
      sign = 1;
    }

    for (j = 0; j < n_theta_zeta; j++) {
      tcos = cosm(j, m) * cosn(j, abs_n) + sinm(j, m) * sinn(j, abs_n) * sign;
      tsin = sinm(j, m) * cosn(j, abs_n) - cosm(j, m) * sinn(j, abs_n) * sign;
      r[j]                +=  tcos * rmnc(jmn, js);
      z[j]                +=  tsin * zmns(jmn, js);
      lambda[j]           +=  tsin * lmns(jmn, js);
      d_lambda_d_theta[j] +=  tcos * lmns(jmn, js) * m;
      d_lambda_d_zeta[j]  += -tcos * lmns(jmn, js) * n;
      if (!asym) continue;
      r[j]                +=  tsin * rmns(jmn, js);
      z[j]                +=  tcos * zmnc(jmn, js);
      lambda[j]           +=  tcos * lmnc(jmn, js);
      d_lambda_d_theta[j] += -tsin * lmnc(jmn, js) * m;
      d_lambda_d_zeta[j]  +=  tsin * lmnc(jmn, js) * n;
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
  w = 0;
  d_w_d_theta = 0;
  d_w_d_zeta = 0;
  for (jmn = 0; jmn < mnmax_nyq; jmn++) {
    m = xm_nyq[jmn];
    n = xn_nyq[jmn] / nfp;
    abs_n = std::abs(n);
    if (n < 0) {
      sign = -1;
    } else {
      sign = 1;
    }

    for (j = 0; j < n_theta_zeta; j++) {
      tcos = cosm_nyq(j, m) * cosn_nyq(j, abs_n) + sinm_nyq(j, m) * sinn_nyq(j, abs_n) * sign;
      tsin = sinm_nyq(j, m) * cosn_nyq(j, abs_n) - cosm_nyq(j, m) * sinn_nyq(j, abs_n) * sign;
      w[j]           +=  tsin * pmns(jmn, js_b);
      d_w_d_theta[j] +=  tcos * pmns(jmn, js_b) * m;
      d_w_d_zeta[j]  += -tcos * pmns(jmn, js_b) * n;
      bmod[j]        +=  tcos * bmnc(jmn, js);
      if (!asym) continue;
      w[j]           +=  tcos * pmnc(jmn, js_b);
      d_w_d_theta[j] += -tsin * pmnc(jmn, js_b) * m;
      d_w_d_zeta[j]  +=  tsin * pmnc(jmn, js_b) * n;
      bmod[j]        +=  tsin * bmns(jmn, js);
    }
  }
  // End of the part corresponding to vcooords_w() in fortran.
}
