#include <iostream>
#include <fstream>
#include <iomanip>
#include "booz_xform.hpp"
#include "init_trig.hpp"

using namespace booz_xform;

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

  int jmn, m, n, abs_n, j;

  // Next come a bunch of work arrays that need to be independent on
  // each thread, so threads do not over-write each other's
  // arrays. There is some time spent allocating the memory here, on
  // each surface, rather than allocating just once in init(). However
  // then it is trickier to have these arrays independent on each
  // thread. The time spent is small compared to the loops that come
  // later. So for simplicity, the work arrays are just allocated
  // here.
  Matrix cosm_b; //!< Stores cos(m*theta_Boozer) for xm_b
  Matrix cosn_b; //!< Stores cos(n*zeta_Boozer) for xn_b
  Matrix sinm_b; //!< Stores sin(m*theta_Boozer) for xm_b
  Matrix sinn_b; //!< Stores sin(n*zeta_Boozer) for xn_b
  Vector r;
  Vector z;
  Vector lambda;
  Vector d_lambda_d_theta;
  Vector d_lambda_d_zeta;
  Vector w;
  Vector d_w_d_theta;
  Vector d_w_d_zeta;
  Vector bmod;
  Vector theta_diff;
  Vector nu;
  Vector d_nu_d_theta;
  Vector d_nu_d_zeta;
  Vector theta_Boozer_grid;
  Vector zeta_Boozer_grid;
  Vector wmns;
  Vector wmnc; //!< Right-hand side of eq (10) in Fourier space
  Vector d_Boozer_d_vmec; //!< The Jacobian in eq (12)
  Vector boozer_jacobian; //!< (G + iota * I) / B^2 on the (theta, zeta) grid.

  cosm_b.setZero(n_theta_zeta, mboz + 1);
  sinm_b.setZero(n_theta_zeta, mboz + 1);
  cosn_b.setZero(n_theta_zeta, nboz + 1);
  sinn_b.setZero(n_theta_zeta, nboz + 1);

  wmns.setZero(mnmax_nyq); // Note mnmax_nyq instead of mnboz for this one.
  if (asym) {
    wmnc.setZero(mnmax_nyq); // Note mnmax_nyq instead of mnboz for this one.
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
  nu.setZero(n_theta_zeta);
  d_nu_d_theta.setZero(n_theta_zeta);
  d_nu_d_zeta.setZero(n_theta_zeta);
  zeta_Boozer_grid.setZero(n_theta_zeta);
  theta_Boozer_grid.setZero(n_theta_zeta);
  d_Boozer_d_vmec.setZero(n_theta_zeta);
  boozer_jacobian.setZero(n_theta_zeta);


  // This next bit corresponds to transpmn.f in the fortran version.

  // Compute the part of p that is independent of lambda (the
  // covariant source terms in eq (10)). The net p in real space is
  //
  // p(FINAL) = { SUM(m,n)[wmns(m,n)*SIN(arg) + wmnc(m,n)*COS(arg)] - I * lambda } / D
  //
  // WHERE arg = m * theta - n * zeta and D = G + iota * I
  for (jmn = 0; jmn < mnmax_nyq; jmn++) {
    m = xm_nyq[jmn];
    n = xn_nyq[jmn];

    if (m != 0) {
      wmns(jmn) = bsubumnc(jmn, js) / m;
      if (asym) wmnc(jmn) = -bsubumns(jmn, js) / m;
    } else if (n != 0) {
      wmns(jmn) = -bsubvmnc(jmn, js) / n;
      if (asym) wmnc(jmn) = bsubvmns(jmn, js) / n;
    } else {
      // So m = n = 0.
      Boozer_I[js_b] = bsubumnc(jmn, js);
      Boozer_G[js_b] = bsubvmnc(jmn, js);
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
  r.setZero();
  z.setZero();
  lambda.setZero();
  d_lambda_d_theta.setZero();
  d_lambda_d_zeta.setZero();
  //std::cout << "lmns: ";
  for (jmn = 0; jmn < mnmax; jmn++) {
    //std::cout << " " << lmns(jmn, js);
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
  //std::cout << std::endl;
  //if (js_b == 1) std::cout << std::setprecision(15) << "r:" << std::endl << r << std::endl << "z:" << std::endl << z << std::endl << "lambda:" << std::endl << lambda << std::endl;
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
  w.setZero();
  d_w_d_theta.setZero();
  d_w_d_zeta.setZero();
  bmod.setZero();
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
      w[j]           +=  tsin * wmns(jmn);
      d_w_d_theta[j] +=  tcos * wmns(jmn) * m;
      d_w_d_zeta[j]  += -tcos * wmns(jmn) * n;
      bmod[j]        +=  tcos * bmnc(jmn, js);
      if (!asym) continue;
      w[j]           +=  tcos * wmnc(jmn);
      d_w_d_theta[j] += -tsin * wmnc(jmn) * m;
      d_w_d_zeta[j]  +=  tsin * wmnc(jmn) * n;
      bmod[j]        +=  tsin * bmns(jmn, js);
    }
  }
  // End of the part corresponding to vcooords_w() in fortran.

  // The next part is located in harfun.f in the fortran version.
  boozfloat this_iota = iota[js];
  boozfloat one_over_GI = 1.0 / (Boozer_G[js_b] + this_iota * Boozer_I[js_b]);

  // Get total nu from eq (10):
  nu = one_over_GI * (w - Boozer_I[js_b] * lambda);

  // Get Boozer poloidal and toroidal angles from eq (3):
  theta_Boozer_grid = theta_grid + lambda + this_iota * nu;
  zeta_Boozer_grid = zeta_grid + nu;

  // Derivatives of Eq (10) with respect to theta or zeta:
  d_nu_d_zeta  = one_over_GI * (d_w_d_zeta  - Boozer_I[js_b] * d_lambda_d_zeta);
  d_nu_d_theta = one_over_GI * (d_w_d_theta - Boozer_I[js_b] * d_lambda_d_theta);

  // Eq (12):
  d_Boozer_d_vmec = (1.0 + d_lambda_d_theta) * (1.0 + d_nu_d_zeta)
    + (this_iota - d_lambda_d_zeta) * d_nu_d_theta;

  if (false && js_b == 1) {
    std::cout << std::setprecision(15) << "G: " << Boozer_G[js_b];
    std::cout << "  iota: " << this_iota;
    std::cout << "  I: " << Boozer_I[js_b] << std::endl;
    std::cout << "r:" << std::endl << r << std::endl;
    std::cout << "nu:" << std::endl << nu << std::endl;
    std::cout << "theta_Boozer:" << std::endl << theta_Boozer_grid << std::endl;
    std::cout << "zeta_Boozer:" << std::endl << zeta_Boozer_grid << std::endl;
    std::cout << "d_nu_d_theta:" << std::endl << d_nu_d_theta << std::endl;
    std::cout << "d_nu_d_zeta:" << std::endl << d_nu_d_zeta << std::endl;
    std::cout << "d_lambda_d_theta:" << std::endl << d_lambda_d_theta << std::endl;
    std::cout << "d_lambda_d_zeta:" << std::endl << d_lambda_d_zeta << std::endl;
    std::cout << "d_Boozer_d_vmec:" << std::endl << d_Boozer_d_vmec << std::endl;
  }

  // Consider removing this test in the next line for speed:
  //if (d_Boozer_d_vmec.max() * d_Boozer_d_vmec.min() <= 0)
  //  throw std::runtime_error("d_Boozer_d_vmec crosses through 0!");

  // End of the part located in harfun.f in the fortran version.

  // Now come lines from boozer.f
  // The goal now is to evaluate eq (11).

  init_trig(theta_Boozer_grid, zeta_Boozer_grid,
	    cosm_b, sinm_b, cosn_b, sinn_b, mboz, nboz, nfp);

  if (!asym) {
    // Only integrate in theta half-way around
    int i = nzeta * (nu2_b - 1) + 1;
    int imax = i - 1 + nzeta;
    for (int jmn = 0; jmn <= mboz; jmn++) {
      for (j = 0; j < nzeta; j++) {
	cosm_b(j, jmn) = 0.5 * cosm_b(j, jmn);   // theta = 0
	sinm_b(j, jmn) = 0.5 * sinm_b(j, jmn);   // Should be 0
      }
      for (j = i - 1; j < imax; j++) { // i-1 because C is 0-based.
	cosm_b(j, jmn) = 0.5 * cosm_b(j, jmn); // theta = pi
	sinm_b(j, jmn) = 0.5 * sinm_b(j, jmn); // Should be 0
      }
    }
  }

  // Jacobian of (psi, theta_Boozer, zeta_Boozer), where 2 * pi * psi is toroidal flux:
  boozer_jacobian = (Boozer_G[js_b] + this_iota * Boozer_I[js_b]) / (bmod * bmod);

  /*
  bmncb = 0;
  rmncb = 0;
  zmnsb = 0;
  pmnsb = 0;
  gmncb = 0;
  if (asym) {
    bmnsb = 0;
    rmnsb = 0;
    zmncb = 0;
    pmncb = 0;
    gmnsb = 0;
  }
  */
  boozfloat fourier_factor0, fourier_factor;
  // This "fourier_factor0" quantity is called "fac" and "scl" in the
  // fortran version.
  if (asym) {
    fourier_factor0 = 2.0 / (ntheta * nzeta);
  } else {
    fourier_factor0 = 2.0 / ((nu2_b - 1) * nzeta);
    // Equivalently, fac = 4.0 / (ntheta * nzeta)
  }
  // There is a factor of 1/2 for the m=0 element which now appears at
  // the end of surface_solve().
  for (jmn = 0; jmn < mnboz; jmn++) {
    m = xm_b[jmn];
    n = xn_b[jmn];
    abs_n = std::abs(n / nfp);
    if (n < 0) {
      sign = -1;
    } else {
      sign = 1;
    }

    fourier_factor = fourier_factor0;
    if (jmn == 0) fourier_factor *= 0.5;

    for (j = 0; j < n_theta_zeta; j++) {
      tcos = (cosm_b(j, m) * cosn_b(j, abs_n)
	      + sinm_b(j, m) * sinn_b(j, abs_n) * sign)
	* d_Boozer_d_vmec[j] * fourier_factor;

      tsin = (sinm_b(j, m) * cosn_b(j, abs_n)
	      - cosm_b(j, m) * sinn_b(j, abs_n) * sign)
	* d_Boozer_d_vmec[j] * fourier_factor;

      bmnc_b(jmn, js_b) +=  tcos * bmod[j];
      rmnc_b(jmn, js_b) +=  tcos * r[j];
      zmns_b(jmn, js_b) +=  tsin * z[j];
      numns_b(jmn, js_b) += tsin * nu[j];
      gmnc_b(jmn, js_b) +=  tcos * boozer_jacobian[j];
      if (!asym) continue;
      bmns_b(jmn, js_b) +=  tsin * bmod[j];
      rmns_b(jmn, js_b) +=  tsin * r[j];
      zmnc_b(jmn, js_b) +=  tcos * z[j];
      numnc_b(jmn, js_b) += tcos * nu[j];
      gmns_b(jmn, js_b) +=  tsin * boozer_jacobian[j];
    }
  }

  // End of the part from boozer.f

  if (verbose > 0) check_accuracy(js, js_b, bmod,
				  theta_Boozer_grid, zeta_Boozer_grid,
				  cosm_b, cosn_b, sinm_b, sinn_b);

  if (false && js_b == 0) {
    std::ofstream output_file;

    output_file.open("r_0");
    output_file << std::setprecision(15) << r;
    output_file.close();
  }
  if (false && js_b == 1) {
    std::ofstream output_file;

    output_file.open("bmod");
    output_file << std::setprecision(15) << bmod;
    output_file.close();

    output_file.open("cosm_b");
    output_file << std::setprecision(15) << cosm_b;
    output_file.close();

    output_file.open("cosn_b");
    output_file << std::setprecision(15) << cosn_b;
    output_file.close();

    output_file.open("sinm_b");
    output_file << std::setprecision(15) << sinm_b;
    output_file.close();

    output_file.open("sinn_b");
    output_file << std::setprecision(15) << sinn_b;
    output_file.close();

    output_file.open("r");
    output_file << std::setprecision(15) << r;
    output_file.close();

    output_file.open("d_Boozer_d_vmec");
    output_file << std::setprecision(15) << d_Boozer_d_vmec;
    output_file.close();

    output_file.open("xm_b");
    output_file << std::setprecision(15) << xm_b;
    output_file.close();

    output_file.open("xn_b");
    output_file << std::setprecision(15) << xn_b;
    output_file.close();

    std::cout << "bmnc_b:" << std::endl << std::setprecision(15);
    for (jmn = 0; jmn < mnboz; jmn++) std::cout << " " << bmnc_b(jmn, js_b);
    std::cout << std::endl;

    std::cout << "rmnc_b:" << std::endl << std::setprecision(15);
    for (jmn = 0; jmn < mnboz; jmn++) std::cout << " " << rmnc_b(jmn, js_b);
    std::cout << std::endl;

    std::cout << "zmns_b:" << std::endl << std::setprecision(15);
    for (jmn = 0; jmn < mnboz; jmn++) std::cout << " " << zmns_b(jmn, js_b);
    std::cout << std::endl;

    output_file.open("bmnc_b");
    for (jmn = 0; jmn < mnboz; jmn++) output_file << std::setprecision(15) << " " << bmnc_b(jmn, js_b);
    output_file.close();

    output_file.open("gmnc_b");
    for (jmn = 0; jmn < mnboz; jmn++) output_file << std::setprecision(15) << " " << gmnc_b(jmn, js_b);
    output_file.close();

    output_file.open("rmnc_b");
    for (jmn = 0; jmn < mnboz; jmn++) output_file << std::setprecision(15) << " " << rmnc_b(jmn, js_b);
    output_file.close();

    output_file.open("zmns_b");
    for (jmn = 0; jmn < mnboz; jmn++) output_file << std::setprecision(15) << " " << zmns_b(jmn, js_b);
    output_file.close();

    output_file.open("numns_b");
    for (jmn = 0; jmn < mnboz; jmn++) output_file << std::setprecision(15) << " " << numns_b(jmn, js_b);
    output_file.close();
  }
}
