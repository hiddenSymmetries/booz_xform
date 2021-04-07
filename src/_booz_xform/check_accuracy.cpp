#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#ifdef OPENMP
#include <omp.h>
#endif
#include "booz_xform.hpp"
#include "init_trig.hpp"

using namespace booz_xform;

/** Check the accuracy of the conversion to Boozer coordinates by comparing |B| at 4 points.
 *
 *  This subroutine corresponds to parts of boozer_coords.f and harfun.f::modbooz().
 */
void Booz_xform::check_accuracy(int js, int js_b,
				Vector& bmod, Vector& theta_Boozer_grid, Vector& zeta_Boozer_grid,
				Matrix& cosm_b, Matrix& cosn_b, Matrix& sinm_b, Matrix& sinn_b) {

  const int n_check = 4;
  Vector bmod_vmec(n_check), bmod_boozer(n_check), bmod_err(n_check);
  Vector theta_b_test(n_check), zeta_b_test(n_check);
  int index, thread;

  // Recall nu2_b = ntheta / 2 + 1 regardless of asym.
  // int nv2_b = nzeta / 2 + 1; // Index of zeta = pi. Note nzeta is even so there is no rounding.

  // First, store |B|, theta_Boozer, and zeta_Boozer at the 4 points
  // given by theta_vmec = (0, pi) and zeta_vmec = (0, pi / nfp):
  
  // zeta is the fast dimension, with range nzeta.
  
  // (theta = 0, zeta = 0):
  index = 0;
  bmod_vmec[0] = bmod[index];
  theta_b_test[0] = theta_Boozer_grid[index];
  zeta_b_test[0] = zeta_Boozer_grid[index];
  assert (std::abs(theta_grid[index] - 0) < 1.0e-13);
  assert (std::abs(zeta_grid[index] - 0) < 1.0e-13);
  
  // (theta = pi, zeta = 0):
  index = (ntheta / 2) * nzeta; // Note ntheta is even so there is no rounding.
  bmod_vmec[1] = bmod[index];
  theta_b_test[1] = theta_Boozer_grid[index];
  zeta_b_test[1] = zeta_Boozer_grid[index];
  assert (std::abs(theta_grid[index] - pi) < 1.0e-13);
  assert (std::abs(zeta_grid[index] - 0) < 1.0e-13);

  // (theta = 0, zeta = pi / nfp):
  index = nzeta / 2;
  // Above, for nonaxisymmetry, nzeta is even so there is no rounding.
  // For axisymmetry, nzeta = 1 so integer division gives index = 0;
  bmod_vmec[2] = bmod[index];
  theta_b_test[2] = theta_Boozer_grid[index];
  zeta_b_test[2] = zeta_Boozer_grid[index];
  assert (std::abs(theta_grid[index] - 0) < 1.0e-13);
  if (nzeta > 1) {
    assert (std::abs(zeta_grid[index] - pi / nfp) < 1.0e-13);
  } else {
    assert (std::abs(zeta_grid[index]) < 1.0e-13);
  }

  // (theta = 0, zeta = pi / nfp):
  index = (ntheta / 2) * nzeta + (nzeta / 2);
  // Above, for nonaxisymmetry, nzeta is even so there is no rounding.
  // For axisymmetry, nzeta = 1 so integer division gives index = 0;
  // ntheta is always even so there is no rounding in ntheta / 2.
  bmod_vmec[3] = bmod[index];
  theta_b_test[3] = theta_Boozer_grid[index];
  zeta_b_test[3] = zeta_Boozer_grid[index];
  assert (std::abs(theta_grid[index] - pi) < 1.0e-13);
  if (nzeta > 1) {
    assert (std::abs(zeta_grid[index] - pi / nfp) < 1.0e-13);
  } else {
    assert (std::abs(zeta_grid[index]) < 1.0e-13);
  }

  // End of parts from boozer_coords.f.
  // Now comes code equivalent to harfun.f::modbooz().

  // The goal is to sum up the bmnc_b (and possibly bmns_b) series to
  // reconstruct |B| at the 4 special points.
  
  // For the next few steps we will re-use cosm_b, sinm_b, cosn_b, and
  // sinn_b. They are larger than necessary, but this is not a
  // problem (except perhaps for memory locality?)
  init_trig(theta_b_test, zeta_b_test,
	    cosm_b, sinm_b, cosn_b, sinn_b,
	    mboz, nboz, nfp);

  int j, jmn, m, n, abs_n, sign;
  boozfloat tcos, tsin;
  for (j = 0; j < n_check; j++) bmod_boozer[j] = 0.0;
  for (jmn = 0; jmn < mnboz; jmn++) {
    m = xm_b[jmn];
    n = xn_b[jmn];
    abs_n = std::abs(n / nfp);
    if (n < 0) {
      sign = -1;
    } else {
      sign = 1;
    }

    for (j = 0; j < n_check; j++) {
      tcos = cosm_b(j, m) * cosn_b(j, abs_n) + sinm_b(j, m) * sinn_b(j, abs_n) * sign;
      bmod_boozer[j] += tcos * bmnc_b(jmn, js_b);
      if (!asym) continue;
      tsin = sinm_b(j, m) * cosn_b(j, abs_n) - cosm_b(j, m) * sinn_b(j, abs_n) * sign;
      bmod_boozer[j] += tsin * bmns_b(jmn, js_b);
    }
  }

  thread = 0;
#ifdef OPENMP
  thread = omp_get_thread_num();
#endif
  
  // Here we copy the error definition used in boozer_coords.f, which is
  // err = ABS(bmodb - bmodv)/MAX(bmodb,bmodv)
  for (j = 0; j < n_check; j++) {
    bmod_err[j] = std::abs(bmod_boozer[j] - bmod_vmec[j])
      / std::max(bmod_boozer[j], bmod_vmec[j]);
  }
  // Only 1 thread should print at a time:
  #pragma omp critical
  {
  std::cout << std::setprecision(3) << std::scientific
	    << std::setw(4) << thread << std::setw(6) << js_b
	    << std::setw(4) << js
	    << "   0" << std::setprecision(3) << std::scientific
	    << std::setw(11) << bmod_vmec[0]
	    << std::setw(11) << bmod_boozer[0]
	    << std::setw(11) << bmod_err[0]
	    << std::setw(11) << bmod_vmec[1]
	    << std::setw(11) << bmod_boozer[1]
	    << std::setw(11) << bmod_err[1]
	    << std::endl;
  std::cout << "                pi" << std::setprecision(3)
	    << std::setw(11) << bmod_vmec[2]
	    << std::setw(11) << bmod_boozer[2]
	    << std::setw(11) << bmod_err[2]
	    << std::setw(11) << bmod_vmec[3]
	    << std::setw(11) << bmod_boozer[3]
	    << std::setw(11) << bmod_err[3]
            << std::endl;
  }
}
