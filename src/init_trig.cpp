#include "init_trig.hpp"

using namespace booz_xform;

/** Evaluate sin(m*theta), cos(m*theta), sin(n*zeta), and cos(n*zeta) for a range
 *  of (m, n, theta, zeta) values.
 *
 *  This function is called trigfunc in the fortran version.
 *
 *  This function assumes that all of the matrices have been allocated
 *  to the appropriate size beforehand.
 *
 *  This function is called both with the 1st two entries (theta_grid, zeta_grid), and
 *  with with the 1st two entries (theta_Boozer_grid, zeta_Boozer_grid).
 */
void booz_xform::init_trig(Vector& theta_grid_x, Vector& zeta_grid_x,
			   Matrix& cosmx, Matrix& sinmx, Matrix& cosnx, Matrix& sinnx,
			   int mmax, int nmax, int nfpx) {
  int j, m, n;
  int n_theta_zeta_x = theta_grid_x.size();

  // Initialize modes with (m,n) = 0 or 1:
  for (j = 0; j < n_theta_zeta_x; j++) {
    cosmx(j, 0) = 1.0;
    cosnx(j, 0) = 1.0;
    cosmx(j, 1) = cos(theta_grid_x[j]);
    sinmx(j, 1) = sin(theta_grid_x[j]);
  }
  if (nmax >= 1) {
    for (j = 0; j < n_theta_zeta_x; j++) {
      cosnx(j, 1) = cos(nfpx * zeta_grid_x[j]);
      sinnx(j, 1) = sin(nfpx * zeta_grid_x[j]);
    }
  }

  // Evaluate the rest of the m values:
  for (m = 2; m <= mmax; m++) {
    for (j = 0; j < n_theta_zeta_x; j++) {
      cosmx(j, m) = cosmx(j, m - 1) * cosmx(j, 1) - sinmx(j, m - 1) * sinmx(j, 1);
      sinmx(j, m) = sinmx(j, m - 1) * cosmx(j, 1) + cosmx(j, m - 1) * sinmx(j, 1);
    }
  }
  
  // Evaluate the rest of the n values:
  for (n = 2; n <= nmax; n++) {
    for (j = 0; j < n_theta_zeta_x; j++) {
      cosnx(j, n) = cosnx(j, n - 1) * cosnx(j, 1) - sinnx(j, n - 1) * sinnx(j, 1);
      sinnx(j, n) = sinnx(j, n - 1) * cosnx(j, 1) + cosnx(j, n - 1) * sinnx(j, 1);
    }
  }
}
