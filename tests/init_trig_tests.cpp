#include "doctest.h"
#include "init_trig.hpp"

using namespace booz_xform;
using doctest::Approx;

/** Verify that the init_trig function correctly populates the
 *  matrices sinm, cosm, sinn, and cosn with sin(m*theta) etc.
 */
TEST_CASE("init_trig") {
  int j, k, m, n;
  int n_theta_zeta = 13;
  Vector theta_grid(n_theta_zeta), zeta_grid(n_theta_zeta);
  int mpol = 5;
  int ntor = 7;
  int nfp = 3;
  for (j = 0; j < n_theta_zeta; j++) {
    // The values of theta and zeta on the grid should not matter, so
    // we just pick some arbitrary values here.
    theta_grid[j] = 0.1 * j;
    zeta_grid[j] = 0.22 * j;
  }
  Matrix cosm, cosn, sinm, sinn;
  cosm.setZero(n_theta_zeta, mpol + 1);
  sinm.setZero(n_theta_zeta, mpol + 1);
  cosn.setZero(n_theta_zeta, ntor + 1);
  sinn.setZero(n_theta_zeta, ntor + 1);

  init_trig(theta_grid, zeta_grid,
	    cosm, sinm, cosn, sinn,
	    mpol, ntor, nfp);

  for (j = 0; j < n_theta_zeta; j++) {
    CAPTURE(j);
    for (m = 0; m <= mpol; m++) {
      CAPTURE(m);
      CHECK(cosm(j, m) == Approx(cos(m * theta_grid[j])));
      CHECK(sinm(j, m) == Approx(sin(m * theta_grid[j])));
    }
  }
	
  for (j = 0; j < n_theta_zeta; j++) {
    CAPTURE(j);
    for (n = 0; n <= ntor; n++) {
      CAPTURE(n);
      CHECK(cosn(j, n) == Approx(cos(n * nfp * zeta_grid[j])));
      CHECK(sinn(j, n) == Approx(sin(n * nfp * zeta_grid[j])));
    }
  }
	
}
