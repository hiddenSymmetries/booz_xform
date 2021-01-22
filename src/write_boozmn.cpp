#include <iostream>
#include <vector>
#include "booz_xform.hpp"
#include "netcdf_writer.hpp"

using namespace booz_xform;

void Booz_xform::write_boozmn(std::string filename) {
  if (verbose > 0) std::cout << "Writing output to " << filename << std::endl;
  NetCDFWriter nc(filename);

  // Define dimensions
  dim_id_type mn_modes_dim, compute_surfs_dim;
  mn_modes_dim = nc.dim("mn_modes", mnboz);
  compute_surfs_dim = nc.dim("compute_surfs", jlist.size());
  
  // Scalars
  std::string version = "C++/Python booz_xform v0.0.1";
  nc.put("version", version, "");
  
  int asym_int = (int) asym;
  nc.put("lasym__logical__", asym_int, "0 if the configuration is stellarator-symmetric, 1 if not", "");

  nc.put("ns_b", ns, "The number of radial grid points in the equilibrium before the transformation to Boozer coordinates", "dimensionless");
  nc.put("nfp_b", nfp, "The number of identical field periods", "dimensionless");
  nc.put("mboz_b", mboz, "Maximum poloidal mode number m for which the Fourier amplitudes rmnc, bmnc etc are stored", "dimensionless");
  nc.put("nboz_b", nboz, "Maximum toroidal mode number n for which the Fourier amplitudes rmnc, bmnc etc are stored", "dimensionless");
  nc.put("mnboz_b", mnboz, "The total number of (m,n) pairs for which Fourier amplitudes rmnc, bmnc etc are stored.", "dimensionless");

  // 1D arrays
  nc.put(compute_surfs_dim, "jlist", jlist, "1-based radial indices of the original vmec solution for which the transformation to Boozer coordinates was computed. 2 corresponds to the first half-grid point.", "dimensionless");
  nc.put(mn_modes_dim, "ixm_b", xmb, "Poloidal mode numbers m for which the Fourier amplitudes rmnc, bmnc etc are stored", "dimensionless");
  nc.put(mn_modes_dim, "ixn_b", xnb, "Toroidal mode numbers n for which the Fourier amplitudes rmnc, bmnc etc are stored", "dimensionless");

  // ND arrays for N > 1:
  std::vector<dim_id_type> bmnc_dim;
  bmnc_dim.push_back(mn_modes_dim);
  bmnc_dim.push_back(compute_surfs_dim);

  nc.put(bmnc_dim, "bmnc_b", &bmnc_b(0, 0),
	 "cos(m * theta_Boozer - n * zeta_Boozer) Fourier amplitudes of the magnetic field strength", "Tesla");
  
  nc.put(bmnc_dim, "rmnc_b", &rmnc_b(0, 0),
	 "cos(m * theta_Boozer - n * zeta_Boozer) Fourier amplitudes of the major radius R", "meter");
  
  nc.put(bmnc_dim, "zmns_b", &zmns_b(0, 0),
	 "sin(m * theta_Boozer - n * zeta_Boozer) Fourier amplitudes of the Cartesian coordinate Z", "meter");
  
  nc.put(bmnc_dim, "pmns_b", &pmns_b(0, 0),
	 "sin(m * theta_Boozer - n * zeta_Boozer) Fourier amplitudes of the angle difference zeta_VMEC - zeta_Boozer", "dimensionless");
  
  nc.put(bmnc_dim, "gmn_b", &gmnc_b(0, 0),
	 "cos(m * theta_Boozer - n * zeta_Boozer) Fourier amplitudes of the Boozer coordinate Jacobian (G + iota * I) / B^2", "meter/Tesla");
  
  // Done defining the NetCDF data.
  nc.write_and_close();
  
}
