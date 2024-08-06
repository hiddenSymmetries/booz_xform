#include <iostream>
#include <vector>
#include "booz_xform.hpp"
#include "netcdf_writer.hpp"

using namespace booz_xform;

void Booz_xform::write_boozmn(std::string filename) {
  int j;
  if (verbose > 0) std::cout << "Writing output to " << filename << std::endl;
  NetCDFWriter nc(filename);

  // Define dimensions
  dim_id_type mn_modes_dim, mn_mode_dim, compute_surfs_dim, pack_rad_dim, radius_dim;
  int ns_in_plus_1 = ns_in + 1;
  radius_dim = nc.dim("radius", ns_in_plus_1);
  // The boozmn.nc format has these 2 dimensions, mn_mode and
  // mn_modes, which are the same. It's silly, but we'll do the same
  // here for backwards-compatibility.
  mn_mode_dim = nc.dim("mn_mode", mnboz);
  mn_modes_dim = nc.dim("mn_modes", mnboz);
  // The boozmn.nc format has these 2 dimensions, comput_surfs and
  // pack_rad, which are the same. It's silly, but we'll do the same
  // here for backwards-compatibility.
  compute_surfs_dim = nc.dim("comput_surfs", compute_surfs.size());
  pack_rad_dim = nc.dim("pack_rad", compute_surfs.size());

  // Scalars
  std::string long_version = std::string("C++/Python booz_xform ") + version;
  nc.put("version", long_version, "");

  int asym_int = (int) asym;
  nc.put("lasym__logical__", asym_int, "0 if the configuration is stellarator-symmetric, 1 if not", "");

  int ns = ns_in + 1;
  nc.put("ns_b", ns, "The number of radial grid points in the equilibrium before the transformation to Boozer coordinates", "dimensionless");
  nc.put("nfp_b", nfp, "The number of identical field periods", "dimensionless");
  nc.put("mboz_b", mboz, "Maximum poloidal mode number m for which the Fourier amplitudes rmnc, bmnc etc are stored", "dimensionless");
  nc.put("nboz_b", nboz, "Maximum toroidal mode number n for which the Fourier amplitudes rmnc, bmnc etc are stored", "dimensionless");
  nc.put("mnboz_b", mnboz, "The total number of (m,n) pairs for which Fourier amplitudes rmnc, bmnc etc are stored.", "dimensionless");
  nc.put("aspect_b", aspect, "Aspect ratio, if provided", "dimensionless");

  // For rmax_b, rmin_b, betaxis_b, just store 0 for now. It is hard
  // to calculate these without assuming the input is from VMEC. It is
  // better to store something rather that nothing, since some codes
  // that read boozmn.nc format may crash if the fields are not
  // present.
  boozfloat rmax_b = 0, rmin_b = 0, betaxis_b = 0;
  std::string placeholder_str = "Not implemented in this version of booz_xform. Just zero.";
  nc.put("rmax_b", rmax_b, placeholder_str, "dimensionless");
  nc.put("rmin_b", rmin_b, placeholder_str, "dimensionless");
  nc.put("betaxis_b", betaxis_b, placeholder_str, "dimensionless");

  // 1D arrays
  IntVector jlist(ns_b);
  for (j = 0; j < ns_b; j++) jlist[j] = compute_surfs[j] + 2;
  nc.put(compute_surfs_dim, "jlist", jlist, "1-based radial indices of the original vmec solution for which the transformation to Boozer coordinates was computed. 2 corresponds to the first half-grid point.", "dimensionless");
  nc.put(mn_modes_dim, "ixm_b", xm_b, "Poloidal mode numbers m for which the Fourier amplitudes rmnc, bmnc etc are stored", "dimensionless");
  nc.put(mn_modes_dim, "ixn_b", xn_b, "Toroidal mode numbers n for which the Fourier amplitudes rmnc, bmnc etc are stored", "dimensionless");

  // Insert a 0 at the start of several arrays
  Vector iota_b(ns_in + 1), buco_b(ns_in + 1), bvco_b(ns_in + 1);
  iota_b[0] = 0;
  buco_b[0] = 0;
  bvco_b[0] = 0;
  for (j = 0; j < ns_in; j++) {
    iota_b[j + 1] = iota[j];
    buco_b[j + 1] = Boozer_I_all[j];
    bvco_b[j + 1] = Boozer_G_all[j];
  }
  std::string radius_dim_comment = " The radial grid corresponds to all surfaces on which input data were available, and a 0 is prepended";
  nc.put(radius_dim, "iota_b", iota_b, "Rotational transform." + radius_dim_comment, "dimensionless");
  nc.put(radius_dim, "buco_b", buco_b,
	 "Coefficient multiplying grad theta_Boozer in the covariant representation of the magnetic field vector, often denoted I(psi)."
	 + radius_dim_comment, "Tesla * meter");
  nc.put(radius_dim, "bvco_b", bvco_b,
	 "Coefficient multiplying grad zeta_Boozer in the covariant representation of the magnetic field vector, often denoted G(psi)."
	 + radius_dim_comment, "Tesla * meter");

  // For some arrays in the boozmn.nc format, correct profiles are not
  // available here. Just write vectors of zeros for
  // backwards-compatibility of the boozmn.nc files.
  Vector pres_b(ns_in + 1), beta_b(ns_in + 1);
  beta_b.setZero();
  nc.put(radius_dim, "beta_b", pres_b, placeholder_str, "dimensionless");
  Vector phip_dummy(ns_in + 1), chi_dummy(ns_in + 1), phi_dummy(ns_in+1);
  phip_dummy.setZero();
  chi_dummy.setZero();
  phi_dummy.setZero();
  if (phip.size()==0) {
      pres_b.setZero();
      nc.put(radius_dim, "phip_b", phip_dummy, placeholder_str, "Tesla * meter^2");
      nc.put(radius_dim, "chi_b", chi_dummy, placeholder_str, "Tesla * meter^2");
      nc.put(radius_dim, "pres_b", pres_b, placeholder_str, "dimensionless");
      nc.put(radius_dim, "phi_b", phi_dummy, placeholder_str, "Tesla * meter^2");
  } else {
      phip[0] = 0;
      nc.put(radius_dim, "phip_b", phip, "Derivative of toroidal flux (not divided by (2*pi)) with respect to s. This grid generally does not correspond to the radial grid used for other quantities in this file!", "Tesla * meter^2");
      nc.put(radius_dim, "chi_b", chi, "Uniformly spaced grid going from 0 to the boundary poloidal flux (not divided by (2*pi)). This grid generally does not correspond to the radial grid used for other quantities in this file!", "Tesla * meter^2");
      nc.put(radius_dim, "pres_b", pres, "Pressure on full vmec grid. This grid generally does not correspond to the radial grid used for other quantities in this file!", "Pascal");
      nc.put(radius_dim, "phi_b", phi, "Uniformly spaced grid going from 0 to the boundary toroidal flux (not divided by (2*pi)). This grid generally does not correspond to the radial grid used for other quantities in this file!", "Tesla * meter^2");
  }

  // ND arrays for N > 1:
  std::vector<dim_id_type> bmnc_dim;
  bmnc_dim.push_back(mn_modes_dim);
  bmnc_dim.push_back(pack_rad_dim);

  nc.put(bmnc_dim, "bmnc_b", &bmnc_b(0, 0),
	 "cos(m * theta_Boozer - n * zeta_Boozer) Fourier amplitudes of the magnetic field strength", "Tesla");

  nc.put(bmnc_dim, "rmnc_b", &rmnc_b(0, 0),
	 "cos(m * theta_Boozer - n * zeta_Boozer) Fourier amplitudes of the major radius R", "meter");

  nc.put(bmnc_dim, "zmns_b", &zmns_b(0, 0),
	 "sin(m * theta_Boozer - n * zeta_Boozer) Fourier amplitudes of the Cartesian coordinate Z", "meter");

  Matrix pmns_b = -numns_b;
  Matrix pmnc_b; // pmnc_b must be declared outside the "if (asym)" block so it does not go out of scope before nc.write_and_close()
  nc.put(bmnc_dim, "pmns_b", &pmns_b(0, 0),
	 "sin(m * theta_Boozer - n * zeta_Boozer) Fourier amplitudes of the angle difference zeta_VMEC - zeta_Boozer", "dimensionless");

  nc.put(bmnc_dim, "gmn_b", &gmnc_b(0, 0),
	 "cos(m * theta_Boozer - n * zeta_Boozer) Fourier amplitudes of the Boozer coordinate Jacobian (G + iota * I) / B^2", "meter/Tesla");

  if (asym) {
    // Stellarator-asymmetric modes:

    nc.put(bmnc_dim, "bmns_b", &bmns_b(0, 0),
	   "sin(m * theta_Boozer - n * zeta_Boozer) Fourier amplitudes of the magnetic field strength", "Tesla");

    nc.put(bmnc_dim, "rmns_b", &rmns_b(0, 0),
	   "sin(m * theta_Boozer - n * zeta_Boozer) Fourier amplitudes of the major radius R", "meter");

    nc.put(bmnc_dim, "zmnc_b", &zmnc_b(0, 0),
	   "cos(m * theta_Boozer - n * zeta_Boozer) Fourier amplitudes of the Cartesian coordinate Z", "meter");

    pmnc_b = -numnc_b;
    nc.put(bmnc_dim, "pmnc_b", &pmnc_b(0, 0),
	   "cos(m * theta_Boozer - n * zeta_Boozer) Fourier amplitudes of the angle difference zeta_VMEC - zeta_Boozer", "dimensionless");

    nc.put(bmnc_dim, "gmns_b", &gmns_b(0, 0),
	   "sin(m * theta_Boozer - n * zeta_Boozer) Fourier amplitudes of the Boozer coordinate Jacobian (G + iota * I) / B^2", "meter/Tesla");
  }

  // Done defining the NetCDF data.
  nc.write_and_close();

}
