#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include "booz_xform.hpp"
#include "netcdf_reader.hpp"

using namespace booz_xform;

void Booz_xform::read_boozmn(std::string filename) {
  int j;

  if (verbose > 0) std::cout << "About to try reading boozmn netcdf file " << filename << std::endl;
  booz_xform::NetCDFReader nc(filename);

  int asym_int;
  nc.get("lasym__logical__", asym_int);
  asym = (bool) asym_int;

  nc.get("nfp_b", nfp);

  nc.get("mboz_b", mboz);
  nc.get("nboz_b", nboz);
  nc.get("mnboz_b", mnboz);

  xm_b.resize(mnboz);
  xn_b.resize(mnboz);
  nc.get("ixm_b", xm_b);
  nc.get("ixn_b", xn_b);

  int radius = nc.getdim("radius");
  Vector iota_in, Boozer_G_in, Boozer_I_in;
  iota_in.resize(radius);
  Boozer_G_in.resize(radius);
  Boozer_I_in.resize(radius);
  ns_in = radius - 1;
  iota.resize(ns_in);
  Boozer_G_all.resize(ns_in);
  Boozer_I_all.resize(ns_in);
  phi.resize(ns_in+1);
  phip.resize(ns_in+1);
  chi.resize(ns_in+1);
  pres.resize(ns_in+1);
  nc.get("iota_b", iota_in);
  nc.get("bvco_b", Boozer_G_in);
  nc.get("buco_b", Boozer_I_in);
  nc.get("phi_b", phi);
  nc.get("phip_b", phip);
  try {
      nc.get("chi_b", chi);
  } catch (std::runtime_error error) {
  }
  nc.get("pres_b", pres);
  for (j = 0; j < ns_in; j++) {
      iota[j] = iota_in[j+1];
      Boozer_G_all[j] = Boozer_G_in[j+1];
      Boozer_I_all[j] = Boozer_I_in[j+1];
  }

  toroidal_flux = phi[ns_in];
  
  ns_b = nc.getdim("comput_surfs");
  if (verbose > 0) std::cout << "Read mboz=" << mboz << ", nboz=" << nboz <<
		     ", mnboz=" << mnboz << ", ns_b=" << ns_b << std::endl;

  compute_surfs.resize(ns_b);
  nc.get("jlist", compute_surfs);
  // Fortran jlist is 1-based and includes an extra 1 from the 0 at
  // the first entry of half-grid quantities. compute_surfs is 0-based.
  for (j = 0; j < ns_b; j++) compute_surfs[j] -= 2;
  if (verbose > 0) std::cout << "Read compute_surfs=" << compute_surfs << std::endl;

  bmnc_b.resize(mnboz, ns_b);
  rmnc_b.resize(mnboz, ns_b);
  zmns_b.resize(mnboz, ns_b);
  numns_b.resize(mnboz, ns_b);
  gmnc_b.resize(mnboz, ns_b);
  nc.get("bmnc_b", bmnc_b);
  nc.get("rmnc_b", rmnc_b);
  nc.get("zmns_b", zmns_b);
  nc.get("pmns_b", numns_b);
  numns_b = -numns_b; // p in the boozmn format is -nu in this new booz_xform.
  nc.get("gmn_b", gmnc_b);

  if (asym) {
    bmns_b.resize(mnboz, ns_b);
    rmns_b.resize(mnboz, ns_b);
    zmnc_b.resize(mnboz, ns_b);
    numnc_b.resize(mnboz, ns_b);
    gmns_b.resize(mnboz, ns_b);
    nc.get("bmns_b", bmns_b);
    nc.get("rmns_b", rmns_b);
    nc.get("zmnc_b", zmnc_b);
    nc.get("pmnc_b", numnc_b);
    numnc_b = -numnc_b; // p in the boozmn format is -nu in this new booz_xform.
    nc.get("gmns_b", gmns_b);

  } else {
    // Stellarator-symmetric.

    bmns_b.resize(0, 0);
    rmns_b.resize(0, 0);
    zmnc_b.resize(0, 0);
    numnc_b.resize(0, 0);
    gmns_b.resize(0, 0);
  }

  // Set up the s grid for output quantites.
  int ns_vmec;
  nc.get("ns_b", ns_vmec);
  boozfloat hs = 1.0 / (ns_vmec - 1.0);
  s_b.resize(ns_b);
  for (j = 0; j < ns_b; j++) s_b[j] = hs * (compute_surfs[j] + 0.5);
  if (verbose > 0) std::cout << "s_b=" << s_b << std::endl;

  nc.close();

}
