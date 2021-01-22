#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include "booz_xform.hpp"
#include "netcdf_reader.hpp"

using namespace booz_xform;

void Booz_xform::read_boozmn(std::string filename) {
  if (verbose > 0) std::cout << "About to try reading boozmn netcdf file " << filename << std::endl;
  booz_xform::NetCDFReader nc(filename);

  nc.get("mboz_b", mboz);
  nc.get("nboz_b", nboz);

  if (verbose > 0) std::cout << "Read mboz=" << mboz << " , nboz=" << nboz << std::endl;
  
  nc.close();
  
}

void Booz_xform::read_wout(std::string filename) {
  int j, k;
  if (verbose > 0) std::cout << "About to try reading VMEC wout file " << filename << std::endl;
  booz_xform::NetCDFReader nc(filename);

  int asym_int;
  nc.get("lasym__logical__", asym_int);
  asym = (bool) asym_int;
  
  nc.get("nfp", nfp);
  nc.get("ns", ns);
  nc.get("mpol", mpol);
  nc.get("ntor", ntor);
  nc.get("mnmax", mnmax);
  nc.get("mnmax_nyq", mnmax_nyq);
  
  iotas.resize(ns);
  nc.get("iotas", iotas);

  jlist.resize(ns - 1);
  for (j = 0; j < ns - 1; j++) jlist[j] = j + 2;
  
  xm.resize(mnmax);
  xn.resize(mnmax);
  nc.get("xm", xm);
  nc.get("xn", xn);

  xm_nyq.resize(mnmax_nyq);
  xn_nyq.resize(mnmax_nyq);
  nc.get("xm_nyq", xm_nyq);
  nc.get("xn_nyq", xn_nyq);
  // Maximum values of m_nyq and n_nyq are in the last elements:
  mpol_nyq = xm_nyq[mnmax_nyq - 1];
  ntor_nyq = xn_nyq[mnmax_nyq - 1] / nfp;

  // Non-Nyquist quantities:
  rmnc.resize(mnmax, ns);
  zmns.resize(mnmax, ns);
  lmns.resize(mnmax, ns);
  nc.get("rmnc", rmnc);
  nc.get("zmns", zmns);
  nc.get("lmns", lmns);

  // Nyquist quantities:
  bmnc.resize(mnmax_nyq, ns);
  bsubumnc.resize(mnmax_nyq, ns);
  bsubvmnc.resize(mnmax_nyq, ns);
  nc.get("bmnc", bmnc);
  nc.get("bsubumnc", bsubumnc);
  nc.get("bsubvmnc", bsubvmnc);

  if (asym) {
    // Non-Nyquist quantities:
    rmns.resize(mnmax, ns);
    zmnc.resize(mnmax, ns);
    lmnc.resize(mnmax, ns);
    nc.get("rmns", rmns);
    nc.get("zmnc", zmnc);
    nc.get("lmnc", lmnc);

    // Nyquist quantities:
    bmns.resize(mnmax_nyq, ns);
    bsubumns.resize(mnmax_nyq, ns);
    bsubvmns.resize(mnmax_nyq, ns);
    nc.get("bmns", bmns);
    nc.get("bsubumns", bsubumns);
    nc.get("bsubvmns", bsubvmns);
  }
  
  if (verbose > 0) {
    std::cout << "Read ns=" << ns << ", mpol=" << mpol << ", ntor=" << ntor
	      << ", mnmax=" << mnmax << ", mnmax_nyq=" << mnmax_nyq << std::endl;

    if (verbose > 1) {
      std::cout << "iotas = " << iotas << std::endl;
      std::cout << "xm = " << xm << std::endl;
      std::cout << "xn = " << xn << std::endl;
    }
    
    std::cout << "rmnc, increasing ns index:";
    for (j = 0; j < 4; j++) std::cout << " " << rmnc(0, j);
    std::cout << std::endl;
    
    std::cout << "rmnc, ncreasing mnmax index:";
    for (j = 0; j < 4; j++) std::cout << " " << rmnc(j, 0);
    std::cout << std::endl;
  }

  // Handle radial interpolation of the full-grid quantities rmnc,
  // zmns, rmns, zmnc to the half mesh. In the fortran version, this is done inside
  // boozer_coords.f, vcoords.f::vcoords_rz(), and booz_jac.f::booz_rzhalf().

  // The radial interpolation for R and Z is done as follows.  For
  // "even parity" modes, meaning those with even m, we just take the
  // mean of values on the two surrounding full-mesh points. For "odd
  // parity" modes, meaning those with odd m, we linearly interpolate
  // rmnc/sqrt(s) and zmns/sqrt(s) in s. Care must be taken for the
  // m=1 modes, since rmnc/sqrt(s) and zmns/sqrt(s) have a finite value
  // at s=0 for the m=1 modes. For these modes, the value at s=0 is
  // found by linear extrapolation using the next two full-grid points.
  
  // Copy the full-grid arrays:
  Matrix rmnc_orig, rmns_orig, zmnc_orig, zmns_orig;
  rmnc_orig.resize(mnmax, ns);
  zmns_orig.resize(mnmax, ns);
  rmnc_orig = rmnc;
  zmns_orig = zmns;
  if (asym) {
    rmns_orig.resize(mnmax, ns);
    zmnc_orig.resize(mnmax, ns);
    rmns_orig = rmns;
    zmnc_orig = zmnc;
  }
  
  // We will need sqrt(s) on the full and half grid:
  hs = 1.0 / (ns - 1.0);
  Vector sqrt_s_full, sqrt_s_half;
  sqrt_s_full.resize(ns);
  sqrt_s_half.resize(ns - 1);
  for (j = 0; j < ns; j++) sqrt_s_full[j] = sqrt(hs * j);
  for (j = 0; j < ns - 1; j++) sqrt_s_half[j] = sqrt(hs * (j + 0.5));
  // To avoid divide-by-zero when we divide by sqrt_s:
  sqrt_s_full[0] = 1.0;

  // Do the interpolation for each (m,n) mode:
  int m;
  for (j = 0; j < mnmax; j++) {
    // Zero out the column that is always zero for half-grid quantities:
    rmnc(j, 0) = 0.0;
    zmns(j, 0) = 0.0;
    if (asym) {
      rmns(j, 0) = 0.0;
      zmnc(j, 0) = 0.0;
    }
    
    m = xm[j];
    if (m % 2 == 0) {
      // Even m modes. Just interpolate rmnc and zmns in s.
      for (k = 0; k < ns - 1; k++) {
	rmnc(j, k + 1) = 0.5 * (rmnc_orig(j, k) + rmnc_orig(j, k + 1));
	zmns(j, k + 1) = 0.5 * (zmns_orig(j, k) + zmns_orig(j, k + 1));
	if (asym) {
	  rmns(j, k + 1) = 0.5 * (rmns_orig(j, k) + rmns_orig(j, k + 1));
	  zmnc(j, k + 1) = 0.5 * (zmnc_orig(j, k) + zmnc_orig(j, k + 1));
	}
      }
    } else {
      // Odd m modes. Interpolate rmnc/sqrt(s) and zmns/sqrt(s) in s.
      // For m > 0, rmnc, zmns etc vanish exactly at s=0. When we
      // divide by the "wrong" value sqrt_s_full=1.0 at s=0, the fact
      // that rmnc(s=0)=0 means that we get a ratio rmnc/sqrt(s)=0,
      // which is correct for modes with m > 1. However this is not
      // correct for the modes with m=1. We will correct these modes a
      // few lines later.
      for (k = 0; k < ns - 1; k++) {
	rmnc(j, k + 1) = 0.5 * (rmnc_orig(j, k) / sqrt_s_full[k]
				+ rmnc_orig(j, k + 1) / sqrt_s_full[k + 1]) * sqrt_s_half[k];
	
	zmns(j, k + 1) = 0.5 * (zmns_orig(j, k) / sqrt_s_full[k]
				+ zmns_orig(j, k + 1) / sqrt_s_full[k + 1]) * sqrt_s_half[k];
	if (asym) {
	  rmns(j, k + 1) = 0.5 * (rmns_orig(j, k) / sqrt_s_full[k]
				  + rmns_orig(j, k + 1) / sqrt_s_full[k + 1]) * sqrt_s_half[k];
	  
	  zmnc(j, k + 1) = 0.5 * (zmnc_orig(j, k) / sqrt_s_full[k]
				  + zmnc_orig(j, k + 1) / sqrt_s_full[k + 1]) * sqrt_s_half[k];
	}
      }
      if (m == 1) {
	// For m=1, rmnc/sqrt(s) has a finite value at s = 0. Get this
	// value by extrapolating from the values at the first two
	// full grid points after s = 0.

	// Extrapolate:
	//val_over_sqrt_s_on_axis = 2 * rmnc_orig(j, 1) / sqrt_s_full[1] - rmnc_orig(j, 2) / sqrt_s_full[2];
	// Interpolate:
	//rmnc(j, 1) = 0.5 * (val_over_sqrt_s_on_axis + rmnc_orig(j, 1) / sqrt_s_full[1]) * sqrt_s_half[0];
	
	rmnc(j, 1) = (1.5 * rmnc_orig(j, 1) / sqrt_s_full[1]
		      - 0.5 * rmnc_orig(j, 2) / sqrt_s_full[2]) * sqrt_s_half[0];
	
	zmns(j, 1) = (1.5 * zmns_orig(j, 1) / sqrt_s_full[1]
		      - 0.5 * zmns_orig(j, 2) / sqrt_s_full[2]) * sqrt_s_half[0];
	if (asym) {
	  rmns(j, 1) = (1.5 * rmns_orig(j, 1) / sqrt_s_full[1]
			- 0.5 * rmns_orig(j, 2) / sqrt_s_full[2]) * sqrt_s_half[0];
	
	  zmnc(j, 1) = (1.5 * zmnc_orig(j, 1) / sqrt_s_full[1]
			- 0.5 * zmnc_orig(j, 2) / sqrt_s_full[2]) * sqrt_s_half[0];
	}
      }
    }
  }
  std::ofstream output_file;
  output_file.open("rmnc_1");
  for (j = 0; j < mnmax; j++) output_file << std::setprecision(15) << " " << rmnc(j, 1);
  output_file.close();
  output_file.open("rmnc_2");
  for (j = 0; j < mnmax; j++) output_file << std::setprecision(15) << " " << rmnc(j, 2);
  output_file.close();
  // End of radial interpolation.
  
  // Set a guess for the Fourier resolution:
  mboz = 6 * mpol;
  nboz = std::max(2 * ntor - 1, 0);
  
  nc.close();
  
}
