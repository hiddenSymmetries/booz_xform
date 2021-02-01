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
  
  int ns;
  nc.get("ns", ns);
  // The Boozer transformation can be done on the half-grid surfaces,
  // which are 1 fewer than the full-grid surfaces, of which there are
  // ns:
  ns_in = ns - 1;

  nc.get("nfp", nfp);
  nc.get("mpol", mpol);
  nc.get("ntor", ntor);
  nc.get("mnmax", mnmax);
  nc.get("mnmax_nyq", mnmax_nyq);

  Vector iotas;
  iotas.resize(ns);
  iota.resize(ns_in);
  nc.get("iotas", iotas);
  for (j = 0; j < ns_in; j++) iota[j] = iotas[j + 1];

  // By default, prepare to do the Boozer transformation at all
  // half-grid surfaces:
  compute_surfs.resize(ns_in);
  for (j = 0; j < ns_in; j++) compute_surfs[j] = j;
  
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
  Matrix rmnc0(mnmax, ns);
  Matrix zmns0(mnmax, ns);
  Matrix lmns0(mnmax, ns);
  Matrix rmns0(mnmax, ns);
  Matrix zmnc0(mnmax, ns);
  rmnc.resize(mnmax, ns_in);
  zmns.resize(mnmax, ns_in);
  lmns.resize(mnmax, ns_in);
  nc.get("rmnc", rmnc0);
  nc.get("zmns", zmns0);
  nc.get("lmns", lmns0);
  // Keep all but the first radial entry, which is all 0. Selecting
  // the subset of the matrix like this could perhaps be done without
  // a copy, or a copy could be done faster and more elegantly. But
  // for now this works.
  for (j = 0; j < ns_in; j++) {
    for (k = 0; k < mnmax; k++) {
      lmns(k, j) = lmns0(k, j + 1);
    }
  }

  // Nyquist quantities:
  Matrix bmnc0(mnmax_nyq, ns);
  Matrix bsubumnc0(mnmax_nyq, ns);
  Matrix bsubvmnc0(mnmax_nyq, ns);
  bmnc.resize(mnmax_nyq, ns_in);
  bsubumnc.resize(mnmax_nyq, ns_in);
  bsubvmnc.resize(mnmax_nyq, ns_in);
  nc.get("bmnc", bmnc0);
  nc.get("bsubumnc", bsubumnc0);
  nc.get("bsubvmnc", bsubvmnc0);
  // Keep all but the first radial entry, which is all 0:
  for (j = 0; j < ns_in; j++) {
    for (k = 0; k < mnmax_nyq; k++) {
      bmnc(k, j) = bmnc0(k, j + 1);
      bsubumnc(k, j) = bsubumnc0(k, j + 1);
      bsubvmnc(k, j) = bsubvmnc0(k, j + 1);
    }
  }

  if (asym) {
    // Non-Nyquist quantities:
    Matrix lmnc0(mnmax, ns);
    rmns.resize(mnmax, ns_in);
    zmnc.resize(mnmax, ns_in);
    lmnc.resize(mnmax, ns_in);
    nc.get("rmns", rmns0);
    nc.get("zmnc", zmnc0);
    nc.get("lmnc", lmnc0);
    // Keep all but the first radial entry, which is all 0:
    for (j = 0; j < ns_in; j++) {
      for (k = 0; k < mnmax; k++) {
	lmnc(k, j) = lmnc0(k, j + 1);
      }
    }

    // Nyquist quantities:
    Matrix bmns0(mnmax_nyq, ns);
    Matrix bsubumns0(mnmax_nyq, ns);
    Matrix bsubvmns0(mnmax_nyq, ns);
    bmns.resize(mnmax_nyq, ns_in);
    bsubumns.resize(mnmax_nyq, ns_in);
    bsubvmns.resize(mnmax_nyq, ns_in);
    nc.get("bmns", bmns0);
    nc.get("bsubumns", bsubumns0);
    nc.get("bsubvmns", bsubvmns0);
    for (j = 0; j < ns_in; j++) {
      for (k = 0; k < mnmax_nyq; k++) {
	bmns(k, j) = bmns0(k, j + 1);
	bsubumns(k, j) = bsubumns0(k, j + 1);
	bsubvmns(k, j) = bsubvmns0(k, j + 1);
      }
    }
  }
  
  if (verbose > 0) {
    std::cout << "Read ns=" << ns << ", mpol=" << mpol << ", ntor=" << ntor
	      << ", mnmax=" << mnmax << ", mnmax_nyq=" << mnmax_nyq << std::endl;

    if (verbose > 1) {
      std::cout << "iotas = " << iotas << std::endl;
      std::cout << "iota = " << iotas << std::endl;
      std::cout << "xm = " << xm << std::endl;
      std::cout << "xn = " << xn << std::endl;
    }
    
    std::cout << "rmnc, increasing ns index:";
    for (j = 0; j < 4; j++) std::cout << " " << rmnc0(0, j);
    std::cout << std::endl;
    
    std::cout << "rmnc, ncreasing mnmax index:";
    for (j = 0; j < 4; j++) std::cout << " " << rmnc0(j, 0);
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
  
  // We will need sqrt(s) on the full and half grid:
  boozfloat hs = 1.0 / (ns - 1.0);
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
    m = xm[j];
    if (m % 2 == 0) {
      // Even m modes. Just interpolate rmnc and zmns in s.
      for (k = 0; k < ns_in; k++) {
	rmnc(j, k) = 0.5 * (rmnc0(j, k) + rmnc0(j, k + 1));
	zmns(j, k) = 0.5 * (zmns0(j, k) + zmns0(j, k + 1));
	if (asym) {
	  rmns(j, k) = 0.5 * (rmns0(j, k) + rmns0(j, k + 1));
	  zmnc(j, k) = 0.5 * (zmnc0(j, k) + zmnc0(j, k + 1));
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
      for (k = 0; k < ns_in; k++) {
	rmnc(j, k) = 0.5 * (rmnc0(j, k) / sqrt_s_full[k]
			    + rmnc0(j, k + 1) / sqrt_s_full[k + 1]) * sqrt_s_half[k];
	
	zmns(j, k) = 0.5 * (zmns0(j, k) / sqrt_s_full[k]
			    + zmns0(j, k + 1) / sqrt_s_full[k + 1]) * sqrt_s_half[k];
	if (asym) {
	  rmns(j, k) = 0.5 * (rmns0(j, k) / sqrt_s_full[k]
			      + rmns0(j, k + 1) / sqrt_s_full[k + 1]) * sqrt_s_half[k];
	  
	  zmnc(j, k) = 0.5 * (zmnc0(j, k) / sqrt_s_full[k]
			      + zmnc0(j, k + 1) / sqrt_s_full[k + 1]) * sqrt_s_half[k];
	}
      }
      if (m == 1) {
	// For m=1, rmnc/sqrt(s) has a finite value at s = 0. Get this
	// value by extrapolating from the values at the first two
	// full grid points after s = 0.

	// Extrapolate:
	//val_over_sqrt_s_on_axis = 2 * rmnc0(j, 1) / sqrt_s_full[1] - rmnc0(j, 2) / sqrt_s_full[2];
	// Interpolate:
	//rmnc(j, 1) = 0.5 * (val_over_sqrt_s_on_axis + rmnc0(j, 1) / sqrt_s_full[1]) * sqrt_s_half[0];
	
	rmnc(j, 0) = (1.5 * rmnc0(j, 1) / sqrt_s_full[1]
		      - 0.5 * rmnc0(j, 2) / sqrt_s_full[2]) * sqrt_s_half[0];
	
	zmns(j, 0) = (1.5 * zmns0(j, 1) / sqrt_s_full[1]
		      - 0.5 * zmns0(j, 2) / sqrt_s_full[2]) * sqrt_s_half[0];
	if (asym) {
	  rmns(j, 0) = (1.5 * rmns0(j, 1) / sqrt_s_full[1]
			- 0.5 * rmns0(j, 2) / sqrt_s_full[2]) * sqrt_s_half[0];
	
	  zmnc(j, 0) = (1.5 * zmnc0(j, 1) / sqrt_s_full[1]
			- 0.5 * zmnc0(j, 2) / sqrt_s_full[2]) * sqrt_s_half[0];
	}
      }
    }
  }
  /*
  std::ofstream output_file;
  output_file.open("rmnc_1");
  for (j = 0; j < mnmax; j++) output_file << std::setprecision(15) << " " << rmnc(j, 1);
  output_file.close();
  output_file.open("rmnc_2");
  for (j = 0; j < mnmax; j++) output_file << std::setprecision(15) << " " << rmnc(j, 2);
  output_file.close();
  */
  
  // End of radial interpolation.
  
  // Set a guess for the Fourier resolution:
  mboz = 6 * mpol;
  nboz = std::max(2 * ntor - 1, 0);
  
  nc.close();
  
}
