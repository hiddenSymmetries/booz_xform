#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include "booz_xform.hpp"
#include "netcdf_reader.hpp"

using namespace booz_xform;

void Booz_xform::init_from_vmec(int ns,
				Vector& iotas,
				Matrix& rmnc0,
				Matrix& rmns0,
				Matrix& zmnc0,
				Matrix& zmns0,
				Matrix& lmnc0,
				Matrix& lmns0,
				Matrix& bmnc0,
				Matrix& bmns0,
				Matrix& bsubumnc0,
				Matrix& bsubumns0,
				Matrix& bsubvmnc0,
				Matrix& bsubvmns0,
        Vector& phip0,
        Vector& chi0,
        Vector& pres0,
        Vector& phi0) {
  int j, k;
  ns_in = ns - 1;

  // Do some validation.

  if (ns < 2) throw std::runtime_error("ns must be at least 2");
  if (nfp < 1) throw std::runtime_error("nfp must be at least 1");
  if (iotas.size() != ns) throw std::runtime_error("iotas.size() is not ns");
  bool skip_phip0 = (phip0.size()==0);
  if (!skip_phip0) {
      if (phip0.size() != ns) throw std::runtime_error("phip0.size() is not ns");
      if (chi0.size() != ns) throw std::runtime_error("chi0.size() is not ns");
      if (pres0.size() != ns) throw std::runtime_error("pres0.size() is not ns");
      if (phi0.size() != ns) throw std::runtime_error("phi0.size() is not ns");
  }
  if (xm.size() != mnmax) throw std::runtime_error("Size of xm is not mnmax");
  if (xn.size() != mnmax) throw std::runtime_error("Size of xn is not mnmax");
  if (xm_nyq.size() != mnmax_nyq) throw std::runtime_error("Size of xm_nyq is not mnmax_nyq");
  if (xn_nyq.size() != mnmax_nyq) throw std::runtime_error("Size of xn_nyq is not mnmax_nyq");

  // if (mnmax != (ntor * 2 + 1) * mpol + ntor + 1) throw std::runtime_error("ntor, mpol, and mnmax are not consistent");
  // if (mnmax_nyq != (ntor_nyq * 2 + 1) * mpol_nyq + ntor_nyq + 1) throw std::runtime_error("ntor, mpol, and mnmax are not consistent");
  if (xm[0] != 0) throw std::runtime_error("xm does not seem correct");
  if (xn[0] != 0) throw std::runtime_error("xn does not seem correct");
  if (xm_nyq[0] != 0) throw std::runtime_error("xm_nyq does not seem correct");
  if (xn_nyq[0] != 0) throw std::runtime_error("xn_nyq does not seem correct");
  // if (xm[mnmax - 1] != mpol) throw std::runtime_error("xm does not seem correct");
  if (xn[mnmax - 1] != nfp * ntor) throw std::runtime_error("xn does not seem correct");
  // if (xm_nyq[mnmax_nyq - 1] != mpol_nyq) throw std::runtime_error("xm_nyq does not seem correct");
  if (xn_nyq[mnmax_nyq - 1] != nfp * ntor_nyq) throw std::runtime_error("xn_nyq does not seem correct");

  if (rmnc0.cols() != ns) throw std::runtime_error("rmnc0 has wrong number of cols");
  if (zmns0.cols() != ns) throw std::runtime_error("zmns0 has wrong number of cols");
  if (lmns0.cols() != ns) throw std::runtime_error("lmns0 has wrong number of cols");
  if (bmnc0.cols() != ns) throw std::runtime_error("bmnc0 has wrong number of cols");
  if (bsubumnc0.cols() != ns) throw std::runtime_error("bsubumnc0 has wrong number of cols");
  if (bsubvmnc0.cols() != ns) throw std::runtime_error("bsubvmnc0 has wrong number of cols");

  if (rmnc0.rows() != mnmax) throw std::runtime_error("rmnc0 has the wrong number of rows");
  if (zmns0.rows() != mnmax) throw std::runtime_error("zmns0 has the wrong number of rows");
  if (lmns0.rows() != mnmax) throw std::runtime_error("lmns0 has the wrong number of rows");
  if (bmnc0.rows() != mnmax_nyq) throw std::runtime_error("bmnc0 has the wrong number of rows");
  if (bsubumnc0.rows() != mnmax_nyq) throw std::runtime_error("bsubumnc0 has the wrong number of rows");
  if (bsubvmnc0.rows() != mnmax_nyq) throw std::runtime_error("bsubvmnc0 has the wrong number of rows");

  if (asym) {
    if (rmns0.cols() != ns) throw std::runtime_error("rmns0 has wrong number of cols");
    if (zmnc0.cols() != ns) throw std::runtime_error("zmnc0 has wrong number of cols");
    if (lmnc0.cols() != ns) throw std::runtime_error("lmnc0 has wrong number of cols");
    if (bmns0.cols() != ns) throw std::runtime_error("bmns0 has wrong number of cols");
    if (bsubumns0.cols() != ns) throw std::runtime_error("bsubumns0 has wrong number of cols");
    if (bsubvmns0.cols() != ns) throw std::runtime_error("bsubvmns0 has wrong number of cols");

    if (rmns0.rows() != mnmax) throw std::runtime_error("rmns0 has the wrong number of rows");
    if (zmnc0.rows() != mnmax) throw std::runtime_error("zmnc0 has the wrong number of rows");
    if (lmnc0.rows() != mnmax) throw std::runtime_error("lmnc0 has the wrong number of rows");
    if (bmns0.rows() != mnmax_nyq) throw std::runtime_error("bmns0 has the wrong number of rows");
    if (bsubumns0.rows() != mnmax_nyq) throw std::runtime_error("bsubumns0 has the wrong number of rows");
    if (bsubvmns0.rows() != mnmax_nyq) throw std::runtime_error("bsubvmns0 has the wrong number of rows");
  }

  for (j = 0; j < compute_surfs.size(); j++) {
    if (compute_surfs[j] < 0) throw std::runtime_error("compute_surfs cannot be negative");
    if (compute_surfs[j] >= ns - 1) throw std::runtime_error("compute_surfs has an entry that is too large for the given ns");
  }

  // Done with validation.

  iota.resize(ns_in);
  for (j = 0; j < ns_in; j++)  {
      iota[j] = iotas[j + 1];
  }
  if (!skip_phip0) {
      phip.resize(ns_in + 1);
      chi.resize(ns_in + 1);
      pres.resize(ns_in + 1);
      phi.resize(ns_in + 1);
      for (j = 0; j < ns_in + 1; j++)  {
          phip[j] = -phip0[j]/(2*pi);
          chi[j] = chi0[j];
          pres[j] = pres0[j];
          phi[j] = phi0[j];
      }
      toroidal_flux = phi[ns - 1];
  } else {
     phip.resize(0);
     chi.resize(0);
     pres.resize(0);
     phi.resize(0);
     toroidal_flux = 0.0;
  }
  // By default, prepare to do the Boozer transformation at all
  // half-grid surfaces:
  compute_surfs.resize(ns_in);
  for (j = 0; j < ns_in; j++) compute_surfs[j] = j;

  // Non-Nyquist quantities:
  rmnc.resize(mnmax, ns_in);
  zmns.resize(mnmax, ns_in);
  lmns.resize(mnmax, ns_in);
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
  bmnc.resize(mnmax_nyq, ns_in);
  bsubumnc.resize(mnmax_nyq, ns_in);
  bsubvmnc.resize(mnmax_nyq, ns_in);
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
    rmns.resize(mnmax, ns_in);
    zmnc.resize(mnmax, ns_in);
    lmnc.resize(mnmax, ns_in);
    // Keep all but the first radial entry, which is all 0:
    for (j = 0; j < ns_in; j++) {
      for (k = 0; k < mnmax; k++) {
	lmnc(k, j) = lmnc0(k, j + 1);
      }
    }

    // Nyquist quantities:
    bmns.resize(mnmax_nyq, ns_in);
    bsubumns.resize(mnmax_nyq, ns_in);
    bsubvmns.resize(mnmax_nyq, ns_in);
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

      std::cout << "rmnc, increasing ns index:";
      for (j = 0; j < 4; j++) std::cout << " " << rmnc0(0, j);
      std::cout << std::endl;

      std::cout << "rmnc, ncreasing mnmax index:";
      for (j = 0; j < 4; j++) std::cout << " " << rmnc0(j, 0);
      std::cout << std::endl;
    }
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
  s_in.resize(ns - 1);
  for (j = 0; j < ns; j++) sqrt_s_full[j] = sqrt(hs * j);
  for (j = 0; j < ns - 1; j++) {
    s_in[j] = hs * (j + 0.5);
    sqrt_s_half[j] = sqrt(s_in[j]);
  }
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

}
