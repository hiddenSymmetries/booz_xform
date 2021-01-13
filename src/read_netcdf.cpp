#include <iostream>
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
  
  iotas.resize(ns, 0.0);
  nc.get("iotas", iotas);

  jlist.resize(ns - 1, 0);
  for (int j = 0; j < ns - 1; j++) jlist[j] = j + 2;
  
  xm.resize(mnmax, 0.0);
  xn.resize(mnmax, 0.0);
  nc.get("xm", xm);
  nc.get("xn", xn);

  xm_nyq.resize(mnmax_nyq, 0.0);
  xn_nyq.resize(mnmax_nyq, 0.0);
  nc.get("xm_nyq", xm_nyq);
  nc.get("xn_nyq", xn_nyq);
  // Maximum values of m_nyq and n_nyq are in the last elements:
  mpol_nyq = xm_nyq[mnmax_nyq - 1];
  ntor_nyq = xn_nyq[mnmax_nyq - 1] / nfp;

  // Non-Nyquist quantities:
  rmnc.resize(mnmax, ns, 0.0);
  zmns.resize(mnmax, ns, 0.0);
  lmns.resize(mnmax, ns, 0.0);
  nc.get("rmnc", rmnc);
  nc.get("zmns", zmns);
  nc.get("lmns", lmns);

  // Nyquist quantities:
  bmnc.resize(mnmax_nyq, ns, 0.0);
  bsubumnc.resize(mnmax_nyq, ns, 0.0);
  bsubvmnc.resize(mnmax_nyq, ns, 0.0);
  nc.get("bmnc", bmnc);
  nc.get("bsubumnc", bsubumnc);
  nc.get("bsubvmnc", bsubvmnc);

  if (asym) {
    // Non-Nyquist quantities:
    rmns.resize(mnmax, ns, 0.0);
    zmnc.resize(mnmax, ns, 0.0);
    lmnc.resize(mnmax, ns, 0.0);
    nc.get("rmns", rmns);
    nc.get("zmnc", zmnc);
    nc.get("lmnc", lmnc);

    // Nyquist quantities:
    bmns.resize(mnmax_nyq, ns, 0.0);
    bsubumns.resize(mnmax_nyq, ns, 0.0);
    bsubvmns.resize(mnmax_nyq, ns, 0.0);
    nc.get("bmns", bmns);
    nc.get("bsubumns", bsubumns);
    nc.get("bsubvmns", bsubvmns);
  }
  
  if (verbose > 0) {
    std::cout << "Read ns=" << ns << ", mpol=" << mpol << ", ntor=" << ntor
	      << ", mnmax=" << mnmax << ", mnmax_nyq=" << mnmax_nyq << std::endl;
    
    std::cout << "iotas = " << iotas << std::endl;
    std::cout << "xm = " << xm << std::endl;
    std::cout << "xn = " << xn << std::endl;
    
    std::cout << "rmnc, increasing ns index:";
    for (int j = 0; j < 4; j++) std::cout << " " << rmnc(0, j);
    std::cout << std::endl;
    
    std::cout << "rmnc, ncreasing mnmax index:";
    for (int j = 0; j < 4; j++) std::cout << " " << rmnc(j, 0);
    std::cout << std::endl;
  }

  // Set a guess for the Fourier resolution:
  mboz = 6 * mpol;
  nboz = std::max(2 * ntor - 1, 0);
  
  nc.close();
  
}
