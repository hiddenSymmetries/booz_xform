#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include "booz_xform.hpp"
#include "netcdf_reader.hpp"

using namespace booz_xform;

void Booz_xform::read_wout(std::string filename, bool flux) {
  int j, k;
  if (verbose > 0) std::cout << "About to try reading VMEC wout file " << filename << std::endl;
  booz_xform::NetCDFReader nc(filename);

  int asym_int;
  nc.get("lasym__logical__", asym_int);
  asym = (bool) asym_int;

  int ns;
  nc.get("ns", ns);

  nc.get("nfp", nfp);
  nc.get("mpol", mpol);
  nc.get("ntor", ntor);
  nc.get("mnmax", mnmax);
  nc.get("mnmax_nyq", mnmax_nyq);
  nc.get("aspect", aspect);

  Vector phip0;
  Vector chi0;
  Vector pres0;
  Vector phi0;
  if (flux) {
      phip0.resize(ns);
      chi0.resize(ns);
      pres0.resize(ns);
      phi0.resize(ns);
      nc.get("chi", chi0);
      nc.get("phipf", phip0);
      nc.get("pres",pres0);
      nc.get("phi", phi0);
  }

  Vector iotas;
  iotas.resize(ns);

  nc.get("iotas", iotas);

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
  Matrix rmns0;
  Matrix zmnc0;
  Matrix lmnc0;
  nc.get("rmnc", rmnc0);
  nc.get("zmns", zmns0);
  nc.get("lmns", lmns0);

  // Nyquist quantities:
  Matrix bmnc0(mnmax_nyq, ns);
  Matrix bsubumnc0(mnmax_nyq, ns);
  Matrix bsubvmnc0(mnmax_nyq, ns);
  Matrix bmns0;
  Matrix bsubumns0;
  Matrix bsubvmns0;
  nc.get("bmnc", bmnc0);
  nc.get("bsubumnc", bsubumnc0);
  nc.get("bsubvmnc", bsubvmnc0);

  if (asym) {
    // Non-Nyquist quantities:
    rmns0.resize(mnmax, ns);
    zmnc0.resize(mnmax, ns);
    lmnc0.resize(mnmax, ns);
    nc.get("rmns", rmns0);
    nc.get("zmnc", zmnc0);
    nc.get("lmnc", lmnc0);

    // Nyquist quantities:
    bmns0.resize(mnmax_nyq, ns);
    bsubumns0.resize(mnmax_nyq, ns);
    bsubvmns0.resize(mnmax_nyq, ns);
    nc.get("bmns", bmns0);
    nc.get("bsubumns", bsubumns0);
    nc.get("bsubvmns", bsubvmns0);
  }

  nc.close();

  if (flux) {
      init_from_vmec(ns, iotas, rmnc0, rmns0, zmnc0, zmns0,
           lmnc0, lmns0, bmnc0, bmns0,
           bsubumnc0, bsubumns0, bsubvmnc0, bsubvmns0, phip0, chi0, pres0, phi0);
  } else {
      init_from_vmec(ns, iotas, rmnc0, rmns0, zmnc0, zmns0,
           lmnc0, lmns0, bmnc0, bmns0,
           bsubumnc0, bsubumns0, bsubvmnc0, bsubvmns0);
  }

}
