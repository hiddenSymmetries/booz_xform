// Tangent (jvp) and adjoint (vjp) of the transformation to Boozer coordinates.
//
// The derivatives are hand-written and differentiate exactly the discrete
// computation of surface_state() and surface_solve(), with mboz, nboz, the
// VMEC grid and compute_surfs held fixed. They reuse the loop structure of
// the primal code: the Fourier syntheses on the VMEC grid are linear, so
// their tangents are the same loops applied to the tangent coefficients and
// their adjoints are the transposed sums. The only nonlinear parts are the
// pointwise algebra of eqs (10)-(12) and the trigonometric functions of the
// Boozer angles inside the Fourier integrals of eq (11).

#include <cmath>
#include <stdexcept>
#include <string>
#include "booz_xform.hpp"

using namespace booz_xform;

namespace {

  // Return m if it has the expected shape, a zero matrix if m is empty, and
  // throw otherwise.
  Matrix or_zero(const Matrix& m, int rows, int cols, const char* name) {
    if (m.size() == 0) return Matrix::Zero(rows, cols);
    if (m.rows() != rows || m.cols() != cols)
      throw std::runtime_error(std::string(name) + " has shape ("
			       + std::to_string(m.rows()) + ", " + std::to_string(m.cols())
			       + "), expected (" + std::to_string(rows) + ", "
			       + std::to_string(cols) + ")");
    return m;
  }

  Vector or_zero(const Vector& v, int n, const char* name) {
    if (v.size() == 0) return Vector::Zero(n);
    if (v.size() != n)
      throw std::runtime_error(std::string(name) + " has size "
			       + std::to_string(v.size()) + ", expected "
			       + std::to_string(n));
    return v;
  }

  // Bring every member of an InputArrays to the full shape.
  InputArrays complete_inputs(const InputArrays& a, int mnmax, int mnmax_nyq, int ns) {
    InputArrays b;
    b.iota = or_zero(a.iota, ns, "iota");
    b.rmnc = or_zero(a.rmnc, mnmax, ns, "rmnc");
    b.rmns = or_zero(a.rmns, mnmax, ns, "rmns");
    b.zmnc = or_zero(a.zmnc, mnmax, ns, "zmnc");
    b.zmns = or_zero(a.zmns, mnmax, ns, "zmns");
    b.lmnc = or_zero(a.lmnc, mnmax, ns, "lmnc");
    b.lmns = or_zero(a.lmns, mnmax, ns, "lmns");
    b.bmnc = or_zero(a.bmnc, mnmax_nyq, ns, "bmnc");
    b.bmns = or_zero(a.bmns, mnmax_nyq, ns, "bmns");
    b.bsubumnc = or_zero(a.bsubumnc, mnmax_nyq, ns, "bsubumnc");
    b.bsubumns = or_zero(a.bsubumns, mnmax_nyq, ns, "bsubumns");
    b.bsubvmnc = or_zero(a.bsubvmnc, mnmax_nyq, ns, "bsubvmnc");
    b.bsubvmns = or_zero(a.bsubvmns, mnmax_nyq, ns, "bsubvmns");
    return b;
  }

  OutputArrays zero_outputs(int mnboz, int ns_b) {
    OutputArrays o;
    o.bmnc_b.setZero(mnboz, ns_b);
    o.bmns_b.setZero(mnboz, ns_b);
    o.rmnc_b.setZero(mnboz, ns_b);
    o.rmns_b.setZero(mnboz, ns_b);
    o.zmnc_b.setZero(mnboz, ns_b);
    o.zmns_b.setZero(mnboz, ns_b);
    o.numnc_b.setZero(mnboz, ns_b);
    o.numns_b.setZero(mnboz, ns_b);
    o.gmnc_b.setZero(mnboz, ns_b);
    o.gmns_b.setZero(mnboz, ns_b);
    o.Boozer_G.setZero(ns_b);
    o.Boozer_I.setZero(ns_b);
    return o;
  }

  // Drop the members that do not exist for a stellarator-symmetric run, so
  // that the outputs mirror the attributes of Booz_xform.
  void drop_asym_outputs(OutputArrays& o) {
    o.bmns_b.resize(0, 0);
    o.rmns_b.resize(0, 0);
    o.zmnc_b.resize(0, 0);
    o.numnc_b.resize(0, 0);
    o.gmns_b.resize(0, 0);
  }

  void drop_asym_inputs(InputArrays& a) {
    a.rmns.resize(0, 0);
    a.zmnc.resize(0, 0);
    a.lmnc.resize(0, 0);
    a.bmns.resize(0, 0);
    a.bsubumns.resize(0, 0);
    a.bsubvmns.resize(0, 0);
  }

} // namespace

// The derivatives need the grids and trigonometric tables set up by run()
// for the current resolution and surfaces.
void Booz_xform::check_run_state(const char* caller) {
  bool ok = completed && ns_b > 0 && ns_b == compute_surfs.size()
    && mnboz == (2 * nboz + 1) * (mboz - 1) + nboz + 1
    && bmnc_b.rows() == mnboz && bmnc_b.cols() == ns_b
    && n_theta_zeta > 0 && theta_grid.size() == n_theta_zeta
    && cosm.rows() == n_theta_zeta && cosm.cols() == mpol
    && cosm_nyq.rows() == n_theta_zeta && cosm_nyq.cols() == mpol_nyq + 1;
  for (int j = 0; ok && j < ns_b; j++)
    ok = compute_surfs[j] >= 0 && compute_surfs[j] < ns_in;
  if (!ok)
    throw std::runtime_error(std::string(caller) + "() must be called after run(), "
			     "with mboz, nboz and compute_surfs unchanged since then");
}

void Booz_xform::surface_jvp(int js_b, const InputArrays& d, OutputArrays& d_out) {
  SurfaceState st;
  surface_state(js_b, st);
  const int js = st.js;
  int jmn, m, n, abs_n, j, sign;
  boozfloat tcos, tsin;

  // Tangent of transpmn.f:
  Vector d_wmns = Vector::Zero(mnmax_nyq), d_wmnc = Vector::Zero(mnmax_nyq);
  boozfloat d_I = 0, d_G = 0;
  for (jmn = 0; jmn < mnmax_nyq; jmn++) {
    m = xm_nyq[jmn];
    n = xn_nyq[jmn];
    if (m != 0) {
      d_wmns(jmn) = d.bsubumnc(jmn, js) / m;
      if (asym) d_wmnc(jmn) = -d.bsubumns(jmn, js) / m;
    } else if (n != 0) {
      d_wmns(jmn) = -d.bsubvmnc(jmn, js) / n;
      if (asym) d_wmnc(jmn) = d.bsubvmns(jmn, js) / n;
    } else {
      d_I = d.bsubumnc(jmn, js);
      d_G = d.bsubvmnc(jmn, js);
    }
  }
  d_out.Boozer_I[js_b] = d_I;
  d_out.Boozer_G[js_b] = d_G;

  // Tangent of vcoords_rz():
  Vector d_r = Vector::Zero(n_theta_zeta), d_z = Vector::Zero(n_theta_zeta);
  Vector d_lambda = Vector::Zero(n_theta_zeta);
  Vector d_lt = Vector::Zero(n_theta_zeta), d_lz = Vector::Zero(n_theta_zeta);
  for (jmn = 0; jmn < mnmax; jmn++) {
    m = xm[jmn];
    n = xn[jmn];
    abs_n = std::abs(n) / nfp;
    sign = (n < 0) ? -1 : 1;
    for (j = 0; j < n_theta_zeta; j++) {
      tcos = cosm(j, m) * cosn(j, abs_n) + sinm(j, m) * sinn(j, abs_n) * sign;
      tsin = sinm(j, m) * cosn(j, abs_n) - cosm(j, m) * sinn(j, abs_n) * sign;
      d_r[j]      +=  tcos * d.rmnc(jmn, js);
      d_z[j]      +=  tsin * d.zmns(jmn, js);
      d_lambda[j] +=  tsin * d.lmns(jmn, js);
      d_lt[j]     +=  tcos * d.lmns(jmn, js) * m;
      d_lz[j]     += -tcos * d.lmns(jmn, js) * n;
      if (!asym) continue;
      d_r[j]      +=  tsin * d.rmns(jmn, js);
      d_z[j]      +=  tcos * d.zmnc(jmn, js);
      d_lambda[j] +=  tcos * d.lmnc(jmn, js);
      d_lt[j]     += -tsin * d.lmnc(jmn, js) * m;
      d_lz[j]     +=  tsin * d.lmnc(jmn, js) * n;
    }
  }

  // Tangent of vcoords_w():
  Vector d_w = Vector::Zero(n_theta_zeta), d_wt = Vector::Zero(n_theta_zeta);
  Vector d_wz = Vector::Zero(n_theta_zeta), d_bmod = Vector::Zero(n_theta_zeta);
  for (jmn = 0; jmn < mnmax_nyq; jmn++) {
    m = xm_nyq[jmn];
    n = xn_nyq[jmn];
    abs_n = std::abs(n / nfp);
    sign = (n < 0) ? -1 : 1;
    for (j = 0; j < n_theta_zeta; j++) {
      tcos = cosm_nyq(j, m) * cosn_nyq(j, abs_n) + sinm_nyq(j, m) * sinn_nyq(j, abs_n) * sign;
      tsin = sinm_nyq(j, m) * cosn_nyq(j, abs_n) - cosm_nyq(j, m) * sinn_nyq(j, abs_n) * sign;
      d_w[j]    +=  tsin * d_wmns(jmn);
      d_wt[j]   +=  tcos * d_wmns(jmn) * m;
      d_wz[j]   += -tcos * d_wmns(jmn) * n;
      d_bmod[j] +=  tcos * d.bmnc(jmn, js);
      if (!asym) continue;
      d_w[j]    +=  tcos * d_wmnc(jmn);
      d_wt[j]   += -tsin * d_wmnc(jmn) * m;
      d_wz[j]   +=  tsin * d_wmnc(jmn) * n;
      d_bmod[j] +=  tsin * d.bmns(jmn, js);
    }
  }

  // Tangent of harfun.f, eqs (10), (3) and (12):
  const boozfloat d_iota = d.iota[js];
  const boozfloat d_GI = d_G + d_iota * st.I + st.iota * d_I;
  const boozfloat one_over_GI = 1.0 / st.GI;
  Vector d_nu = one_over_GI * (d_w - d_I * st.lambda - st.I * d_lambda - d_GI * st.nu);
  Vector d_nut = one_over_GI * (d_wt - d_I * st.d_lambda_d_theta - st.I * d_lt
				- d_GI * st.d_nu_d_theta);
  Vector d_nuz = one_over_GI * (d_wz - d_I * st.d_lambda_d_zeta - st.I * d_lz
				- d_GI * st.d_nu_d_zeta);
  Vector d_thetaB = d_lambda + d_iota * st.nu + st.iota * d_nu;
  Vector d_zetaB = d_nu;
  Vector d_D = d_lt * (1.0 + st.d_nu_d_zeta) + (1.0 + st.d_lambda_d_theta) * d_nuz
    + (d_iota - d_lz) * st.d_nu_d_theta + (st.iota - st.d_lambda_d_zeta) * d_nut;
  Vector d_jac = d_GI / (st.bmod * st.bmod)
    - 2.0 * st.GI * d_bmod / (st.bmod * st.bmod * st.bmod);

  // Tangent of the Fourier integrals of eq (11). With c = cos(arg) and
  // s = sin(arg), arg = m theta_B - n zeta_B, the integrand c D f has the
  // tangent -s (m d_thetaB - n d_zetaB) D f + c (d_D f + D d_f).
  const Vector& D = st.d_Boozer_d_vmec;
  boozfloat fac, c, s, d_arg, dc, ds;
  for (jmn = 0; jmn < mnboz; jmn++) {
    m = xm_b[jmn];
    n = xn_b[jmn];
    abs_n = std::abs(n / nfp);
    sign = (n < 0) ? -1 : 1;
    fac = fourier_factor(jmn);
    for (j = 0; j < n_theta_zeta; j++) {
      c = st.cosm_b(j, m) * st.cosn_b(j, abs_n) + st.sinm_b(j, m) * st.sinn_b(j, abs_n) * sign;
      s = st.sinm_b(j, m) * st.cosn_b(j, abs_n) - st.cosm_b(j, m) * st.sinn_b(j, abs_n) * sign;
      d_arg = m * d_thetaB[j] - n * d_zetaB[j];
      // Tangents of c D and s D, times fac:
      dc = fac * ((-s * d_arg) * D[j] + c * d_D[j]);
      ds = fac * ((c * d_arg) * D[j] + s * d_D[j]);
      c *= fac * D[j];
      s *= fac * D[j];
      d_out.bmnc_b(jmn, js_b)  += dc * st.bmod[j] + c * d_bmod[j];
      d_out.rmnc_b(jmn, js_b)  += dc * st.r[j] + c * d_r[j];
      d_out.zmns_b(jmn, js_b)  += ds * st.z[j] + s * d_z[j];
      d_out.numns_b(jmn, js_b) += ds * st.nu[j] + s * d_nu[j];
      d_out.gmnc_b(jmn, js_b)  += dc * st.boozer_jacobian[j] + c * d_jac[j];
      if (!asym) continue;
      d_out.bmns_b(jmn, js_b)  += ds * st.bmod[j] + s * d_bmod[j];
      d_out.rmns_b(jmn, js_b)  += ds * st.r[j] + s * d_r[j];
      d_out.zmnc_b(jmn, js_b)  += dc * st.z[j] + c * d_z[j];
      d_out.numnc_b(jmn, js_b) += dc * st.nu[j] + c * d_nu[j];
      d_out.gmns_b(jmn, js_b)  += ds * st.boozer_jacobian[j] + s * d_jac[j];
    }
  }
}

void Booz_xform::surface_vjp(int js_b, const OutputArrays& bar, InputArrays& bar_in) {
  SurfaceState st;
  surface_state(js_b, st);
  const int js = st.js;
  int jmn, m, n, abs_n, j, sign;
  boozfloat tcos, tsin;
  const Vector& D = st.d_Boozer_d_vmec;

  // Adjoint of the Fourier integrals of eq (11). For each real-space field
  // f, accumulate A_f = sum_mn fac (cbar c + sbar s), which gives
  // fbar = D A_f and Dbar = sum_f f A_f, and accumulate the derivative of
  // the integrand with respect to arg = m theta_B - n zeta_B,
  // sum_f f sum_mn fac (sbar c - cbar s) times m or -n.
  Vector A_b = Vector::Zero(n_theta_zeta), A_r = Vector::Zero(n_theta_zeta);
  Vector A_z = Vector::Zero(n_theta_zeta), A_nu = Vector::Zero(n_theta_zeta);
  Vector A_g = Vector::Zero(n_theta_zeta);
  Vector E_theta = Vector::Zero(n_theta_zeta), E_zeta = Vector::Zero(n_theta_zeta);
  boozfloat fac, c, s, e;
  boozfloat cb_b, cb_r, cb_z = 0, cb_nu = 0, cb_g, sb_b = 0, sb_r = 0, sb_z, sb_nu, sb_g = 0;
  for (jmn = 0; jmn < mnboz; jmn++) {
    m = xm_b[jmn];
    n = xn_b[jmn];
    abs_n = std::abs(n / nfp);
    sign = (n < 0) ? -1 : 1;
    fac = fourier_factor(jmn);
    // Cotangents of the cosine and sine coefficients of each field, times fac:
    cb_b = fac * bar.bmnc_b(jmn, js_b);
    cb_r = fac * bar.rmnc_b(jmn, js_b);
    sb_z = fac * bar.zmns_b(jmn, js_b);
    sb_nu = fac * bar.numns_b(jmn, js_b);
    cb_g = fac * bar.gmnc_b(jmn, js_b);
    if (asym) {
      sb_b = fac * bar.bmns_b(jmn, js_b);
      sb_r = fac * bar.rmns_b(jmn, js_b);
      cb_z = fac * bar.zmnc_b(jmn, js_b);
      cb_nu = fac * bar.numnc_b(jmn, js_b);
      sb_g = fac * bar.gmns_b(jmn, js_b);
    }
    for (j = 0; j < n_theta_zeta; j++) {
      c = st.cosm_b(j, m) * st.cosn_b(j, abs_n) + st.sinm_b(j, m) * st.sinn_b(j, abs_n) * sign;
      s = st.sinm_b(j, m) * st.cosn_b(j, abs_n) - st.cosm_b(j, m) * st.sinn_b(j, abs_n) * sign;
      if (!asym) {
	A_b[j]  += cb_b * c;
	A_r[j]  += cb_r * c;
	A_z[j]  += sb_z * s;
	A_nu[j] += sb_nu * s;
	A_g[j]  += cb_g * c;
	e = -s * (cb_b * st.bmod[j] + cb_r * st.r[j] + cb_g * st.boozer_jacobian[j])
	  + c * (sb_z * st.z[j] + sb_nu * st.nu[j]);
      } else {
	A_b[j]  += cb_b * c + sb_b * s;
	A_r[j]  += cb_r * c + sb_r * s;
	A_z[j]  += cb_z * c + sb_z * s;
	A_nu[j] += cb_nu * c + sb_nu * s;
	A_g[j]  += cb_g * c + sb_g * s;
	e = (sb_b * c - cb_b * s) * st.bmod[j]
	  + (sb_r * c - cb_r * s) * st.r[j]
	  + (sb_z * c - cb_z * s) * st.z[j]
	  + (sb_nu * c - cb_nu * s) * st.nu[j]
	  + (sb_g * c - cb_g * s) * st.boozer_jacobian[j];
      }
      E_theta[j] += m * e;
      E_zeta[j]  -= n * e;
    }
  }

  Vector b_bmod = D * A_b;
  Vector b_r = D * A_r;
  Vector b_z = D * A_z;
  Vector b_nu = D * A_nu;
  Vector b_jac = D * A_g;
  Vector b_D = st.bmod * A_b + st.r * A_r + st.z * A_z + st.nu * A_nu
    + st.boozer_jacobian * A_g;
  Vector b_thetaB = D * E_theta;
  Vector b_zetaB = D * E_zeta;

  // Adjoint of boozer_jacobian = GI / bmod^2:
  boozfloat b_GI = (b_jac / (st.bmod * st.bmod)).sum();
  b_bmod -= 2.0 * st.GI * b_jac / (st.bmod * st.bmod * st.bmod);

  // Adjoint of eq (3):
  boozfloat b_iota = (b_thetaB * st.nu).sum();
  Vector b_lambda = b_thetaB;
  b_nu += b_zetaB + st.iota * b_thetaB;

  // Adjoint of eq (12):
  Vector b_lt = b_D * (1.0 + st.d_nu_d_zeta);
  Vector b_nuz = b_D * (1.0 + st.d_lambda_d_theta);
  b_iota += (b_D * st.d_nu_d_theta).sum();
  Vector b_lz = -b_D * st.d_nu_d_theta;
  Vector b_nut = b_D * (st.iota - st.d_lambda_d_zeta);

  // Adjoint of eq (10) and of its theta and zeta derivatives:
  const boozfloat one_over_GI = 1.0 / st.GI;
  Vector b_w = one_over_GI * b_nu;
  Vector b_wt = one_over_GI * b_nut;
  Vector b_wz = one_over_GI * b_nuz;
  b_lambda -= st.I * b_w;
  b_lt -= st.I * b_wt;
  b_lz -= st.I * b_wz;
  boozfloat b_I = -(b_w * st.lambda + b_wt * st.d_lambda_d_theta
		    + b_wz * st.d_lambda_d_zeta).sum();
  b_GI -= (b_w * st.nu + b_wt * st.d_nu_d_theta + b_wz * st.d_nu_d_zeta).sum();

  // Adjoint of GI = G + iota I, plus the direct outputs G and I:
  boozfloat b_G = b_GI + bar.Boozer_G[js_b];
  b_I += st.iota * b_GI + bar.Boozer_I[js_b];
  b_iota += st.I * b_GI;
  bar_in.iota[js] += b_iota;

  // Adjoint of vcoords_rz(), the transposed sums:
  boozfloat acc_c[3], acc_s[3];
  for (jmn = 0; jmn < mnmax; jmn++) {
    m = xm[jmn];
    n = xn[jmn];
    abs_n = std::abs(n) / nfp;
    sign = (n < 0) ? -1 : 1;
    for (int k = 0; k < 3; k++) acc_c[k] = acc_s[k] = 0;
    for (j = 0; j < n_theta_zeta; j++) {
      tcos = cosm(j, m) * cosn(j, abs_n) + sinm(j, m) * sinn(j, abs_n) * sign;
      tsin = sinm(j, m) * cosn(j, abs_n) - cosm(j, m) * sinn(j, abs_n) * sign;
      // rmnc, zmns, lmns:
      acc_c[0] += tcos * b_r[j];
      acc_s[0] += tsin * b_z[j];
      acc_s[1] += tsin * b_lambda[j] + tcos * (m * b_lt[j] - n * b_lz[j]);
      if (!asym) continue;
      // rmns, zmnc, lmnc:
      acc_s[2] += tsin * b_r[j];
      acc_c[1] += tcos * b_z[j];
      acc_c[2] += tcos * b_lambda[j] - tsin * (m * b_lt[j] - n * b_lz[j]);
    }
    bar_in.rmnc(jmn, js) += acc_c[0];
    bar_in.zmns(jmn, js) += acc_s[0];
    bar_in.lmns(jmn, js) += acc_s[1];
    if (!asym) continue;
    bar_in.rmns(jmn, js) += acc_s[2];
    bar_in.zmnc(jmn, js) += acc_c[1];
    bar_in.lmnc(jmn, js) += acc_c[2];
  }

  // Adjoint of vcoords_w() and of transpmn.f:
  boozfloat b_wmns, b_wmnc, b_bmnc, b_bmns;
  for (jmn = 0; jmn < mnmax_nyq; jmn++) {
    m = xm_nyq[jmn];
    n = xn_nyq[jmn];
    abs_n = std::abs(n / nfp);
    sign = (n < 0) ? -1 : 1;
    b_wmns = b_wmnc = b_bmnc = b_bmns = 0;
    for (j = 0; j < n_theta_zeta; j++) {
      tcos = cosm_nyq(j, m) * cosn_nyq(j, abs_n) + sinm_nyq(j, m) * sinn_nyq(j, abs_n) * sign;
      tsin = sinm_nyq(j, m) * cosn_nyq(j, abs_n) - cosm_nyq(j, m) * sinn_nyq(j, abs_n) * sign;
      b_wmns += tsin * b_w[j] + tcos * (m * b_wt[j] - n * b_wz[j]);
      b_bmnc += tcos * b_bmod[j];
      if (!asym) continue;
      b_wmnc += tcos * b_w[j] - tsin * (m * b_wt[j] - n * b_wz[j]);
      b_bmns += tsin * b_bmod[j];
    }
    bar_in.bmnc(jmn, js) += b_bmnc;
    if (asym) bar_in.bmns(jmn, js) += b_bmns;
    if (m != 0) {
      bar_in.bsubumnc(jmn, js) += b_wmns / m;
      if (asym) bar_in.bsubumns(jmn, js) -= b_wmnc / m;
    } else if (n != 0) {
      bar_in.bsubvmnc(jmn, js) -= b_wmns / n;
      if (asym) bar_in.bsubvmns(jmn, js) += b_wmnc / n;
    } else {
      bar_in.bsubumnc(jmn, js) += b_I;
      bar_in.bsubvmnc(jmn, js) += b_G;
    }
  }
}

OutputArrays Booz_xform::jvp(const InputArrays& d_in) {
  check_run_state("jvp");
  InputArrays d = complete_inputs(d_in, mnmax, mnmax_nyq, ns_in);
  OutputArrays d_out = zero_outputs(mnboz, ns_b);
#pragma omp parallel for
  for (int js_b = 0; js_b < ns_b; js_b++) {
    surface_jvp(js_b, d, d_out);
  }
  if (!asym) drop_asym_outputs(d_out);
  return d_out;
}

InputArrays Booz_xform::vjp(const OutputArrays& bar_out) {
  check_run_state("vjp");
  OutputArrays bar;
  bar.bmnc_b = or_zero(bar_out.bmnc_b, mnboz, ns_b, "bmnc_b");
  bar.bmns_b = or_zero(bar_out.bmns_b, mnboz, ns_b, "bmns_b");
  bar.rmnc_b = or_zero(bar_out.rmnc_b, mnboz, ns_b, "rmnc_b");
  bar.rmns_b = or_zero(bar_out.rmns_b, mnboz, ns_b, "rmns_b");
  bar.zmnc_b = or_zero(bar_out.zmnc_b, mnboz, ns_b, "zmnc_b");
  bar.zmns_b = or_zero(bar_out.zmns_b, mnboz, ns_b, "zmns_b");
  bar.numnc_b = or_zero(bar_out.numnc_b, mnboz, ns_b, "numnc_b");
  bar.numns_b = or_zero(bar_out.numns_b, mnboz, ns_b, "numns_b");
  bar.gmnc_b = or_zero(bar_out.gmnc_b, mnboz, ns_b, "gmnc_b");
  bar.gmns_b = or_zero(bar_out.gmns_b, mnboz, ns_b, "gmns_b");
  bar.Boozer_G = or_zero(bar_out.Boozer_G, ns_b, "Boozer_G");
  bar.Boozer_I = or_zero(bar_out.Boozer_I, ns_b, "Boozer_I");

  InputArrays bar_in = complete_inputs(InputArrays(), mnmax, mnmax_nyq, ns_in);
#pragma omp parallel
  {
    // Each thread accumulates into its own arrays, since compute_surfs may
    // contain a surface more than once.
    InputArrays local = complete_inputs(InputArrays(), mnmax, mnmax_nyq, ns_in);
#pragma omp for
    for (int js_b = 0; js_b < ns_b; js_b++) {
      surface_vjp(js_b, bar, local);
    }
#pragma omp critical
    {
      bar_in.iota += local.iota;
      bar_in.rmnc += local.rmnc;
      bar_in.rmns += local.rmns;
      bar_in.zmnc += local.zmnc;
      bar_in.zmns += local.zmns;
      bar_in.lmnc += local.lmnc;
      bar_in.lmns += local.lmns;
      bar_in.bmnc += local.bmnc;
      bar_in.bmns += local.bmns;
      bar_in.bsubumnc += local.bsubumnc;
      bar_in.bsubumns += local.bsubumns;
      bar_in.bsubvmnc += local.bsubvmnc;
      bar_in.bsubvmns += local.bsubvmns;
    }
  }
  if (!asym) drop_asym_inputs(bar_in);
  return bar_in;
}

// The radial interpolation of init_from_vmec(). For a full-grid quantity
// with poloidal mode number m, the half-grid value at index k is
// sum_i weight_i * full(index_i). The arithmetic of the primal code is
// reproduced exactly in init_from_vmec_jvp() below.
namespace {
  struct HalfGrid {
    Vector sqrt_s_full, sqrt_s_half;
    explicit HalfGrid(int ns) {
      boozfloat hs = 1.0 / (ns - 1.0);
      sqrt_s_full.resize(ns);
      sqrt_s_half.resize(ns - 1);
      for (int j = 0; j < ns; j++) sqrt_s_full[j] = sqrt(hs * j);
      for (int j = 0; j < ns - 1; j++) sqrt_s_half[j] = sqrt(hs * (j + 0.5));
      sqrt_s_full[0] = 1.0;
    }
    // Interpolate one row (mode m) of a full-grid array to the half grid.
    void interp(int m, const Matrix& in, int row, Matrix& out) const {
      int ns_in = out.cols();
      if (m % 2 == 0) {
	for (int k = 0; k < ns_in; k++)
	  out(row, k) = 0.5 * (in(row, k) + in(row, k + 1));
      } else {
	for (int k = 0; k < ns_in; k++)
	  out(row, k) = 0.5 * (in(row, k) / sqrt_s_full[k]
			       + in(row, k + 1) / sqrt_s_full[k + 1]) * sqrt_s_half[k];
	if (m == 1)
	  out(row, 0) = (1.5 * in(row, 1) / sqrt_s_full[1]
			 - 0.5 * in(row, 2) / sqrt_s_full[2]) * sqrt_s_half[0];
      }
    }
    // Transpose of interp().
    void interp_transpose(int m, const Matrix& out, int row, Matrix& in) const {
      int ns_in = out.cols();
      in.row(row).setZero();
      if (m % 2 == 0) {
	for (int k = 0; k < ns_in; k++) {
	  in(row, k) += 0.5 * out(row, k);
	  in(row, k + 1) += 0.5 * out(row, k);
	}
      } else {
	for (int k = (m == 1) ? 1 : 0; k < ns_in; k++) {
	  in(row, k) += 0.5 * out(row, k) * sqrt_s_half[k] / sqrt_s_full[k];
	  in(row, k + 1) += 0.5 * out(row, k) * sqrt_s_half[k] / sqrt_s_full[k + 1];
	}
	if (m == 1) {
	  in(row, 1) += 1.5 * out(row, 0) * sqrt_s_half[0] / sqrt_s_full[1];
	  in(row, 2) -= 0.5 * out(row, 0) * sqrt_s_half[0] / sqrt_s_full[2];
	}
      }
    }
  };

  // Half-grid quantities: drop the first radial index.
  Matrix shift(const Matrix& in) {
    return in.rightCols(in.cols() - 1);
  }
  Matrix shift_transpose(const Matrix& out) {
    Matrix in = Matrix::Zero(out.rows(), out.cols() + 1);
    in.rightCols(out.cols()) = out;
    return in;
  }
}

InputArrays Booz_xform::init_from_vmec_jvp(const InputArrays& d_wout) {
  int ns = ns_in + 1;
  if (ns_in < 1 || rmnc.cols() != ns_in)
    throw std::runtime_error("init_from_vmec_jvp() must be called after init_from_vmec()");
  InputArrays d = complete_inputs(d_wout, mnmax, mnmax_nyq, ns);
  InputArrays r = complete_inputs(InputArrays(), mnmax, mnmax_nyq, ns_in);
  r.iota = d.iota.tail(ns_in);
  r.lmns = shift(d.lmns);
  r.lmnc = shift(d.lmnc);
  r.bmnc = shift(d.bmnc);
  r.bmns = shift(d.bmns);
  r.bsubumnc = shift(d.bsubumnc);
  r.bsubumns = shift(d.bsubumns);
  r.bsubvmnc = shift(d.bsubvmnc);
  r.bsubvmns = shift(d.bsubvmns);
  HalfGrid hg(ns);
  for (int j = 0; j < mnmax; j++) {
    hg.interp(xm[j], d.rmnc, j, r.rmnc);
    hg.interp(xm[j], d.zmns, j, r.zmns);
    hg.interp(xm[j], d.rmns, j, r.rmns);
    hg.interp(xm[j], d.zmnc, j, r.zmnc);
  }
  if (!asym) drop_asym_inputs(r);
  return r;
}

InputArrays Booz_xform::init_from_vmec_vjp(const InputArrays& bar_internal) {
  int ns = ns_in + 1;
  if (ns_in < 1 || rmnc.cols() != ns_in)
    throw std::runtime_error("init_from_vmec_vjp() must be called after init_from_vmec()");
  InputArrays b = complete_inputs(bar_internal, mnmax, mnmax_nyq, ns_in);
  InputArrays r = complete_inputs(InputArrays(), mnmax, mnmax_nyq, ns);
  r.iota.tail(ns_in) = b.iota;
  r.lmns = shift_transpose(b.lmns);
  r.lmnc = shift_transpose(b.lmnc);
  r.bmnc = shift_transpose(b.bmnc);
  r.bmns = shift_transpose(b.bmns);
  r.bsubumnc = shift_transpose(b.bsubumnc);
  r.bsubumns = shift_transpose(b.bsubumns);
  r.bsubvmnc = shift_transpose(b.bsubvmnc);
  r.bsubvmns = shift_transpose(b.bsubvmns);
  HalfGrid hg(ns);
  for (int j = 0; j < mnmax; j++) {
    hg.interp_transpose(xm[j], b.rmnc, j, r.rmnc);
    hg.interp_transpose(xm[j], b.zmns, j, r.zmns);
    hg.interp_transpose(xm[j], b.rmns, j, r.rmns);
    hg.interp_transpose(xm[j], b.zmnc, j, r.zmnc);
  }
  if (!asym) drop_asym_inputs(r);
  return r;
}
