#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "omori_kernel.h"
using namespace Rcpp;

//' Sequential MAP or Bernoulli process labels for bivariate ETAS
//'
//' Walks events in time order. Rows with \code{assignable==0} keep
//' \code{process_id_init}. Assignable rows use already-assigned parents.
//' If \code{sample_bernoulli} is false, assign
//' \code{argmax_k \lambda_k(t_i | H_{t_i})}. If true, draw
//' \code{Z_i ~ Bern(\lambda_1 / (\lambda_0 + \lambda_1))}.
//'
//' @keywords internal
// [[Rcpp::export]]
IntegerVector etas_bivariate_sequential_map_cpp(
    NumericVector t, NumericVector x, NumericVector y, NumericVector mag,
    IntegerVector assignable,
    IntegerVector process_id_init,
    NumericVector W_val_0, NumericVector W_val_1,
    double mu_0, double mu_1,
    double A_00, double alpha_m_00,
    double A_11, double alpha_m_11,
    double A_01, double alpha_m_01,
    double A_10, double alpha_m_10,
    double cc, double p, double D, double gamma_par, double q,
    double m0,
    double areaS_0, double areaS_1,
    double t_trunc = -1.0,
    bool sample_bernoulli = false) {

  const int n = t.size();
  IntegerVector pid(n);
  if (n == 0) return pid;

  const double pi_val = 3.14159265358979323846;
  const bool do_trunc = (t_trunc > 0.0);
  if (areaS_0 <= 0.0) areaS_0 = 1.0;
  if (areaS_1 <= 0.0) areaS_1 = 1.0;

  const double mu_base_0 = mu_0 / areaS_0;
  const double mu_base_1 = mu_1 / areaS_1;
  const double omori_pref = omori_time_prefactor(p, cc, t_trunc, do_trunc);
  const double base_const = omori_pref * (q - 1.0) / pi_val;

  double A_mat[2][2], alpha_mat[2][2];
  A_mat[0][0] = A_00;   alpha_mat[0][0] = alpha_m_00;
  A_mat[0][1] = A_01;   alpha_mat[0][1] = alpha_m_01;
  A_mat[1][0] = A_10;   alpha_mat[1][0] = alpha_m_10;
  A_mat[1][1] = A_11;   alpha_mat[1][1] = alpha_m_11;

  const double* pt = t.begin();
  const double* px = x.begin();
  const double* py = y.begin();
  const double* pmag = mag.begin();
  const double* pW0 = W_val_0.begin();
  const double* pW1 = W_val_1.begin();
  const int* passign = assignable.begin();
  const int* pinit = process_id_init.begin();

  std::vector<double> d_par(n), pk0(n), pk1(n);
  const double inv_cc = 1.0 / cc;

  auto fill_parent_const = [&](int j) {
    const double dm = pmag[j] - m0;
    const double dj = D * std::exp(gamma_par * dm);
    d_par[j] = dj;
    const double inv_dj = 1.0 / dj;
    const int l = pid[j];
    const double A0l = A_mat[0][l];
    const double A1l = A_mat[1][l];
    const double k0 = (A0l < 1e-20) ? 0.0 : A0l * std::exp(alpha_mat[0][l] * dm);
    const double k1 = (A1l < 1e-20) ? 0.0 : A1l * std::exp(alpha_mat[1][l] * dm);
    pk0[j] = base_const * k0 * inv_dj;
    pk1[j] = base_const * k1 * inv_dj;
  };

  for (int i = 0; i < n; ++i) {
    if (passign[i] == 0) {
      const int z = pinit[i];
      pid[i] = (z == 1) ? 1 : 0;
      fill_parent_const(i);
      continue;
    }

    double lam0 = mu_base_0 * pW0[i];
    double lam1 = mu_base_1 * pW1[i];
    const double ti = pt[i], xi = px[i], yi = py[i];

    if (do_trunc) {
      for (int j = i - 1; j >= 0; --j) {
        const double dt = ti - pt[j];
        if (dt > t_trunc) break;
        if (dt <= 0.0) continue;
        const double dx = xi - px[j];
        const double dy = yi - py[j];
        const double r2 = dx * dx + dy * dy;
        const double trig = std::exp(-p * std::log(1.0 + dt * inv_cc)
                                     - q * std::log(1.0 + r2 / d_par[j]));
        lam0 += pk0[j] * trig;
        lam1 += pk1[j] * trig;
      }
    } else {
      for (int j = 0; j < i; ++j) {
        const double dt = ti - pt[j];
        if (dt <= 0.0) continue;
        const double dx = xi - px[j];
        const double dy = yi - py[j];
        const double r2 = dx * dx + dy * dy;
        const double trig = std::exp(-p * std::log(1.0 + dt * inv_cc)
                                     - q * std::log(1.0 + r2 / d_par[j]));
        lam0 += pk0[j] * trig;
        lam1 += pk1[j] * trig;
      }
    }

    if (sample_bernoulli) {
      const double den = lam0 + lam1;
      const double p1 = (den > 0.0) ? (lam1 / den) : 0.5;
      pid[i] = (R::unif_rand() < p1) ? 1 : 0;
    } else {
      pid[i] = (lam1 > lam0) ? 1 : 0;
    }
    fill_parent_const(i);
  }

  return pid;
}
