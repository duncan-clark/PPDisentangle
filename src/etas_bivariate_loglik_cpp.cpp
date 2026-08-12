#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "omori_kernel.h"
using namespace Rcpp;

//' Bivariate ETAS log-likelihood with cross-excitation
//'
//' Computes the joint log-likelihood of a two-component ETAS model where
//' events in each process can trigger offspring in both processes via a
//' 2x2 kernel matrix.
//'
//' The conditional intensity for process k in {0, 1} is:
//' \deqn{\lambda_k(t,x,y) = \frac{\mu_k}{|S_k|} W_k(x,y) +
//'   \sum_{l \in \{0,1\}} \sum_{t_j < t,\, \mathrm{proc}_j = l}
//'     A_{kl}\,e^{\alpha_{kl}(m_j - m_0)}\,g(\Delta t)\,f(\Delta x, \Delta y | m_j)}
//'
//' The four kernel components (child, parent) are:
//'   - (0,0): control self-excitation
//'   - (0,1): treated events exciting control intensity
//'   - (1,0): control events exciting treated intensity
//'   - (1,1): treated self-excitation
//'
//' Temporal and spatial kernels g, f use shared structural parameters
//' (c, p, D, gamma, q), identical to the univariate ETAS formulation.
//'
//' The pairwise product of power-law kernels is evaluated as a single fused
//' exponential, exp(-p*log(1+dt/c) - q*log(1+r2/d)), and all per-parent
//' constants are hoisted out of the inner pair loop. This is mathematically
//' identical to evaluating the two pow() factors separately (results agree
//' to floating-point rounding).
//'
//' @param t          Event times (sorted ascending, window starts at 0).
//' @param x,y        Spatial coordinates.
//' @param mag        Event magnitudes.
//' @param process_id Integer vector: 0 = control, 1 = treated.
//' @param W_val_0    Background weights for control process (0 in treated region).
//' @param W_val_1    Background weights for treated process (0 in control region).
//' @param mu_0,mu_1  Background rates for control and treated.
//' @param A_00,alpha_m_00  Control self-excitation productivity.
//' @param A_11,alpha_m_11  Treated self-excitation productivity.
//' @param A_01,alpha_m_01  Treated-to-control cross-excitation.
//' @param A_10,alpha_m_10  Control-to-treated cross-excitation.
//' @param cc,p,D,gamma_par,q  Shared structural parameters.
//' @param m0         Reference magnitude.
//' @param areaS_0    Active area for control background.
//' @param areaS_1    Active area for treated background.
//' @param t_max      Length of temporal observation window.
//' @param t_trunc    Temporal truncation (negative to disable).
//' @param bg_exposure_0,bg_exposure_1 Time span over which each background is
//'   switched on (compensator charge). Negative (default) means the full
//'   window \code{t_max}. Use these when a policy mask zeroes a background
//'   for part of the window (e.g. treated background only after treatment).
//' @return Scalar joint log-likelihood.
// [[Rcpp::export]]
double etas_bivariate_loglik_cpp(
    NumericVector t, NumericVector x, NumericVector y, NumericVector mag,
    IntegerVector process_id,
    NumericVector W_val_0, NumericVector W_val_1,
    double mu_0, double mu_1,
    double A_00, double alpha_m_00,
    double A_11, double alpha_m_11,
    double A_01, double alpha_m_01,
    double A_10, double alpha_m_10,
    double cc, double p, double D, double gamma_par, double q,
    double m0,
    double areaS_0, double areaS_1,
    double t_max,
    double t_trunc = -1.0,
    double bg_exposure_0 = -1.0,
    double bg_exposure_1 = -1.0) {

  const int n = t.size();
  if (n == 0) return -1e15;

  const double pi_val = 3.14159265358979323846;
  const bool do_trunc = (t_trunc > 0.0);

  if (areaS_0 <= 0.0) areaS_0 = 1.0;
  if (areaS_1 <= 0.0) areaS_1 = 1.0;

  const double mu_base_0 = mu_0 / areaS_0;
  const double mu_base_1 = mu_1 / areaS_1;

  const double omori_pref = omori_time_prefactor(p, cc, t_trunc, do_trunc);
  const double base_const = omori_pref * (q - 1.0) / pi_val;

  // Kernel matrix: A_mat[child][parent], alpha_mat[child][parent]
  double A_mat[2][2], alpha_mat[2][2];
  A_mat[0][0] = A_00;   alpha_mat[0][0] = alpha_m_00;
  A_mat[0][1] = A_01;   alpha_mat[0][1] = alpha_m_01;
  A_mat[1][0] = A_10;   alpha_mat[1][0] = alpha_m_10;
  A_mat[1][1] = A_11;   alpha_mat[1][1] = alpha_m_11;

  const double* pt = t.begin();
  const double* px = x.begin();
  const double* py = y.begin();
  const double* pmag = mag.begin();
  const int*    ppid = process_id.begin();
  const double* pW0 = W_val_0.begin();
  const double* pW1 = W_val_1.begin();

  // Per-parent precomputation. pk{k}[j] = base_const * kappa(k,j) / d(j)
  // so the inner pair loop multiplies by a single fused exponential.
  std::vector<double> d_par(n);
  std::vector<double> kap0(n), kap1(n), pk0(n), pk1(n);
  const double inv_cc = 1.0 / cc;
  for (int j = 0; j < n; ++j) {
    const double dm = pmag[j] - m0;
    const double dj = D * std::exp(gamma_par * dm);
    d_par[j] = dj;
    const double inv_dj = 1.0 / dj;
    const int l = ppid[j];
    const double A0l = A_mat[0][l];
    const double A1l = A_mat[1][l];
    const double k0 = (A0l < 1e-20) ? 0.0 : A0l * std::exp(alpha_mat[0][l] * dm);
    const double k1 = (A1l < 1e-20) ? 0.0 : A1l * std::exp(alpha_mat[1][l] * dm);
    kap0[j] = k0;
    kap1[j] = k1;
    pk0[j] = base_const * k0 * inv_dj;
    pk1[j] = base_const * k1 * inv_dj;
  }

  double loglik = 0.0;

  // --- Sum of log-intensities ---
  for (int i = 0; i < n; ++i) {
    const int k = ppid[i];
    double lambda_i = (k == 0) ? mu_base_0 * pW0[i] : mu_base_1 * pW1[i];
    const double ti = pt[i], xi = px[i], yi = py[i];
    const double* kap = (k == 0) ? kap0.data() : kap1.data();
    const double* pk  = (k == 0) ? pk0.data()  : pk1.data();

    if (do_trunc) {
      for (int j = i - 1; j >= 0; --j) {
        const double dt = ti - pt[j];
        if (dt > t_trunc) break;
        if (kap[j] < 1e-20) continue;
        const double dx = xi - px[j];
        const double dy = yi - py[j];
        const double r2 = dx * dx + dy * dy;
        // (1+dt/c)^{-p} * (1+r2/d)^{-q} fused into one exp
        lambda_i += pk[j] * std::exp(-p * std::log(1.0 + dt * inv_cc)
                                     - q * std::log(1.0 + r2 / d_par[j]));
      }
    } else {
      for (int j = 0; j < i; ++j) {
        const double dt = ti - pt[j];
        if (dt <= 0.0) continue;
        if (kap[j] < 1e-20) continue;
        const double dx = xi - px[j];
        const double dy = yi - py[j];
        const double r2 = dx * dx + dy * dy;
        lambda_i += pk[j] * std::exp(-p * std::log(1.0 + dt * inv_cc)
                                     - q * std::log(1.0 + r2 / d_par[j]));
      }
    }

    if (lambda_i <= 1e-15) lambda_i = 1e-15;
    loglik += std::log(lambda_i);
  }

  // --- Compensator ---
  // Each parent j in process l contributes to both child process compensators:
  //   child 0: A_{0l} * exp(alpha_{0l} * dm_j) * G(h_j) / G(t_trunc)
  //   child 1: A_{1l} * exp(alpha_{1l} * dm_j) * G(h_j) / G(t_trunc)
  // The division by temporal_norm matches the truncation-renormalized kernel
  // used in the intensity (A = expected offspring within t_trunc); without it
  // the compensator undercharges productivity by the factor G(t_trunc).
  double comp_trig = 0.0;
  for (int j = 0; j < n; ++j) {
    double horizon = t_max - pt[j];
    if (do_trunc && horizon > t_trunc) horizon = t_trunc;
    if (horizon <= 0.0) continue;
    const double G_h = omori_time_cdf(p, cc, horizon, t_trunc, do_trunc);
    if (kap0[j] >= 1e-20) comp_trig += kap0[j] * G_h;
    if (kap1[j] >= 1e-20) comp_trig += kap1[j] * G_h;
  }

  // Background compensator: charge each mu_k only over the time span its
  // background is actually on (policy masks can shorten the exposure).
  const double expo_0 = (bg_exposure_0 >= 0.0) ? bg_exposure_0 : t_max;
  const double expo_1 = (bg_exposure_1 >= 0.0) ? bg_exposure_1 : t_max;
  loglik -= (mu_0 * expo_0 + mu_1 * expo_1 + comp_trig);

  if (NumericVector::is_na(loglik) || std::isinf(loglik)) return -1e15;

  return loglik;
}
