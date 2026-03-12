#include <Rcpp.h>
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
    double t_trunc = -1.0) {

  int n = t.size();
  if (n == 0) return -1e15;

  const double pi_val = 3.14159265358979323846;
  bool do_trunc = (t_trunc > 0.0);

  if (areaS_0 <= 0.0) areaS_0 = 1.0;
  if (areaS_1 <= 0.0) areaS_1 = 1.0;

  double mu_base_0 = mu_0 / areaS_0;
  double mu_base_1 = mu_1 / areaS_1;

  double temporal_norm = do_trunc ?
    (1.0 - std::pow(1.0 + t_trunc / cc, -(p - 1.0))) : 1.0;
  if (temporal_norm < 1e-15) temporal_norm = 1e-15;

  double base_const = (p - 1.0) * (q - 1.0) / (pi_val * cc * temporal_norm);

  // Kernel matrix: A_mat[child][parent], alpha_mat[child][parent]
  double A_mat[2][2], alpha_mat[2][2];
  A_mat[0][0] = A_00;   alpha_mat[0][0] = alpha_m_00;
  A_mat[0][1] = A_01;   alpha_mat[0][1] = alpha_m_01;
  A_mat[1][0] = A_10;   alpha_mat[1][0] = alpha_m_10;
  A_mat[1][1] = A_11;   alpha_mat[1][1] = alpha_m_11;

  double loglik = 0.0;

  // --- Sum of log-intensities ---
  for (int i = 0; i < n; ++i) {
    int k = process_id[i];

    double lambda_i = (k == 0) ? mu_base_0 * W_val_0[i]
                                : mu_base_1 * W_val_1[i];

    for (int j = 0; j < i; ++j) {
      double dt = t[i] - t[j];
      if (dt <= 0.0) continue;
      if (do_trunc && dt > t_trunc) continue;

      int l = process_id[j];
      double A_kl = A_mat[k][l];
      if (A_kl < 1e-20) continue;
      double alpha_kl = alpha_mat[k][l];

      double dm_j = mag[j] - m0;
      double kappa_j = A_kl * std::exp(alpha_kl * dm_j);
      double d_j = D * std::exp(gamma_par * dm_j);

      double temporal = std::pow(1.0 + dt / cc, -p);

      double dx = x[i] - x[j];
      double dy = y[i] - y[j];
      double r2 = dx * dx + dy * dy;
      double spatial = std::pow(1.0 + r2 / d_j, -q) / d_j;

      lambda_i += base_const * kappa_j * temporal * spatial;
    }

    if (lambda_i <= 1e-15) lambda_i = 1e-15;
    loglik += std::log(lambda_i);
  }

  // --- Compensator ---
  // Each parent j in process l contributes to both child process compensators:
  //   child 0: A_{0l} * exp(alpha_{0l} * dm_j) * G(h_j)
  //   child 1: A_{1l} * exp(alpha_{1l} * dm_j) * G(h_j)
  double comp_trig = 0.0;
  for (int j = 0; j < n; ++j) {
    double dm_j = mag[j] - m0;
    double horizon = t_max - t[j];
    if (do_trunc && horizon > t_trunc) horizon = t_trunc;
    if (horizon <= 0.0) continue;
    double G_h = 1.0 - std::pow(1.0 + horizon / cc, -(p - 1.0));

    int l = process_id[j];
    for (int k = 0; k < 2; ++k) {
      double A_kl = A_mat[k][l];
      if (A_kl < 1e-20) continue;
      comp_trig += A_kl * std::exp(alpha_mat[k][l] * dm_j) * G_h;
    }
  }

  loglik -= (mu_0 * t_max + mu_1 * t_max + comp_trig);

  if (NumericVector::is_na(loglik) || std::isinf(loglik)) return -1e15;

  return loglik;
}
