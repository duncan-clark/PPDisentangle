#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

//' ETAS log-likelihood for an inhomogeneous spatio-temporal point process
//'
//' Computes the log-likelihood of an Epidemic-Type Aftershock Sequence (ETAS)
//' model on a marked spatio-temporal point pattern.
//'
//' **Conditional intensity:**
//' \deqn{\lambda(t,x,y) = \frac{\mu}{|S|} W(x,y) +
//'   \sum_{t_j < t} \kappa(m_j)\, g(t - t_j)\, f(x - x_j, y - y_j \mid m_j)}
//'
//' **Productivity:**
//' \deqn{\kappa(m) = A \exp\bigl(\alpha_m (m - m_0)\bigr)}
//'
//' **Omori-Utsu temporal kernel (normalised density on \eqn{[0,\infty)}
//' when \eqn{p > 1}):**
//' \deqn{g(\Delta t) = \frac{p - 1}{c}\Bigl(1 + \frac{\Delta t}{c}\Bigr)^{-p}}
//'
//' **Isotropic power-law spatial kernel (Zhuang et al. 2002):**
//' \deqn{f(x,y \mid m) = \frac{q - 1}{\pi\, d(m)}
//'   \Bigl(1 + \frac{r^2}{d(m)}\Bigr)^{-q}, \quad r^2 = x^2 + y^2}
//' with \eqn{d(m) = D \exp\bigl(\gamma (m - m_0)\bigr)}.
//'
//' **Compensator (integral of intensity):**
//' \deqn{\int_0^T \lambda\,\mathrm{d}t \approx \mu T +
//'   \sum_i \kappa(m_i)\bigl[1 - \bigl(1 + h_i / c\bigr)^{-(p-1)}\bigr]}
//' where \eqn{h_i = \min(T - t_i,\; t_{\mathrm{trunc}})} when temporal
//' truncation is active and \eqn{h_i = T - t_i} otherwise.  The spatial
//' kernel integrates to 1 over \eqn{\mathbb{R}^2} (infinite-plane
//' approximation, consistent with the Hawkes implementation).
//'
//' When truncation is active the temporal density is renormalised so that
//' it integrates to 1 over \eqn{[0, t_{\mathrm{trunc}}]}:
//' \deqn{G(t_{\mathrm{trunc}}) = 1 - (1 + t_{\mathrm{trunc}}/c)^{-(p-1)}}
//'
//' The pairwise product of power-law kernels is evaluated as a single fused
//' exponential, exp(-p*log(1+dt/c) - q*log(1+r2/d)), and all per-parent
//' constants are hoisted out of the inner pair loop. This is mathematically
//' identical to evaluating the two pow() factors separately (results agree
//' to floating-point rounding).
//'
//' @param t      Numeric vector of event times (sorted ascending, shifted
//'               so the window starts at 0).
//' @param x      Numeric vector of x-coordinates.
//' @param y      Numeric vector of y-coordinates.
//' @param mag    Numeric vector of event magnitudes.
//' @param W_val  Numeric vector of inhomogeneous background weights (set to
//'               1 for homogeneous, 0 for points in a zero-background region).
//' @param mu     Background rate (events per unit time over the full window).
//' @param A      Productivity scaling constant.
//' @param alpha_m  Magnitude efficiency parameter.
//' @param cc     Omori-Utsu temporal offset (\eqn{c > 0}).  Named \code{cc}
//'               to avoid collision with R's \code{c()} function.
//' @param p      Omori-Utsu temporal exponent (\eqn{p > 1}).
//' @param D      Spatial spread base parameter (\eqn{D > 0}).
//' @param gamma_par  Magnitude-dependent spatial scaling (\eqn{\gamma \ge 0}).
//' @param q      Spatial power-law exponent (\eqn{q > 1}).
//' @param m0     Reference (cutoff) magnitude.
//' @param areaS  Area of the spatial observation window \eqn{|S|}.
//' @param t_max  Length of the temporal observation window \eqn{T}.
//' @param t_trunc  Temporal truncation horizon.  Set to a negative value
//'                 (default \code{-1}) to disable truncation.
//' @return Scalar log-likelihood value.
// [[Rcpp::export]]
double etas_loglik_inhom_cpp(NumericVector t,
                              NumericVector x,
                              NumericVector y,
                              NumericVector mag,
                              NumericVector W_val,
                              double mu, double A, double alpha_m,
                              double cc, double p, double D,
                              double gamma_par, double q, double m0,
                              double areaS, double t_max,
                              double t_trunc = -1.0) {

  const int n = t.size();
  double loglik = 0.0;

  const double pi_val = 3.14159265358979323846;
  const bool do_trunc = (t_trunc > 0.0);

  const double mu_base = mu / areaS;

  // Temporal truncation normalisation factor.
  // G(t_trunc) = 1 - (1 + t_trunc / c)^{-(p-1)}
  // When not truncating, temporal_norm = 1 (the Omori-Utsu density already
  // integrates to 1 over [0, inf) for p > 1).
  double temporal_norm = do_trunc ?
    (1.0 - std::pow(1.0 + t_trunc / cc, -(p - 1.0))) : 1.0;
  if (temporal_norm < 1e-15) temporal_norm = 1e-15;

  // Pre-factor that is constant across all (i, j) pairs:
  //   A * (p-1) * (q-1) / (pi * c * temporal_norm)
  // The magnitude-dependent terms kappa(m_j)/d(m_j) are folded into the
  // per-parent factor pk[j] below.
  const double base_const = A * (p - 1.0) * (q - 1.0) / (pi_val * cc * temporal_norm);
  const double inv_cc = 1.0 / cc;

  const double* pt = t.begin();
  const double* px = x.begin();
  const double* py = y.begin();
  const double* pmag = mag.begin();
  const double* pW = W_val.begin();

  std::vector<double> d_par(n), pk(n), comp_kappa(n);
  for (int j = 0; j < n; ++j) {
    const double dm = pmag[j] - m0;
    const double kappa = std::exp(alpha_m * dm);
    const double dj = D * std::exp(gamma_par * dm);
    d_par[j] = dj;
    pk[j] = base_const * kappa / dj;
    comp_kappa[j] = A * kappa;
  }

  // --- Sum of log-intensities ---
  for (int i = 0; i < n; ++i) {
    double lambda_i = mu_base * pW[i];
    const double ti = pt[i], xi = px[i], yi = py[i];

    if (do_trunc) {
      for (int j = i - 1; j >= 0; --j) {
        const double dt = ti - pt[j];
        if (dt > t_trunc) break;
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

  // --- Compensator (integral of intensity over the observation domain) ---
  // Temporal CDF: G(h) = 1 - (1 + h/c)^{-(p-1)}
  // Contribution from parent i: kappa(m_i) * G(horizon_i) / G(t_trunc).
  // The division by temporal_norm matches the truncation-renormalized kernel
  // used in the intensity (A = expected offspring within t_trunc); without it
  // the compensator undercharges productivity by the factor G(t_trunc).
  double triggering_integral = 0.0;
  for (int i = 0; i < n; ++i) {
    double horizon = t_max - pt[i];
    if (do_trunc && horizon > t_trunc) horizon = t_trunc;
    if (horizon <= 0.0) continue;
    triggering_integral += comp_kappa[i] *
      (1.0 - std::pow(1.0 + horizon / cc, -(p - 1.0))) / temporal_norm;
  }

  loglik -= (mu * t_max + triggering_integral);

  if (NumericVector::is_na(loglik) || std::isinf(loglik)) return -1e15;

  return loglik;
}

//' Conditional ETAS log-likelihood with observed parent history
//'
//' Evaluates the ETAS likelihood only for events in the post-treatment
//' observation window while allowing earlier events to contribute as parents.
//' The triggering compensator is integrated only over the observation window,
//' including the remaining offspring mass of pre-window parents.
//'
//' @param post_t,post_x,post_y Coordinates and times of observed events.
//' @param W_val Background weights for the observed events.
//' @param parent_t,parent_x,parent_y,parent_mag Coordinates, times, and
//'   magnitudes of all possible parents (history plus observed events), sorted
//'   by time.
//' @param t_start,t_end Bounds of the observation window.
//' @inheritParams etas_loglik_inhom_cpp
//' @return Scalar conditional log-likelihood value.
// [[Rcpp::export]]
double etas_loglik_inhom_filtration_cpp(NumericVector post_t,
                                         NumericVector post_x,
                                         NumericVector post_y,
                                         NumericVector W_val,
                                         NumericVector parent_t,
                                         NumericVector parent_x,
                                         NumericVector parent_y,
                                         NumericVector parent_mag,
                                         double mu, double A, double alpha_m,
                                         double cc, double p, double D,
                                         double gamma_par, double q, double m0,
                                         double areaS,
                                         double t_start, double t_end,
                                         double t_trunc = -1.0) {
  const int n_post = post_t.size();
  const int n_parent = parent_t.size();
  if (n_post == 0 || n_parent == 0 || areaS <= 0.0 || t_end <= t_start) {
    return -1e15;
  }

  const double pi_val = 3.14159265358979323846;
  const bool do_trunc = (t_trunc > 0.0);
  const double mu_base = mu / areaS;
  const double inv_cc = 1.0 / cc;

  auto temporal_cdf = [&](double dt) {
    if (dt <= 0.0) return 0.0;
    return 1.0 - std::pow(1.0 + dt * inv_cc, -(p - 1.0));
  };

  double temporal_norm = do_trunc ? temporal_cdf(t_trunc) : 1.0;
  if (temporal_norm < 1e-15) temporal_norm = 1e-15;

  const double base_const =
    A * (p - 1.0) * (q - 1.0) / (pi_val * cc * temporal_norm);

  const double* ppt = parent_t.begin();
  const double* ppx = parent_x.begin();
  const double* ppy = parent_y.begin();
  const double* ppmag = parent_mag.begin();

  std::vector<double> d_parent(n_parent), pair_parent(n_parent),
    compensator_parent(n_parent);
  for (int j = 0; j < n_parent; ++j) {
    const double dm = ppmag[j] - m0;
    const double productivity = std::exp(alpha_m * dm);
    const double spatial_scale = D * std::exp(gamma_par * dm);
    d_parent[j] = spatial_scale;
    pair_parent[j] = base_const * productivity / spatial_scale;
    compensator_parent[j] = A * productivity;
  }

  double loglik = 0.0;
  for (int i = 0; i < n_post; ++i) {
    double lambda_i = mu_base * W_val[i];
    const double ti = post_t[i];
    const double xi = post_x[i];
    const double yi = post_y[i];

    const double* first_ge = std::lower_bound(ppt, ppt + n_parent, ti);
    int j_start = static_cast<int>(first_ge - ppt) - 1;
    for (int j = j_start; j >= 0; --j) {
      const double dt = ti - ppt[j];
      if (do_trunc && dt > t_trunc) break;
      if (dt <= 0.0) continue;
      const double dx = xi - ppx[j];
      const double dy = yi - ppy[j];
      const double r2 = dx * dx + dy * dy;
      lambda_i += pair_parent[j] *
        std::exp(-p * std::log(1.0 + dt * inv_cc)
                 - q * std::log(1.0 + r2 / d_parent[j]));
    }

    if (lambda_i <= 1e-15) lambda_i = 1e-15;
    loglik += std::log(lambda_i);
  }

  double triggering_integral = 0.0;
  for (int j = 0; j < n_parent; ++j) {
    const double parent_time = ppt[j];
    double integration_start = std::max(t_start, parent_time);
    double integration_end = t_end;
    if (do_trunc) {
      integration_end = std::min(integration_end, parent_time + t_trunc);
    }
    if (integration_end <= integration_start) continue;
    triggering_integral += compensator_parent[j] *
      (temporal_cdf(integration_end - parent_time) -
       temporal_cdf(integration_start - parent_time)) / temporal_norm;
  }

  loglik -= mu * (t_end - t_start) + triggering_integral;
  if (NumericVector::is_na(loglik) || std::isinf(loglik)) return -1e15;
  return loglik;
}
