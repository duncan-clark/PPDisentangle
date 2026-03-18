#include <Rcpp.h>
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

  int n = t.size();
  double loglik = 0.0;

  const double pi_val = 3.14159265358979323846;
  bool do_trunc = (t_trunc > 0.0);

  double mu_base = mu / areaS;

  // Temporal truncation normalisation factor.
  // G(t_trunc) = 1 - (1 + t_trunc / c)^{-(p-1)}
  // When not truncating, temporal_norm = 1 (the Omori-Utsu density already
  // integrates to 1 over [0, inf) for p > 1).
  double temporal_norm = do_trunc ?
    (1.0 - std::pow(1.0 + t_trunc / cc, -(p - 1.0))) : 1.0;
  if (temporal_norm < 1e-15) temporal_norm = 1e-15;

  // Pre-factor that is constant across all (i, j) pairs:
  //   A * (p-1) * (q-1) / (pi * c * temporal_norm)
  // The magnitude-dependent terms kappa(m_j)/d(m_j) are applied per parent.
  double base_const = A * (p - 1.0) * (q - 1.0) / (pi_val * cc * temporal_norm);
  NumericVector dm(n), kappa_factor(n), d_parent(n), inv_d_parent(n), comp_kappa(n);
  for (int j = 0; j < n; ++j) {
    dm[j] = mag[j] - m0;
    kappa_factor[j] = std::exp(alpha_m * dm[j]);
    d_parent[j] = D * std::exp(gamma_par * dm[j]);
    inv_d_parent[j] = 1.0 / d_parent[j];
    comp_kappa[j] = A * kappa_factor[j];
  }

  // --- Sum of log-intensities ---
  for (int i = 0; i < n; ++i) {

    double lambda_i = mu_base * W_val[i];

    for (int j = 0; j < i; ++j) {
      double dt = t[i] - t[j];
      if (dt <= 0.0) continue;
      if (do_trunc && dt > t_trunc) continue;

      // Magnitude-dependent quantities for parent j
      double kappa_j = kappa_factor[j];
      double d_j = d_parent[j];

      // Omori-Utsu temporal factor: (1 + dt/c)^{-p}
      double temporal = std::pow(1.0 + dt / cc, -p);

      // Power-law spatial factor: (1 + r^2/d_j)^{-q} / d_j
      double dx = x[i] - x[j];
      double dy = y[i] - y[j];
      double r2 = dx * dx + dy * dy;
      double spatial = std::pow(1.0 + r2 / d_j, -q) * inv_d_parent[j];

      lambda_i += base_const * kappa_j * temporal * spatial;
    }

    if (lambda_i <= 1e-15) lambda_i = 1e-15;
    loglik += std::log(lambda_i);
  }

  // --- Compensator (integral of intensity over the observation domain) ---
  // Temporal CDF: G(h) = 1 - (1 + h/c)^{-(p-1)}
  // Contribution from parent i: kappa(m_i) * G(horizon_i)
  double triggering_integral = 0.0;
  for (int i = 0; i < n; ++i) {
    double kappa_i = comp_kappa[i];
    double horizon = t_max - t[i];
    if (do_trunc && horizon > t_trunc) horizon = t_trunc;
    if (horizon <= 0.0) continue;
    triggering_integral += kappa_i * (1.0 - std::pow(1.0 + horizon / cc, -(p - 1.0)));
  }

  loglik -= (mu * t_max + triggering_integral);

  if (NumericVector::is_na(loglik) || std::isinf(loglik)) return -1e15;

  return loglik;
}
