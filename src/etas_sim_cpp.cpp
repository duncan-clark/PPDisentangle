#include <Rcpp.h>
using namespace Rcpp;

//' Simulate ETAS offspring via BFS branching
//'
//' Starting from a set of parent events, generates all offspring (and their
//' descendants) using the ETAS triggering kernel with Omori-Utsu temporal
//' delays and isotropic power-law spatial displacements.
//'
//' **Temporal delays** are sampled by inverse-CDF from the Omori-Utsu
//' distribution:
//' \deqn{\Delta t = c\bigl[(1 - u)^{-1/(p-1)} - 1\bigr], \quad u \sim
//'   \mathrm{Uniform}(0, G_{\max})}
//' where \eqn{G_{\max} = 1} (no truncation) or
//' \eqn{G_{\max} = 1 - (1 + t_{\mathrm{trunc}}/c)^{-(p-1)}} (truncated).
//'
//' **Spatial displacements** are sampled by inverse-CDF from the isotropic
//' power-law kernel (magnitude-dependent spread):
//' \deqn{r^2 = d(m)\bigl[(1 - u)^{-1/(q-1)} - 1\bigr], \quad
//'   d(m) = D\exp\bigl(\gamma(m - m_0)\bigr)}
//' with direction uniform on \eqn{[0, 2\pi)}.
//'
//' **Offspring magnitudes** are drawn either from a Gutenberg-Richter
//' distribution \eqn{m = m_0 + \mathrm{Exp}(\beta_{GR})} (when
//' \code{beta_gr > 0}) or by resampling uniformly from the supplied
//' \code{mag_pool} vector.
//'
//' @param parent_x,parent_y,parent_t,parent_mag  Numeric vectors of parent
//'   coordinates and magnitudes.
//' @param A         Productivity scaling.
//' @param alpha_m   Magnitude efficiency.
//' @param cc        Omori-Utsu temporal offset (\eqn{c > 0}).
//' @param p         Omori-Utsu temporal exponent (\eqn{p > 1}).
//' @param D         Spatial spread base parameter.
//' @param gamma_par Magnitude-spatial scaling (\eqn{\gamma}).
//' @param q         Spatial power-law exponent (\eqn{q > 1}).
//' @param m0        Reference magnitude.
//' @param beta_gr   Gutenberg-Richter \eqn{\beta = b \ln 10}.  Set to 0 or
//'                  negative to use empirical resampling from \code{mag_pool}.
//' @param t_min,t_max  Temporal observation window.
//' @param x_min,x_max,y_min,y_max  Spatial bounding box.
//' @param t_trunc   Temporal truncation (negative to disable).
//' @param mag_pool  Numeric vector of observed magnitudes for resampling.
//'                  Ignored when \code{beta_gr > 0}.
//' @return A \code{DataFrame} with columns \code{x}, \code{y}, \code{t},
//'   \code{mag} containing all offspring (excluding parents).
// [[Rcpp::export]]
DataFrame sim_etas_children_cpp(NumericVector parent_x,
                                NumericVector parent_y,
                                NumericVector parent_t,
                                NumericVector parent_mag,
                                double A, double alpha_m,
                                double cc, double p,
                                double D, double gamma_par, double q,
                                double m0, double beta_gr,
                                double t_min, double t_max,
                                double x_min, double x_max,
                                double y_min, double y_max,
                                double t_trunc = -1.0,
                                NumericVector mag_pool = NumericVector::create()) {

  bool do_trunc = (t_trunc > 0.0);
  bool use_gr = (beta_gr > 0.0);
  int pool_n = mag_pool.size();

  // CDF at t_trunc for truncated Omori inverse-CDF sampling
  // G(t_trunc) = 1 - (1 + t_trunc/c)^{-(p-1)}
  double cdf_max = do_trunc ?
    (1.0 - std::pow(1.0 + t_trunc / cc, -(p - 1.0))) : 1.0;

  double pm1_inv = 1.0 / (p - 1.0);  // 1/(p-1) for temporal inverse CDF
  double qm1_inv = 1.0 / (q - 1.0);  // 1/(q-1) for spatial inverse CDF

  const double two_pi = 2.0 * 3.14159265358979323846;

  // BFS queue: start with all parents
  std::vector<double> q_x = as<std::vector<double>>(parent_x);
  std::vector<double> q_y = as<std::vector<double>>(parent_y);
  std::vector<double> q_t = as<std::vector<double>>(parent_t);
  std::vector<double> q_m = as<std::vector<double>>(parent_mag);

  // Output: offspring only (not parents)
  std::vector<double> out_x, out_y, out_t, out_m;
  int estimated_n = parent_t.size() * 2;
  if (estimated_n < 100) estimated_n = 100;
  out_x.reserve(estimated_n);
  out_y.reserve(estimated_n);
  out_t.reserve(estimated_n);
  out_m.reserve(estimated_n);

  size_t head = 0;

  while (head < q_x.size()) {
    double px = q_x[head];
    double py = q_y[head];
    double pt = q_t[head];
    double pm = q_m[head];
    head++;

    // Productivity: kappa(m) = A * exp(alpha_m * (m - m0))
    double dm = pm - m0;
    double kappa = A * std::exp(alpha_m * dm);

    int n_kids = R::rpois(kappa);
    if (n_kids == 0) continue;

    // Magnitude-dependent spatial spread: d(m) = D * exp(gamma * (m - m0))
    double d_parent = D * std::exp(gamma_par * dm);

    for (int k = 0; k < n_kids; ++k) {

      // --- Temporal delay: inverse-CDF of Omori-Utsu ---
      // dt = c * ((1 - u)^{-1/(p-1)} - 1),  u ~ Uniform(0, cdf_max)
      double u_t = R::runif(0.0, 1.0) * cdf_max;
      double dt = cc * (std::pow(1.0 - u_t, -pm1_inv) - 1.0);
      double new_t = pt + dt;

      if (new_t > t_max || new_t < t_min) continue;

      // --- Spatial displacement: inverse-CDF of power-law ---
      // r^2 = d_parent * ((1-u)^{-1/(q-1)} - 1),  u ~ Uniform(0,1)
      double u_s = R::runif(0.0, 1.0);
      double r2 = d_parent * (std::pow(1.0 - u_s, -qm1_inv) - 1.0);
      double r = std::sqrt(r2);
      double angle = R::runif(0.0, two_pi);

      double new_x = px + r * std::cos(angle);
      double new_y = py + r * std::sin(angle);

      if (new_x < x_min || new_x > x_max ||
          new_y < y_min || new_y > y_max) continue;

      // --- Offspring magnitude ---
      double new_mag;
      if (use_gr) {
        new_mag = m0 + R::rexp(1.0 / beta_gr);
      } else if (pool_n > 0) {
        int idx = (int)(R::runif(0.0, (double)pool_n));
        if (idx >= pool_n) idx = pool_n - 1;
        new_mag = mag_pool[idx];
      } else {
        new_mag = m0;
      }

      // Accept offspring: add to BFS queue and output
      q_x.push_back(new_x);
      q_y.push_back(new_y);
      q_t.push_back(new_t);
      q_m.push_back(new_mag);

      out_x.push_back(new_x);
      out_y.push_back(new_y);
      out_t.push_back(new_t);
      out_m.push_back(new_mag);
    }
  }

  return List::create(
    Named("x") = out_x,
    Named("y") = out_y,
    Named("t") = out_t,
    Named("mag") = out_m
  );
}
