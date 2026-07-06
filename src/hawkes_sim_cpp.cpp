#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame sim_hawkes_children_cpp(NumericVector parent_x,
                                  NumericVector parent_y,
                                  NumericVector parent_t,
                                  double alpha,
                                  double beta,
                                  double K,
                                  double t_min,
                                  double t_max,
                                  double x_min, double x_max,
                                  double y_min, double y_max,
                                  double t_trunc = -1.0,
                                  int kernel_type = 0,
                                  double cc = 1.0,
                                  double p = 2.0) {

  bool do_trunc = (t_trunc > 0.0);
  bool power_law = (kernel_type == 1);
  if (cc <= 0.0) cc = 1.0;
  if (p <= 1.0) p = 1.000001;

  std::vector<double> out_x;
  std::vector<double> out_y;
  std::vector<double> out_t;

  int estimated_n = parent_t.size() * K * 2;
  if(estimated_n < 100) estimated_n = 100;
  out_x.reserve(estimated_n);
  out_y.reserve(estimated_n);
  out_t.reserve(estimated_n);

  std::vector<double> q_x = as<std::vector<double>>(parent_x);
  std::vector<double> q_y = as<std::vector<double>>(parent_y);
  std::vector<double> q_t = as<std::vector<double>>(parent_t);

  auto temporal_cdf = [&](double dt) {
    if (dt <= 0.0) return 0.0;
    if (power_law) {
      return 1.0 - std::pow(1.0 + dt / cc, 1.0 - p);
    }
    return 1.0 - std::exp(-beta * dt);
  };
  auto temporal_quantile = [&](double u) {
    if (u <= 0.0) return 0.0;
    if (u >= 1.0) u = 1.0 - 1e-15;
    if (power_law) {
      return cc * (std::pow(1.0 - u, -1.0 / (p - 1.0)) - 1.0);
    }
    return -std::log(1.0 - u) / beta;
  };

  // CDF normalisation for the triggering-time law. With truncation the
  // temporal density is renormalised onto [0, t_trunc].
  double cdf_norm = do_trunc ? temporal_cdf(t_trunc) : 1.0;
  if (cdf_norm < 1e-15) cdf_norm = 1e-15;

  size_t head = 0;

  while(head < q_x.size()) {
    double px = q_x[head];
    double py = q_y[head];
    double pt = q_t[head];
    head++;

    double lower_dt = t_min - pt;
    if (lower_dt < 0.0) lower_dt = 0.0;
    double upper_dt = t_max - pt;
    if (do_trunc && upper_dt > t_trunc) upper_dt = t_trunc;
    if (upper_dt <= lower_dt) continue;

    double cdf_lower = temporal_cdf(lower_dt);
    double cdf_upper = temporal_cdf(upper_dt);
    double observable_mass = (cdf_upper - cdf_lower) / cdf_norm;
    if (observable_mass <= 0.0) continue;

    int n_kids = R::rpois(K * observable_mass);
    if(n_kids == 0) continue;

    for(int k = 0; k < n_kids; ++k) {
      double u = R::runif(cdf_lower, cdf_upper);
      double dt = temporal_quantile(u);
      double new_t = pt + dt;

      double r2 = R::rexp(1.0/alpha);
      double dist = std::sqrt(r2);
      double angle = R::runif(0.0, 2.0 * 3.14159265358979323846);

      double new_x = px + dist * cos(angle);
      double new_y = py + dist * sin(angle);

      if(new_x >= x_min && new_x <= x_max &&
         new_y >= y_min && new_y <= y_max) {

        q_x.push_back(new_x);
        q_y.push_back(new_y);
        q_t.push_back(new_t);

        out_x.push_back(new_x);
        out_y.push_back(new_y);
        out_t.push_back(new_t);
      }
    }
  }

  return List::create(
    Named("x") = out_x,
    Named("y") = out_y,
    Named("t") = out_t
  );
}
