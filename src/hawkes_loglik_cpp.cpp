#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
double hawkes_loglik_inhom_cpp(NumericVector t,
                               NumericVector x,
                               NumericVector y,
                               NumericVector W_val,
                               double mu, double alpha, double beta, double K,
                               double areaS, double t_max,
                               double t_trunc = -1.0,
                               int kernel_type = 0,
                               double cc = 1.0,
                               double p = 2.0,
                               int spatial_kernel_type = 0,
                               double spatial_q = 2.0,
                               double spatial_d = -1.0) {

  int n = t.size();
  double loglik = 0.0;

  double pi = 3.14159265358979323846;
  bool do_trunc = (t_trunc > 0.0);
  bool power_law = (kernel_type == 1);
  bool spatial_power_law = (spatial_kernel_type == 1);
  if (cc <= 0.0) cc = 1.0;
  if (p <= 1.0) p = 1.000001;
  if (spatial_q <= 1.5) spatial_q = 1.500001;
  if (spatial_d <= 0.0) spatial_power_law = false;

  double mu_base = mu / areaS;

  auto temporal_cdf = [&](double dt) {
    if (dt <= 0.0) return 0.0;
    if (power_law) {
      return 1.0 - std::pow(1.0 + dt / cc, 1.0 - p);
    }
    return 1.0 - std::exp(-beta * dt);
  };
  auto temporal_density = [&](double dt) {
    if (dt <= 0.0) return 0.0;
    if (power_law) {
      return ((p - 1.0) / cc) * std::pow(1.0 + dt / cc, -p);
    }
    return beta * std::exp(-beta * dt);
  };

  // Triggering constant, renormalized for temporal truncation so K remains
  // the branching ratio over the retained temporal support.
  double temporal_norm = do_trunc ? temporal_cdf(t_trunc) : 1.0;
  if (temporal_norm < 1e-15) temporal_norm = 1e-15;
  double spatial_const = alpha / pi;
  double spatial_power_const = (spatial_q - 1.0) / (pi * spatial_d);
  double dt_cut = power_law ? R_PosInf : 20.0 / beta;
  if (do_trunc && t_trunc < dt_cut) dt_cut = t_trunc;
  auto spatial_density = [&](double r2) {
    if (spatial_power_law) {
      return spatial_power_const * std::pow(1.0 + r2 / spatial_d, -spatial_q);
    }
    if(r2 * alpha > 20.0) return 0.0;
    return spatial_const * std::exp(-alpha * r2);
  };

  for(int i = 0; i < n; ++i) {

    double lambda_i = mu_base * W_val[i];

    for(int j = i - 1; j >= 0; --j) {
      double dt = t[i] - t[j];

      if(dt > dt_cut) break;

      double dx = x[i] - x[j];
      double dy = y[i] - y[j];
      double r2 = dx*dx + dy*dy;

      double s_density = spatial_density(r2);
      if(s_density <= 0.0) continue;
      lambda_i += K * temporal_density(dt) * s_density / temporal_norm;
    }

    if(lambda_i <= 1e-15) lambda_i = 1e-15;

    loglik += std::log(lambda_i);
  }

  // Compensator integral. With truncation the temporal integral from event i
  // saturates: min(T - t_i, t_trunc) replaces (T - t_i).
  double triggering_integral = 0;
  for(int i = 0; i < n; ++i) {
    double horizon = t_max - t[i];
    if(do_trunc && horizon > t_trunc) horizon = t_trunc;
    triggering_integral += K * temporal_cdf(horizon) / temporal_norm;
  }

  loglik -= (mu * t_max + triggering_integral);

  if (NumericVector::is_na(loglik) || std::isinf(loglik)) return -1e15;

  return loglik;
}

// [[Rcpp::export]]
double hawkes_loglik_inhom_filtration_cpp(NumericVector post_t,
                                          NumericVector post_x,
                                          NumericVector post_y,
                                          NumericVector W_val,
                                          NumericVector parent_t,
                                          NumericVector parent_x,
                                          NumericVector parent_y,
                                          double mu, double alpha, double beta, double K,
                                          double areaS,
                                          double t_start, double t_end,
                                          double adjust_factor = 1.0,
                                          double t_trunc = -1.0,
                                          int kernel_type = 0,
                                          double cc = 1.0,
                                          double p = 2.0,
                                          int spatial_kernel_type = 0,
                                          double spatial_q = 2.0,
                                          double spatial_d = -1.0) {
  int n_post = post_t.size();
  int n_parent = parent_t.size();
  if (n_post == 0) return -1e15;

  const double pi = 3.14159265358979323846;
  bool do_trunc = (t_trunc > 0.0);
  bool power_law = (kernel_type == 1);
  bool spatial_power_law = (spatial_kernel_type == 1);
  if (cc <= 0.0) cc = 1.0;
  if (p <= 1.0) p = 1.000001;
  if (spatial_q <= 1.5) spatial_q = 1.500001;
  if (spatial_d <= 0.0) spatial_power_law = false;
  double dt_window = t_end - t_start;
  if (dt_window <= 0.0 || areaS <= 0.0) return -1e15;

  double mu_base = mu / areaS;
  auto temporal_cdf = [&](double dt) {
    if (dt <= 0.0) return 0.0;
    if (power_law) {
      return 1.0 - std::pow(1.0 + dt / cc, 1.0 - p);
    }
    return 1.0 - std::exp(-beta * dt);
  };
  auto temporal_density = [&](double dt) {
    if (dt <= 0.0) return 0.0;
    if (power_law) {
      return ((p - 1.0) / cc) * std::pow(1.0 + dt / cc, -p);
    }
    return beta * std::exp(-beta * dt);
  };

  double temporal_norm = do_trunc ? temporal_cdf(t_trunc) : 1.0;
  if (temporal_norm < 1e-15) temporal_norm = 1e-15;
  double spatial_const = alpha / pi;
  double spatial_power_const = (spatial_q - 1.0) / (pi * spatial_d);
  double dt_cut = power_law ? R_PosInf : 20.0 / beta;
  if (do_trunc && t_trunc < dt_cut) dt_cut = t_trunc;
  auto spatial_density = [&](double r2) {
    if (spatial_power_law) {
      return spatial_power_const * std::pow(1.0 + r2 / spatial_d, -spatial_q);
    }
    if (r2 * alpha > 20.0) return 0.0;
    return spatial_const * std::exp(-alpha * r2);
  };

  double loglik = 0.0;
  bool parent_sorted = true;
  for (int j = 1; j < n_parent; ++j) {
    if (parent_t[j] < parent_t[j - 1]) {
      parent_sorted = false;
      break;
    }
  }

  for (int i = 0; i < n_post; ++i) {
    double lambda_i = mu_base * W_val[i];
    double ti = post_t[i];
    if (parent_sorted) {
      NumericVector::iterator first_ge = std::lower_bound(parent_t.begin(), parent_t.end(), ti);
      int j_start = static_cast<int>(first_ge - parent_t.begin()) - 1;
      for (int j = j_start; j >= 0; --j) {
        double dt = ti - parent_t[j];
        if (dt > dt_cut) break;
        double dx = post_x[i] - parent_x[j];
        double dy = post_y[i] - parent_y[j];
        double r2 = dx * dx + dy * dy;
        double s_density = spatial_density(r2);
        if (s_density <= 0.0) continue;
        lambda_i += K * temporal_density(dt) * s_density / temporal_norm;
      }
    } else {
      for (int j = n_parent - 1; j >= 0; --j) {
        double dt = ti - parent_t[j];
        if (dt <= 0.0) continue;
        if (dt > dt_cut) break;
        double dx = post_x[i] - parent_x[j];
        double dy = post_y[i] - parent_y[j];
        double r2 = dx * dx + dy * dy;
        double s_density = spatial_density(r2);
        if (s_density <= 0.0) continue;
        lambda_i += K * temporal_density(dt) * s_density / temporal_norm;
      }
    }
    if (lambda_i <= 1e-15) lambda_i = 1e-15;
    loglik += std::log(lambda_i);
  }

  double triggering_integral = 0.0;
  for (int j = 0; j < n_parent; ++j) {
    double p_t = parent_t[j];
    double start_h = t_start;
    if (p_t > start_h) start_h = p_t;
    double end_h = t_end;
    if (do_trunc) {
      double trunc_end = p_t + t_trunc;
      if (trunc_end < end_h) end_h = trunc_end;
    }
    if (end_h <= start_h) continue;
    triggering_integral += K * (temporal_cdf(end_h - p_t) - temporal_cdf(start_h - p_t)) / temporal_norm;
  }

  loglik -= (mu * adjust_factor * dt_window + triggering_integral);
  if (NumericVector::is_na(loglik) || std::isinf(loglik)) return -1e15;
  return loglik;
}
