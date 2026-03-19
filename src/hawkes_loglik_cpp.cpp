#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double hawkes_loglik_inhom_cpp(NumericVector t,
                               NumericVector x,
                               NumericVector y,
                               NumericVector W_val,
                               double mu, double alpha, double beta, double K,
                               double areaS, double t_max,
                               double t_trunc = -1.0) {

  int n = t.size();
  double loglik = 0.0;

  double pi = 3.14159265358979323846;
  bool do_trunc = (t_trunc > 0.0);

  double mu_base = mu / areaS;

  // Triggering constant, renormalized for truncation.
  // Untruncated: integral of beta*exp(-beta*dt) from 0 to inf = 1.
  // Truncated to [0, t_trunc]: integral = 1 - exp(-beta*t_trunc).
  // We divide by this factor so the kernel still integrates to K over [0,t_trunc] x R^2.
  double temporal_norm = do_trunc ? (1.0 - std::exp(-beta * t_trunc)) : 1.0;
  if (temporal_norm < 1e-15) temporal_norm = 1e-15;
  double const_val = K * alpha * beta / (pi * temporal_norm);

  for(int i = 0; i < n; ++i) {

    double lambda_i = mu_base * W_val[i];

    for(int j = 0; j < i; ++j) {
      double dt = t[i] - t[j];

      if(do_trunc && dt > t_trunc) continue;

      if(dt * beta > 20.0) continue;

      double dx = x[i] - x[j];
      double dy = y[i] - y[j];
      double r2 = dx*dx + dy*dy;

      if(r2 * alpha > 20.0) continue;

      double excitation = std::exp(-beta * dt - alpha * r2);
      lambda_i += const_val * excitation;
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
    triggering_integral += K * (1.0 - std::exp(-beta * horizon));
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
                                          double t_trunc = -1.0) {
  int n_post = post_t.size();
  int n_parent = parent_t.size();
  if (n_post == 0) return -1e15;

  const double pi = 3.14159265358979323846;
  bool do_trunc = (t_trunc > 0.0);
  double dt_window = t_end - t_start;
  if (dt_window <= 0.0 || areaS <= 0.0) return -1e15;

  double mu_base = mu / areaS;
  double temporal_norm = do_trunc ? (1.0 - std::exp(-beta * t_trunc)) : 1.0;
  if (temporal_norm < 1e-15) temporal_norm = 1e-15;
  double const_val = K * alpha * beta / (pi * temporal_norm);

  double loglik = 0.0;
  for (int i = 0; i < n_post; ++i) {
    double lambda_i = mu_base * W_val[i];
    double ti = post_t[i];
    for (int j = 0; j < n_parent; ++j) {
      double dt = ti - parent_t[j];
      if (dt <= 0.0) continue;
      if (do_trunc && dt > t_trunc) continue;
      if (dt * beta > 20.0) continue;
      double dx = post_x[i] - parent_x[j];
      double dy = post_y[i] - parent_y[j];
      double r2 = dx * dx + dy * dy;
      if (r2 * alpha > 20.0) continue;
      lambda_i += const_val * std::exp(-beta * dt - alpha * r2);
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
    triggering_integral += K * (std::exp(-beta * (start_h - p_t)) - std::exp(-beta * (end_h - p_t)));
  }

  loglik -= (mu * adjust_factor * dt_window + triggering_integral);
  if (NumericVector::is_na(loglik) || std::isinf(loglik)) return -1e15;
  return loglik;
}
