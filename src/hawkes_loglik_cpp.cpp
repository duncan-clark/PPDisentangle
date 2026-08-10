#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// Both likelihoods below evaluate each pair's temporal x spatial kernel
// product as a single fused exponential (exp of summed log-terms) instead of
// separate pow()/exp() factors. This is mathematically identical to the
// original formulation (results agree to floating-point rounding) and
// roughly halves the transcendental cost per pair. The exponential spatial
// kernel's hard cutoff (r2 * alpha > 20 contributes nothing) is preserved
// exactly.

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

  const int n = t.size();

  const double pi = 3.14159265358979323846;
  const bool do_trunc = (t_trunc > 0.0);
  const bool power_law = (kernel_type == 1);
  bool spatial_power_law = (spatial_kernel_type == 1);
  if (cc <= 0.0) cc = 1.0;
  if (p <= 1.0) p = 1.000001;
  if (spatial_q <= 1.5) spatial_q = 1.500001;
  if (spatial_d <= 0.0) spatial_power_law = false;

  const double mu_base = mu / areaS;

  auto temporal_cdf = [&](double dt) {
    if (dt <= 0.0) return 0.0;
    if (power_law) {
      return 1.0 - std::pow(1.0 + dt / cc, 1.0 - p);
    }
    return 1.0 - std::exp(-beta * dt);
  };

  // Triggering constant, renormalized for temporal truncation so K remains
  // the branching ratio over the retained temporal support.
  double temporal_norm = do_trunc ? temporal_cdf(t_trunc) : 1.0;
  if (temporal_norm < 1e-15) temporal_norm = 1e-15;
  const double spatial_const = alpha / pi;
  const double spatial_power_const = (spatial_q - 1.0) / (pi * spatial_d);
  double dt_cut = power_law ? R_PosInf : 20.0 / beta;
  if (do_trunc && t_trunc < dt_cut) dt_cut = t_trunc;

  const double inv_cc = 1.0 / cc;
  const double inv_sd = spatial_power_law ? (1.0 / spatial_d) : 0.0;

  // Constant prefactor of every pair contribution:
  //   K * g_pref * f_pref / temporal_norm, with the dt/r2 dependence fused
  //   into a single exp() inside pair_term.
  const double temp_pref = power_law ? ((p - 1.0) * inv_cc) : beta;
  const double spat_pref = spatial_power_law ? spatial_power_const : spatial_const;
  const double pair_pref = K * temp_pref * spat_pref / temporal_norm;

  auto pair_term = [&](double dt, double r2) -> double {
    double expo = 0.0;
    if (power_law) expo -= p * std::log(1.0 + dt * inv_cc);
    else           expo -= beta * dt;
    if (spatial_power_law) {
      expo -= spatial_q * std::log(1.0 + r2 * inv_sd);
    } else {
      if (r2 * alpha > 20.0) return 0.0;
      expo -= alpha * r2;
    }
    return pair_pref * std::exp(expo);
  };

  const double* pt = t.begin();
  const double* px = x.begin();
  const double* py = y.begin();

  double loglik = 0.0;
  for (int i = 0; i < n; ++i) {
    double lambda_i = mu_base * W_val[i];
    const double ti = pt[i], xi = px[i], yi = py[i];

    for (int j = i - 1; j >= 0; --j) {
      const double dt = ti - pt[j];
      if (dt > dt_cut) break;
      if (dt <= 0.0) continue;  // temporal density is 0 at dt <= 0
      const double dx = xi - px[j];
      const double dy = yi - py[j];
      const double term = pair_term(dt, dx * dx + dy * dy);
      if (term <= 0.0) continue;
      lambda_i += term;
    }

    if (lambda_i <= 1e-15) lambda_i = 1e-15;
    loglik += std::log(lambda_i);
  }

  // Compensator integral. With truncation the temporal integral from event i
  // saturates: min(T - t_i, t_trunc) replaces (T - t_i).
  double triggering_integral = 0;
  for (int i = 0; i < n; ++i) {
    double horizon = t_max - pt[i];
    if (do_trunc && horizon > t_trunc) horizon = t_trunc;
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
  const int n_post = post_t.size();
  const int n_parent = parent_t.size();
  if (n_post == 0) return -1e15;

  const double pi = 3.14159265358979323846;
  const bool do_trunc = (t_trunc > 0.0);
  const bool power_law = (kernel_type == 1);
  bool spatial_power_law = (spatial_kernel_type == 1);
  if (cc <= 0.0) cc = 1.0;
  if (p <= 1.0) p = 1.000001;
  if (spatial_q <= 1.5) spatial_q = 1.500001;
  if (spatial_d <= 0.0) spatial_power_law = false;
  const double dt_window = t_end - t_start;
  if (dt_window <= 0.0 || areaS <= 0.0) return -1e15;

  const double mu_base = mu / areaS;
  auto temporal_cdf = [&](double dt) {
    if (dt <= 0.0) return 0.0;
    if (power_law) {
      return 1.0 - std::pow(1.0 + dt / cc, 1.0 - p);
    }
    return 1.0 - std::exp(-beta * dt);
  };

  double temporal_norm = do_trunc ? temporal_cdf(t_trunc) : 1.0;
  if (temporal_norm < 1e-15) temporal_norm = 1e-15;
  const double spatial_const = alpha / pi;
  const double spatial_power_const = (spatial_q - 1.0) / (pi * spatial_d);
  double dt_cut = power_law ? R_PosInf : 20.0 / beta;
  if (do_trunc && t_trunc < dt_cut) dt_cut = t_trunc;

  const double inv_cc = 1.0 / cc;
  const double inv_sd = spatial_power_law ? (1.0 / spatial_d) : 0.0;
  const double temp_pref = power_law ? ((p - 1.0) * inv_cc) : beta;
  const double spat_pref = spatial_power_law ? spatial_power_const : spatial_const;
  const double pair_pref = K * temp_pref * spat_pref / temporal_norm;

  auto pair_term = [&](double dt, double r2) -> double {
    double expo = 0.0;
    if (power_law) expo -= p * std::log(1.0 + dt * inv_cc);
    else           expo -= beta * dt;
    if (spatial_power_law) {
      expo -= spatial_q * std::log(1.0 + r2 * inv_sd);
    } else {
      if (r2 * alpha > 20.0) return 0.0;
      expo -= alpha * r2;
    }
    return pair_pref * std::exp(expo);
  };

  const double* ppt = parent_t.begin();
  const double* ppx = parent_x.begin();
  const double* ppy = parent_y.begin();

  bool parent_sorted = true;
  for (int j = 1; j < n_parent; ++j) {
    if (ppt[j] < ppt[j - 1]) {
      parent_sorted = false;
      break;
    }
  }

  double loglik = 0.0;
  for (int i = 0; i < n_post; ++i) {
    double lambda_i = mu_base * W_val[i];
    const double ti = post_t[i], xi = post_x[i], yi = post_y[i];
    if (parent_sorted) {
      const double* first_ge = std::lower_bound(ppt, ppt + n_parent, ti);
      int j_start = static_cast<int>(first_ge - ppt) - 1;
      for (int j = j_start; j >= 0; --j) {
        const double dt = ti - ppt[j];
        if (dt > dt_cut) break;
        const double dx = xi - ppx[j];
        const double dy = yi - ppy[j];
        const double term = pair_term(dt, dx * dx + dy * dy);
        if (term <= 0.0) continue;
        lambda_i += term;
      }
    } else {
      for (int j = n_parent - 1; j >= 0; --j) {
        const double dt = ti - ppt[j];
        if (dt <= 0.0) continue;
        if (dt > dt_cut) break;
        const double dx = xi - ppx[j];
        const double dy = yi - ppy[j];
        const double term = pair_term(dt, dx * dx + dy * dy);
        if (term <= 0.0) continue;
        lambda_i += term;
      }
    }
    if (lambda_i <= 1e-15) lambda_i = 1e-15;
    loglik += std::log(lambda_i);
  }

  double triggering_integral = 0.0;
  for (int j = 0; j < n_parent; ++j) {
    const double p_t = ppt[j];
    double start_h = t_start;
    if (p_t > start_h) start_h = p_t;
    double end_h = t_end;
    if (do_trunc) {
      const double trunc_end = p_t + t_trunc;
      if (trunc_end < end_h) end_h = trunc_end;
    }
    if (end_h <= start_h) continue;
    triggering_integral += K * (temporal_cdf(end_h - p_t) - temporal_cdf(start_h - p_t)) / temporal_norm;
  }

  loglik -= (mu * adjust_factor * dt_window + triggering_integral);
  if (NumericVector::is_na(loglik) || std::isinf(loglik)) return -1e15;
  return loglik;
}
