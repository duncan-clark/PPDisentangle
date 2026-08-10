// Batched Hawkes log-likelihood with filtration history over labellings.
// Events share (t,x,y); labellings differ only in process membership.
// Pairwise temporal x spatial kernel weights are computed once per (i,j)
// and reused across K labellings. Each entry matches
// hawkes_loglik_inhom_filtration_cpp for that labelling's post/parent sets
// (to floating-point rounding).
//
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <thread>
#include <algorithm>
using namespace Rcpp;

static const double PI_VAL = 3.14159265358979323846;

//' Batched Hawkes filtration log-likelihood over labellings
//'
//' @param t,x,y Shared event geometry, sorted ascending in time.
//' @param is_observed Integer 0/1 length-n: 1 = contributes a log-intensity
//'   term (post-treatment observation); 0 = parent-only (pre-history).
//' @param member Integer matrix n x K: 1 if event belongs to the process
//'   under that labelling (parents + observed members).
//' @param W_val Length-n background weights (used only for observed events).
//' @param mu,alpha,beta,K Hawkes parameters.
//' @param areaS Active spatial area.
//' @param t_start,t_end Observation window for the compensator.
//' @param adjust_factor Background compensator scale (usually 1).
//' @param t_trunc Temporal truncation (negative disables).
//' @param kernel_type 0 = exponential, 1 = power-law temporal.
//' @param cc,p Power-law temporal parameters.
//' @param spatial_kernel_type 0 = exponential, 1 = power-law spatial.
//' @param spatial_q,spatial_d Power-law spatial parameters.
//' @param n_threads Worker threads (1 = serial).
//' @return Length-K numeric vector of log-likelihoods.
//' @keywords internal
// [[Rcpp::export]]
NumericVector hawkes_loglik_inhom_filtration_batch_cpp(
    NumericVector t, NumericVector x, NumericVector y,
    IntegerVector is_observed,
    IntegerMatrix member,
    NumericVector W_val,
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
    double spatial_d = -1.0,
    int n_threads = 1) {

  const int n = t.size();
  const int Klab = member.ncol();
  NumericVector out(Klab);
  if (n == 0 || Klab < 1) {
    for (int kk = 0; kk < Klab; ++kk) out[kk] = -1e15;
    return out;
  }
  if (n_threads < 1) n_threads = 1;

  const bool do_trunc = (t_trunc > 0.0);
  const bool power_law = (kernel_type == 1);
  bool spatial_power_law = (spatial_kernel_type == 1);
  if (cc <= 0.0) cc = 1.0;
  if (p <= 1.0) p = 1.000001;
  if (spatial_q <= 1.5) spatial_q = 1.500001;
  if (spatial_d <= 0.0) spatial_power_law = false;
  const double dt_window = t_end - t_start;
  if (dt_window <= 0.0 || areaS <= 0.0) {
    for (int kk = 0; kk < Klab; ++kk) out[kk] = -1e15;
    return out;
  }

  auto temporal_cdf = [&](double dt) {
    if (dt <= 0.0) return 0.0;
    if (power_law) return 1.0 - std::pow(1.0 + dt / cc, 1.0 - p);
    return 1.0 - std::exp(-beta * dt);
  };
  double temporal_norm = do_trunc ? temporal_cdf(t_trunc) : 1.0;
  if (temporal_norm < 1e-15) temporal_norm = 1e-15;
  const double spatial_const = alpha / PI_VAL;
  const double spatial_power_const = (spatial_q - 1.0) / (PI_VAL * spatial_d);
  double dt_cut = power_law ? R_PosInf : 20.0 / beta;
  if (do_trunc && t_trunc < dt_cut) dt_cut = t_trunc;

  const double inv_cc = 1.0 / cc;
  const double inv_sd = spatial_power_law ? (1.0 / spatial_d) : 0.0;
  const double temp_pref = power_law ? ((p - 1.0) * inv_cc) : beta;
  const double spat_pref = spatial_power_law ? spatial_power_const : spatial_const;
  const double pair_pref = K * temp_pref * spat_pref / temporal_norm;
  const double mu_base = mu / areaS;
  const double bg_comp = mu * adjust_factor * dt_window;

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

  // Plain buffers for thread safety.
  std::vector<double> vt(t.begin(), t.end()), vx(x.begin(), x.end()),
      vy(y.begin(), y.end()), vW(W_val.begin(), W_val.end());
  std::vector<int> obs(is_observed.begin(), is_observed.end());
  std::vector<int> mem(static_cast<size_t>(n) * Klab);
  for (int kk = 0; kk < Klab; ++kk) {
    for (int i = 0; i < n; ++i) {
      mem[static_cast<size_t>(kk) * n + i] = member(i, kk);
    }
  }

  // Which labellings have at least one observed member?
  std::vector<char> has_obs(Klab, 0);
  for (int kk = 0; kk < Klab; ++kk) {
    const int* mk = mem.data() + static_cast<size_t>(kk) * n;
    for (int i = 0; i < n; ++i) {
      if (obs[i] && mk[i]) { has_obs[kk] = 1; break; }
    }
  }

  const int T = std::min<int>(n_threads, std::max(1, n));
  std::vector<std::vector<double>> partial(T, std::vector<double>(Klab, 0.0));

  auto worker = [&](int tid) {
    std::vector<double> wbuf;
    std::vector<int> jbuf;
    wbuf.reserve(1024); jbuf.reserve(1024);
    double* acc = partial[tid].data();

    for (int i = tid; i < n; i += T) {
      if (!obs[i]) continue;
      const double ti = vt[i], xi = vx[i], yi = vy[i];

      // Pass 1: shared pairwise kernel weights from earlier events.
      wbuf.clear(); jbuf.clear();
      for (int j = i - 1; j >= 0; --j) {
        const double dt = ti - vt[j];
        if (dt <= 0.0) continue;
        if (dt > dt_cut) break;
        const double dx = xi - vx[j];
        const double dy = yi - vy[j];
        const double term = pair_term(dt, dx * dx + dy * dy);
        if (term <= 0.0) continue;
        jbuf.push_back(j);
        wbuf.push_back(term);
      }
      const int m = static_cast<int>(jbuf.size());

      // Pass 2: mix by labelling membership (no transcendentals).
      for (int kk = 0; kk < Klab; ++kk) {
        if (!has_obs[kk]) continue;
        const int* mk = mem.data() + static_cast<size_t>(kk) * n;
        if (!mk[i]) continue;
        double lam = mu_base * vW[i];
        for (int s = 0; s < m; ++s) {
          if (mk[jbuf[s]]) lam += wbuf[s];
        }
        if (lam <= 1e-15) lam = 1e-15;
        acc[kk] += std::log(lam);
      }
    }
  };

  if (T == 1) {
    worker(0);
  } else {
    std::vector<std::thread> pool;
    pool.reserve(T);
    for (int tid = 0; tid < T; ++tid) pool.emplace_back(worker, tid);
    for (auto& th : pool) th.join();
  }

  // Shared compensator CDF increments per event.
  std::vector<double> Gcontrib(n, 0.0);
  for (int j = 0; j < n; ++j) {
    const double p_t = vt[j];
    double start_h = t_start;
    if (p_t > start_h) start_h = p_t;
    double end_h = t_end;
    if (do_trunc) {
      const double trunc_end = p_t + t_trunc;
      if (trunc_end < end_h) end_h = trunc_end;
    }
    if (end_h <= start_h) continue;
    Gcontrib[j] = K * (temporal_cdf(end_h - p_t) - temporal_cdf(start_h - p_t)) / temporal_norm;
  }

  for (int kk = 0; kk < Klab; ++kk) {
    if (!has_obs[kk]) {
      out[kk] = -1e15;
      continue;
    }
    double loglik = 0.0;
    for (int tid = 0; tid < T; ++tid) loglik += partial[tid][kk];
    const int* mk = mem.data() + static_cast<size_t>(kk) * n;
    double trig = 0.0;
    for (int j = 0; j < n; ++j) {
      if (mk[j]) trig += Gcontrib[j];
    }
    loglik -= (bg_comp + trig);
    if (!std::isfinite(loglik)) loglik = -1e15;
    out[kk] = loglik;
  }
  return out;
}
