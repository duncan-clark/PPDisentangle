// Batched bivariate ETAS log-likelihood over multiple labellings.
// Labellings share (t, x, y, mag); only process labels and background masks
// differ. Pairwise temporal/spatial kernel weights are computed once and
// reused. Each returned value matches etas_bivariate_loglik_cpp for that
// labelling (to floating-point rounding).
//
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <thread>
#include <algorithm>
#include "omori_kernel.h"
using namespace Rcpp;

static const double PI_VAL = 3.14159265358979323846;

//' Batched bivariate ETAS log-likelihood over labellings
//'
//' @param t Event times (shared across labellings, window starts at 0).
//' @param x,y Spatial coordinates (shared).
//' @param mag Event magnitudes (shared).
//' @param process_ids Integer matrix n x K: 0 = control, 1 = treated.
//' @param W0s,W1s Numeric matrices n x K of background weights, or n x 1
//'   to recycle the same mask across all labellings.
//' @param mu_0,mu_1 Background rates.
//' @param A_00,alpha_m_00,A_11,alpha_m_11,A_01,alpha_m_01,A_10,alpha_m_10
//'   Productivity parameters.
//' @param cc,p,D,gamma_par,q Shared structural parameters.
//' @param m0 Reference magnitude.
//' @param areaS_0,areaS_1 Active areas.
//' @param t_max Observation window length.
//' @param t_trunc Temporal truncation (negative to disable).
//' @param bg_exposure_0,bg_exposure_1 Time span over which each background is
//'   switched on (compensator charge). Negative (default) means the full
//'   window \code{t_max}.
//' @param n_threads Number of worker threads (1 = serial).
//' @return Length-K numeric vector of log-likelihoods.
//' @keywords internal
// [[Rcpp::export]]
NumericVector etas_bivariate_loglik_batch_cpp(
    NumericVector t, NumericVector x, NumericVector y, NumericVector mag,
    IntegerMatrix process_ids,
    NumericMatrix W0s, NumericMatrix W1s,
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
    double bg_exposure_1 = -1.0,
    int n_threads = 1) {

  const int n = t.size();
  const int K = process_ids.ncol();
  NumericVector out(K);
  if (n == 0) {
    for (int kk = 0; kk < K; ++kk) out[kk] = -1e15;
    return out;
  }
  if (n_threads < 1) n_threads = 1;

  const bool do_trunc = (t_trunc > 0.0);
  if (areaS_0 <= 0.0) areaS_0 = 1.0;
  if (areaS_1 <= 0.0) areaS_1 = 1.0;
  const double mu_base_0 = mu_0 / areaS_0;
  const double mu_base_1 = mu_1 / areaS_1;

  const double omori_pref = omori_time_prefactor(p, cc, t_trunc, do_trunc);
  const double base_const = omori_pref * (q - 1.0) / PI_VAL;

  double A_mat[2][2], alpha_mat[2][2];
  A_mat[0][0] = A_00;   alpha_mat[0][0] = alpha_m_00;
  A_mat[0][1] = A_01;   alpha_mat[0][1] = alpha_m_01;
  A_mat[1][0] = A_10;   alpha_mat[1][0] = alpha_m_10;
  A_mat[1][1] = A_11;   alpha_mat[1][1] = alpha_m_11;

  // Plain buffers: worker threads must not touch R memory management.
  // W0s/W1s may be n x K or n x 1 (shared masks recycled across labellings).
  const int w0_cols = W0s.ncol();
  const int w1_cols = W1s.ncol();
  const int w0_stride = (w0_cols == 1) ? 0 : n;
  const int w1_stride = (w1_cols == 1) ? 0 : n;
  std::vector<double> vt(t.begin(), t.end()), vx(x.begin(), x.end()),
      vy(y.begin(), y.end());
  std::vector<int> pid(static_cast<size_t>(n) * K);
  std::vector<double> W0(static_cast<size_t>(w0_cols == 1 ? n : n * K));
  std::vector<double> W1(static_cast<size_t>(w1_cols == 1 ? n : n * K));
  for (int kk = 0; kk < K; ++kk) {
    for (int i = 0; i < n; ++i) {
      pid[static_cast<size_t>(kk) * n + i] = process_ids(i, kk);
    }
  }
  if (w0_cols == 1) {
    for (int i = 0; i < n; ++i) W0[static_cast<size_t>(i)] = W0s(i, 0);
  } else {
    for (int kk = 0; kk < K; ++kk) {
      for (int i = 0; i < n; ++i) {
        W0[static_cast<size_t>(kk) * n + i] = W0s(i, kk);
      }
    }
  }
  if (w1_cols == 1) {
    for (int i = 0; i < n; ++i) W1[static_cast<size_t>(i)] = W1s(i, 0);
  } else {
    for (int kk = 0; kk < K; ++kk) {
      for (int i = 0; i < n; ++i) {
        W1[static_cast<size_t>(kk) * n + i] = W1s(i, kk);
      }
    }
  }

  std::vector<double> d_par(n);
  std::vector<double> E00(n), E01(n), E10(n), E11(n);
  const double inv_cc = 1.0 / cc;
  for (int j = 0; j < n; ++j) {
    const double dm = mag[j] - m0;
    d_par[j] = D * std::exp(gamma_par * dm);
    E00[j] = (A_mat[0][0] < 1e-20) ? 0.0 : A_mat[0][0] * std::exp(alpha_mat[0][0] * dm);
    E01[j] = (A_mat[0][1] < 1e-20) ? 0.0 : A_mat[0][1] * std::exp(alpha_mat[0][1] * dm);
    E10[j] = (A_mat[1][0] < 1e-20) ? 0.0 : A_mat[1][0] * std::exp(alpha_mat[1][0] * dm);
    E11[j] = (A_mat[1][1] < 1e-20) ? 0.0 : A_mat[1][1] * std::exp(alpha_mat[1][1] * dm);
    if (E00[j] < 1e-20) E00[j] = 0.0;
    if (E01[j] < 1e-20) E01[j] = 0.0;
    if (E10[j] < 1e-20) E10[j] = 0.0;
    if (E11[j] < 1e-20) E11[j] = 0.0;
  }
  const double* Echan[4] = { E00.data(), E01.data(), E10.data(), E11.data() };

  const int T = std::min<int>(n_threads, std::max(1, n));
  std::vector<std::vector<double>> partial(T, std::vector<double>(K, 0.0));

  auto worker = [&](int tid) {
    std::vector<double> wbuf;
    std::vector<int> jbuf;
    wbuf.reserve(1024); jbuf.reserve(1024);
    double* acc = partial[tid].data();

    for (int i = tid; i < n; i += T) {
      const double ti = vt[i], xi = vx[i], yi = vy[i];

      wbuf.clear(); jbuf.clear();
      if (do_trunc) {
        for (int j = i - 1; j >= 0; --j) {
          const double dt = ti - vt[j];
          if (dt > t_trunc) break;
          const double dx = xi - vx[j];
          const double dy = yi - vy[j];
          const double r2 = dx * dx + dy * dy;
          const double w = std::exp(-p * std::log(1.0 + dt * inv_cc)
                                    - q * std::log(1.0 + r2 / d_par[j])) / d_par[j];
          jbuf.push_back(j);
          wbuf.push_back(base_const * w);
        }
      } else {
        for (int j = 0; j < i; ++j) {
          const double dt = ti - vt[j];
          if (dt <= 0.0) continue;
          const double dx = xi - vx[j];
          const double dy = yi - vy[j];
          const double r2 = dx * dx + dy * dy;
          const double w = std::exp(-p * std::log(1.0 + dt * inv_cc)
                                    - q * std::log(1.0 + r2 / d_par[j])) / d_par[j];
          jbuf.push_back(j);
          wbuf.push_back(base_const * w);
        }
      }
      const int m = static_cast<int>(jbuf.size());

      for (int kk = 0; kk < K; ++kk) {
        const int* pidk = pid.data() + static_cast<size_t>(kk) * n;
        const int k = pidk[i];
        double lam = (k == 0)
          ? mu_base_0 * W0[static_cast<size_t>(kk) * w0_stride + i]
          : mu_base_1 * W1[static_cast<size_t>(kk) * w1_stride + i];
        for (int s = 0; s < m; ++s) {
          const int j = jbuf[s];
          lam += Echan[k * 2 + pidk[j]][j] * wbuf[s];
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

  // Compensator temporal mass per parent: G(h) / G(t_trunc). The division by
  // temporal_norm matches the truncation-renormalized kernel used in the
  // intensity (A = expected offspring within t_trunc); without it the
  // compensator undercharges productivity by the factor G(t_trunc).
  std::vector<double> G_h(n);
  for (int j = 0; j < n; ++j) {
    double horizon = t_max - vt[j];
    if (do_trunc && horizon > t_trunc) horizon = t_trunc;
    G_h[j] = (horizon <= 0.0) ? 0.0
      : omori_time_cdf(p, cc, horizon, t_trunc, do_trunc);
  }

  // Background compensator: charge each mu_k only over the time span its
  // background is actually on (policy masks can shorten the exposure).
  const double expo_0 = (bg_exposure_0 >= 0.0) ? bg_exposure_0 : t_max;
  const double expo_1 = (bg_exposure_1 >= 0.0) ? bg_exposure_1 : t_max;

  for (int kk = 0; kk < K; ++kk) {
    double loglik = 0.0;
    for (int tid = 0; tid < T; ++tid) loglik += partial[tid][kk];
    const int* pidk = pid.data() + static_cast<size_t>(kk) * n;
    double comp_trig = 0.0;
    for (int j = 0; j < n; ++j) {
      if (G_h[j] <= 0.0) continue;
      const int l = pidk[j];
      const double e0 = Echan[0 * 2 + l][j];
      const double e1 = Echan[1 * 2 + l][j];
      if (e0 >= 1e-20) comp_trig += e0 * G_h[j];
      if (e1 >= 1e-20) comp_trig += e1 * G_h[j];
    }
    loglik -= (mu_0 * expo_0 + mu_1 * expo_1 + comp_trig);
    if (!std::isfinite(loglik)) loglik = -1e15;
    out[kk] = loglik;
  }
  return out;
}
