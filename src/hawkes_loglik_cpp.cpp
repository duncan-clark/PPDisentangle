#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double hawkes_loglik_inhom_cpp(NumericVector t,
                               NumericVector x,
                               NumericVector y,
                               NumericVector W_val,
                               double mu, double alpha, double beta, double K,
                               double areaS, double t_max) {

  int n = t.size();
  double loglik = 0.0;

  double pi = 3.14159265358979323846;

  double mu_base = mu / areaS;

  double const_val = K * (alpha / pi) * beta;

  double dt_cutoff = 15.0 / beta;

  for(int i = 0; i < n; ++i) {

    double lambda_i = mu_base * W_val[i];

    double t_lo = t[i] - dt_cutoff;

    int j_start = 0;
    if(t_lo > t[0]) {
      int lo = 0, hi = i;
      while(lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if(t[mid] < t_lo) lo = mid + 1;
        else hi = mid;
      }
      j_start = lo;
    }

    for(int j = j_start; j < i; ++j) {
      double dt = t[i] - t[j];

      double dx = x[i] - x[j];
      double dy = y[i] - y[j];
      double r2 = dx*dx + dy*dy;

      if(r2 * alpha > 15.0) continue;

      double excitation = exp(-beta * dt - alpha * r2);
      lambda_i += const_val * excitation;
    }

    if(lambda_i <= 1e-10) lambda_i = 1e-10;

    loglik += log(lambda_i);
  }

  return loglik;
}
