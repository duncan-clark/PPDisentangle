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

  // Base background density (scalar part)
  // mu is the total background rate over the whole window (events per unit time).
  // The background intensity is lambda_0(x,y) = (mu / areaS) * W(x,y)
  // where W is normalized such that its average value is 1.
  double mu_base = mu / areaS;

  // Triggering constant
  // The kernel is h(dt,x,y) = const_val * exp(-beta*dt - alpha*r2)
  // Integral of exp(-beta*dt) from 0 to inf is 1/beta
  // Integral of exp(-alpha*r2) over R^2 is pi/alpha
  // So integral of triggering kernel is const_val * (1/beta) * (pi/alpha)
  // We want this integral to be K.
  // So const_val = K * alpha * beta / pi
  double const_val = K * alpha * beta / pi;

  for(int i = 0; i < n; ++i) {

    // 1. Background Intensity at this specific point
    double lambda_i = mu_base * W_val[i];

    // 2. Triggering Intensity (Sum over past events)
    for(int j = 0; j < i; ++j) {
      double dt = t[i] - t[j];

      // Optimization: Time truncation
      if(dt * beta > 20.0) continue;

      double dx = x[i] - x[j];
      double dy = y[i] - y[j];
      double r2 = dx*dx + dy*dy;

      // Optimization: Space truncation
      if(r2 * alpha > 20.0) continue;

      double excitation = std::exp(-beta * dt - alpha * r2);
      lambda_i += const_val * excitation;
    }

    // Safety check to avoid log(0) or log(negative)
    if(lambda_i <= 1e-15) lambda_i = 1e-15;

    loglik += std::log(lambda_i);
  }

  // INTEGRAL TERM (Negative part of log-likelihood)
  // loglik = sum(log(lambda(t_i))) - integral(lambda(t) dt)
  // The integral of lambda(t,x,y) over [0,T] x S is:
  // mu * T + sum_{i=1}^N K * (1 - exp(-beta * (T - t_i)))
  
  double triggering_integral = 0;
  for(int i = 0; i < n; ++i) {
    triggering_integral += K * (1.0 - std::exp(-beta * (t_max - t[i])));
  }

  loglik -= (mu * t_max + triggering_integral);

  if (NumericVector::is_na(loglik) || std::isinf(loglik)) return -1e15;

  return loglik;
}
