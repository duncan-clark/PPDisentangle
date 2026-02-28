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
                                  double t_trunc = -1.0) {

  bool do_trunc = (t_trunc > 0.0);

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

  // CDF at t_trunc for truncated exponential inverse-CDF sampling
  double cdf_max = do_trunc ? (1.0 - std::exp(-beta * t_trunc)) : 0.0;

  size_t head = 0;

  while(head < q_x.size()) {
    double px = q_x[head];
    double py = q_y[head];
    double pt = q_t[head];
    head++;

    int n_kids = R::rpois(K);
    if(n_kids == 0) continue;

    for(int k = 0; k < n_kids; ++k) {
      double dt;
      if(do_trunc) {
        double u = R::runif(0.0, 1.0);
        dt = -std::log(1.0 - u * cdf_max) / beta;
      } else {
        dt = R::rexp(1.0/beta);
      }
      double new_t = pt + dt;

      if(new_t > t_max) continue;
      if(new_t < t_min) continue;

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
