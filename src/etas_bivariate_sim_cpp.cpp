#include <Rcpp.h>
using namespace Rcpp;

//' Simulate bivariate ETAS offspring via BFS branching
//'
//' Each parent in process k generates offspring in BOTH processes via
//' the 2x2 kernel matrix.  Offspring counts are:
//'   n_{k->l} ~ Poisson(A_{lk} * exp(alpha_m_{lk} * (m - m0)))
//' for each child process l in {0, 1}.
//'
//' Temporal delays, spatial displacements, and magnitude assignment
//' are identical to the univariate ETAS simulation (shared structural
//' parameters c, p, D, gamma, q).
//'
//' @param parent_x,parent_y,parent_t,parent_mag  Parent coordinates.
//' @param parent_process  Integer vector: 0 = control, 1 = treated.
//' @param A_00,alpha_m_00  Control self-excitation.
//' @param A_11,alpha_m_11  Treated self-excitation.
//' @param A_01,alpha_m_01  Treated-to-control cross-excitation.
//' @param A_10,alpha_m_10  Control-to-treated cross-excitation.
//' @param cc,p,D,gamma_par,q  Shared structural parameters.
//' @param m0         Reference magnitude.
//' @param beta_gr    Gutenberg-Richter beta (negative to use mag_pool).
//' @param t_min,t_max  Temporal window.
//' @param x_min,x_max,y_min,y_max  Spatial bounding box.
//' @param t_trunc    Temporal truncation (negative to disable).
//' @param mag_pool   Magnitude pool for resampling.
//' @return List with x, y, t, mag, process_id of all offspring.
// [[Rcpp::export]]
List sim_etas_bivariate_children_cpp(
    NumericVector parent_x, NumericVector parent_y,
    NumericVector parent_t, NumericVector parent_mag,
    IntegerVector parent_process,
    double A_00, double alpha_m_00,
    double A_11, double alpha_m_11,
    double A_01, double alpha_m_01,
    double A_10, double alpha_m_10,
    double cc, double p, double D, double gamma_par, double q,
    double m0, double beta_gr,
    double t_min, double t_max,
    double x_min, double x_max, double y_min, double y_max,
    double t_trunc = -1.0,
    NumericVector mag_pool = NumericVector::create()) {

  bool do_trunc = (t_trunc > 0.0);
  bool use_gr = (beta_gr > 0.0);
  int pool_n = mag_pool.size();

  double cdf_max = do_trunc ?
    (1.0 - std::pow(1.0 + t_trunc / cc, -(p - 1.0))) : 1.0;

  double pm1_inv = 1.0 / (p - 1.0);
  double qm1_inv = 1.0 / (q - 1.0);
  const double two_pi = 2.0 * 3.14159265358979323846;

  // Kernel matrix: A_mat[child][parent], alpha_mat[child][parent]
  double A_mat[2][2], alpha_mat[2][2];
  A_mat[0][0] = A_00;   alpha_mat[0][0] = alpha_m_00;
  A_mat[0][1] = A_01;   alpha_mat[0][1] = alpha_m_01;
  A_mat[1][0] = A_10;   alpha_mat[1][0] = alpha_m_10;
  A_mat[1][1] = A_11;   alpha_mat[1][1] = alpha_m_11;

  // BFS queue
  std::vector<double> q_x = as<std::vector<double>>(parent_x);
  std::vector<double> q_y = as<std::vector<double>>(parent_y);
  std::vector<double> q_t = as<std::vector<double>>(parent_t);
  std::vector<double> q_m = as<std::vector<double>>(parent_mag);
  std::vector<int>    q_p = as<std::vector<int>>(parent_process);

  std::vector<double> out_x, out_y, out_t, out_m;
  std::vector<int>    out_p;
  int est = parent_t.size() * 4;
  if (est < 100) est = 100;
  out_x.reserve(est); out_y.reserve(est);
  out_t.reserve(est); out_m.reserve(est);
  out_p.reserve(est);

  size_t head = 0;

  while (head < q_x.size()) {
    double px = q_x[head];
    double py = q_y[head];
    double pt = q_t[head];
    double pm = q_m[head];
    int    pp = q_p[head];
    head++;

    double dm = pm - m0;
    double d_parent = D * std::exp(gamma_par * dm);

    // Generate offspring for each child process
    for (int child = 0; child < 2; ++child) {
      double A_kl = A_mat[child][pp];
      if (A_kl < 1e-20) continue;
      double alpha_kl = alpha_mat[child][pp];
      double kappa = A_kl * std::exp(alpha_kl * dm);

      int n_kids = R::rpois(kappa);
      if (n_kids == 0) continue;

      for (int k = 0; k < n_kids; ++k) {
        double u_t = R::runif(0.0, 1.0) * cdf_max;
        double dt = cc * (std::pow(1.0 - u_t, -pm1_inv) - 1.0);
        double new_t = pt + dt;
        if (new_t > t_max || new_t < t_min) continue;

        double u_s = R::runif(0.0, 1.0);
        double r2 = d_parent * (std::pow(1.0 - u_s, -qm1_inv) - 1.0);
        double r = std::sqrt(r2);
        double angle = R::runif(0.0, two_pi);
        double new_x = px + r * std::cos(angle);
        double new_y = py + r * std::sin(angle);

        if (new_x < x_min || new_x > x_max ||
            new_y < y_min || new_y > y_max) continue;

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

        q_x.push_back(new_x); q_y.push_back(new_y);
        q_t.push_back(new_t); q_m.push_back(new_mag);
        q_p.push_back(child);

        out_x.push_back(new_x); out_y.push_back(new_y);
        out_t.push_back(new_t); out_m.push_back(new_mag);
        out_p.push_back(child);
      }
    }
  }

  return List::create(
    Named("x") = out_x,
    Named("y") = out_y,
    Named("t") = out_t,
    Named("mag") = out_m,
    Named("process_id") = out_p
  );
}
