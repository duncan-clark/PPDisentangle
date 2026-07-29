#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector tile_index_rect_cpp(NumericVector x,
                                  NumericVector y,
                                  NumericVector xgrid,
                                  NumericVector ygrid) {
  int n = x.size();
  int nx = xgrid.size() - 1;
  int ny = ygrid.size() - 1;
  IntegerVector out(n);

  if (nx <= 0 || ny <= 0) {
    std::fill(out.begin(), out.end(), NA_INTEGER);
    return out;
  }

  double x_min = xgrid[0];
  double x_max = xgrid[nx];
  double y_min = ygrid[0];
  double y_max = ygrid[ny];
  double dx = (x_max - x_min) / nx;
  double dy = (y_max - y_min) / ny;

  for (int i = 0; i < n; ++i) {
    double xi = x[i];
    double yi = y[i];
    if (!R_finite(xi) || !R_finite(yi) ||
        xi < x_min || xi > x_max || yi < y_min || yi > y_max) {
      out[i] = NA_INTEGER;
      continue;
    }

    int col = (xi == x_max) ? nx : static_cast<int>(std::floor((xi - x_min) / dx)) + 1;
    int row_from_bottom = (yi == y_max) ? ny : static_cast<int>(std::floor((yi - y_min) / dy)) + 1;
    int row = ny - row_from_bottom + 1;
    if (col < 1 || col > nx || row < 1 || row > ny) {
      out[i] = NA_INTEGER;
    } else {
      out[i] = (row - 1) * nx + col;
    }
  }

  return out;
}
