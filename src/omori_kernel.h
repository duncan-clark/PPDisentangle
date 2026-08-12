#pragma once
#include <cmath>

// Truncated Omori-Utsu helpers. Untruncated Omori is a density on [0, inf)
// only for p > 1. On a finite interval [0, t_trunc] the kernel
// (1 + dt/c)^{-p} is integrable for every p (c > 0); p = 1 is logarithmic.
//
// Intensity: g(dt) = pref * (1 + dt/c)^{-p}
// CDF: F(h) = I(min(h, T)) / I(T) when truncated, else the improper G(h).

inline bool omori_trunc_active(double t_trunc) {
  return t_trunc > 0.0;
}

inline double omori_time_prefactor(double p, double c, double t_trunc,
                                   bool do_trunc) {
  if (!do_trunc) return (p - 1.0) / c;
  const double uT = 1.0 + t_trunc / c;
  if (std::fabs(p - 1.0) < 1e-8) {
    const double logu = std::log(uT);
    if (!(logu > 0.0)) return 0.0;
    return 1.0 / (c * logu);
  }
  double G = 1.0 - std::pow(uT, -(p - 1.0));
  if (std::fabs(G) < 1e-15) G = (G < 0.0) ? -1e-15 : 1e-15;
  return (p - 1.0) / (c * G);
}

inline double omori_time_cdf(double p, double c, double h, double t_trunc,
                             bool do_trunc) {
  if (h <= 0.0) return 0.0;
  if (do_trunc && h > t_trunc) h = t_trunc;
  const double uh = 1.0 + h / c;
  if (std::fabs(p - 1.0) < 1e-8) {
    if (!do_trunc) return 0.0;
    const double den = std::log(1.0 + t_trunc / c);
    if (!(den > 0.0)) return 0.0;
    return std::log(uh) / den;
  }
  const double Gh = 1.0 - std::pow(uh, -(p - 1.0));
  if (!do_trunc) return Gh;
  double GT = 1.0 - std::pow(1.0 + t_trunc / c, -(p - 1.0));
  if (std::fabs(GT) < 1e-15) GT = (GT < 0.0) ? -1e-15 : 1e-15;
  return Gh / GT;
}

// Inverse-CDF sample. u is Uniform(0, 1). Truncated draws land in [0, t_trunc].
inline double omori_sample_delay(double p, double c, double t_trunc,
                                 bool do_trunc, double u) {
  if (u <= 0.0) return 0.0;
  if (u >= 1.0) u = 1.0 - 1e-15;
  if (do_trunc && std::fabs(p - 1.0) < 1e-8) {
    return c * (std::pow(1.0 + t_trunc / c, u) - 1.0);
  }
  double a = 0.0;
  if (do_trunc) a = std::pow(1.0 + t_trunc / c, 1.0 - p);
  const double inner = 1.0 + u * (a - 1.0);
  if (!(inner > 0.0)) return do_trunc ? t_trunc : 1e15;
  return c * (std::pow(inner, 1.0 / (1.0 - p)) - 1.0);
}
