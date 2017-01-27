#pragma once

#include <cmath>

#include "geo_primitives.h"

namespace statement {

extern const double eps;

extern const Rectangle area;

inline double F(double x, double y) {
  double tmp = 1 + x * y;
  return (x * x + y * y) / tmp / tmp;
}

inline double Phi(double x, double y) { return log(1 + x * y); }
}  // namespace statement
