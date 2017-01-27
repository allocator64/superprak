#pragma once

#include "matrix.h"
#include "range_mapper.h"
#include "statement.h"

namespace math {
inline double Laplas(const Matrix& m, const RangeMapper& x_grid,
                     const RangeMapper& y_grid, int x, int y) {
  double ldx = (m.at(x, y) - m.at(x - 1, y)) / x_grid.step;
  double rdx = (m.at(x + 1, y) - m.at(x, y)) / x_grid.step;
  double ldy = (m.at(x, y) - m.at(x, y - 1)) / y_grid.step;
  double rdy = (m.at(x, y + 1) - m.at(x, y)) / y_grid.step;
  double dx = (ldx - rdx) / x_grid.step;
  double dy = (ldy - rdy) / y_grid.step;
  return dx + dy;
}

inline void ApplyR(const Matrix& p, const RangeMapper& x_grid,
                   const RangeMapper& y_grid, Matrix* r) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 1; i < r->x_size - 1; ++i) {
    for (int j = 1; j < r->y_size - 1; ++j) {
      r->at(i, j) =
          Laplas(p, x_grid, y_grid, i, j) - statement::F(x_grid[i], y_grid[j]);
    }
  }
}

inline double Tau(const Matrix& r, const Matrix& g, const RangeMapper& x_grid,
                  const RangeMapper& y_grid, double* dst_up = NULL,
                  double* dst_down = NULL) {
  double up = 0;
  double down = 0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+ : up, down)
#endif
  for (int i = 1; i < r.x_size - 1; ++i) {
    for (int j = 1; j < r.y_size - 1; ++j) {
      double tmp = g.at(i, j);
      up += r.at(i, j) * tmp;
      down += Laplas(g, x_grid, y_grid, i, j) * tmp;
    }
  }
  if (dst_up) {
    *dst_up = up;
  }
  if (dst_down) {
    *dst_down = down;
  }
  return up / down;
}

inline double ApplyP(const Matrix& g, double tau, Matrix* p) {
  double diff = 0;
#ifdef USE_OPENMP
// #pragma omp parallel for reduction(max : diff)
#endif
  for (int i = 1; i < p->x_size - 1; ++i) {
    for (int j = 1; j < p->y_size - 1; ++j) {
      double tmp = p->at(i, j) - tau * g.at(i, j);
      diff = std::max(diff, fabs(tmp - p->at(i, j)));
      p->at(i, j) = tmp;
    }
  }
  return diff;
}

inline double Alpha(const Matrix& r, const Matrix& g, const RangeMapper& x_grid,
                    const RangeMapper& y_grid, double* dst_up = NULL,
                    double* dst_down = NULL) {
  double up = 0;
  double down = 0;
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+ : up, down)
#endif
  for (int i = 1; i < r.x_size - 1; ++i) {
    for (int j = 1; j < r.y_size - 1; ++j) {
      up += Laplas(r, x_grid, y_grid, i, j) * g.at(i, j);
      down += Laplas(g, x_grid, y_grid, i, j) * g.at(i, j);
    }
  }
  if (dst_up) {
    *dst_up = up;
  }
  if (dst_down) {
    *dst_down = down;
  }
  return up / down;
}

inline void ApplyG(const Matrix& r, double alpha, Matrix* g) {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 1; i < r.x_size - 1; ++i) {
    for (int j = 1; j < r.y_size - 1; ++j) {
      g->at(i, j) = r.at(i, j) - alpha * g->at(i, j);
    }
  }
}

}  // namespace math
