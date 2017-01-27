#include <fstream>

#include "logger.h"
#include "math.h"
#include "matrix.h"
#include "range_mapper.h"
#include "solution.h"
#include "state.h"
#include "statement.h"

namespace solution {

void Simple() {
  RangeMapper x_grid(statement::area.x, state::x_points_num);
  RangeMapper y_grid(statement::area.y, state::y_points_num);

  Matrix p(state::x_points_num, state::y_points_num);
  Matrix r(state::x_points_num, state::y_points_num);

  for (int i = 0; i < state::x_points_num; ++i) {
    p.at(i, 0) = statement::Phi(x_grid[i], y_grid[0]);
    p.at(i, state::y_points_num - 1) =
        statement::Phi(x_grid[i], y_grid[state::y_points_num - 1]);
  }

  for (int j = 0; j < state::y_points_num; ++j) {
    p.at(0, j) = statement::Phi(x_grid[0], y_grid[j]);
    p.at(state::x_points_num - 1, j) =
        statement::Phi(x_grid[state::x_points_num - 1], y_grid[j]);
  }

  math::ApplyR(p, x_grid, y_grid, &r);
  double tau = math::Tau(r, r, x_grid, y_grid);
  math::ApplyP(r, tau, &p);

  Matrix g = r;
  int counter = 0;
  for (double diff = 1e100; diff > statement::eps;) {
    math::ApplyR(p, x_grid, y_grid, &r);
    double alpha = math::Alpha(r, g, x_grid, y_grid);
    math::ApplyG(r, alpha, &g);
    tau = math::Tau(r, g, x_grid, y_grid);
    diff = math::ApplyP(g, tau, &p);
    // DEBUG(counter++);
    // DEBUG(diff);
  }

  double err = 0;
  for (int i = 1; i < p.x_size - 1; ++i) {
    for (int j = 1; j < p.y_size - 1; ++j) {
      err = std::max(err,
                     fabs(statement::Phi(x_grid[i], y_grid[j]) - p.at(i, j)));
    }
  }

  DEBUG(err);

  if (!state::output_path.empty()) {
    std::ofstream out(state::output_path.c_str());
    for (int i = 0; i < p.x_size; ++i) {
      for (int j = 0; j < p.y_size; ++j) {
        out << x_grid[i] << " " << y_grid[j] << " " << p.at(i, j) << std::endl;
      }
    }
  }
}

}  // namespace solution
