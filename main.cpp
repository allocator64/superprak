// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com

#include <mpi.h>
#include <stdio.h>
#include <exception>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <cstdlib>

struct Point {
  double x;
  double y;
};

struct Range {
  double first;
  double second;

  Range(double first, double second) : first(first), second(second) {}

  double GetLength() { return second - first; }
};

struct Rectangle {
  Range x;
  Range y;

  Rectangle(Range x, Range y) : x(x), y(y) {}
};

namespace statement {
static double eps = 1e-4; 
Rectangle area(Range(0, 3), Range(0, 3));
double F(double x, double y) {
  double tmp = 1 + x * y;
  return (x * x + y * y) / tmp / tmp;
}

double Phi(double x, double y) { return log(1 + x * y); }
}

namespace state {
static int rank;
static int process_num;
static int x_points_num;
static int y_points_num;
}

#define LOG_ERROR \
  std::cerr << __FILE__ << ":" << __LINE__ << " (" << state::rank << ") "

#define LOG_INFO \
  std::cout << __FILE__ << ":" << __LINE__ << " (" << state::rank << ") "

#define DEBUG(var) LOG_INFO << #var << ": " << var << std::endl

#define CALL_MPI(func, args...) mpi::Call(#func, func(args))

namespace mpi {
void Call(const char* func_name, int func_result);
double Time() {
  CALL_MPI(MPI_Barrier, MPI_COMM_WORLD);
  return MPI_Wtime();
}

void Call(const char* func_name, int func_result) {
  if (func_result != MPI_SUCCESS) {
    std::stringstream ss;
    ss << "MPI failed: " << func_name << " code: " << func_result;
    throw std::runtime_error(ss.str());
  }
}
}  // namespace mpi

class RangeMapper {
 public:
  RangeMapper(Range range, int parts)
      : step((range.second - range.first) / parts),
        range_(range),
        parts_(parts)
         {}

  double operator[](int idx) {
    if (idx < 0 || idx > parts_) {
      std::stringstream ss;
      ss << "idx: " << idx << " parts: " << parts_;
      throw std::out_of_range(ss.str());
    }
    return step * idx + range_.first;
  }

  const double step;

 private:
  Range range_;
  int parts_;
};

class Matrix {
 public:
  Matrix(int x_size, int y_size)
      : x_size(x_size), y_size(y_size), data_(x_size * y_size, 0.0) {}

  double& at(int x, int y) { return data_.at(x * y_size + y); }

  const int x_size;
  const int y_size;

 private:
  std::vector<double> data_;
};

namespace math {
double Laplas(Matrix& m, RangeMapper x_grid, RangeMapper y_grid, int x, int y) {
  double ldx = (m.at(x, y) - m.at(x - 1, y)) / x_grid.step;
  double rdx = (m.at(x + 1, y) - m.at(x, y)) / x_grid.step;
  double ldy = (m.at(x, y) - m.at(x, y - 1)) / y_grid.step;
  double rdy = (m.at(x, y + 1) - m.at(x, y)) / y_grid.step;
  double dx = (ldx - rdx) / x_grid.step;
  double dy = (ldy - rdy) / y_grid.step;
  return dx + dy;
}

void ApplyR(Matrix& p, RangeMapper x_grid, RangeMapper y_grid, Matrix* r) {
  for (int i = 1; i < r->x_size - 1; ++i) {
    for (int j = 1; j < r->y_size - 1; ++j) {
      r->at(i, j) = Laplas(p, x_grid, y_grid, i, j) - statement::F(x_grid[i], y_grid[j]);
    }
  }
}

double Tau(Matrix& r, Matrix& g, RangeMapper x_grid, RangeMapper y_grid) {
  double up = 0;
  double down = 0;
  for (int i = 1; i < r.x_size - 1; ++i) {
    for (int j = 1; j < r.y_size - 1; ++j) {
      double tmp = g.at(i, j);
      up += r.at(i, j) * tmp;
      down += Laplas(g, x_grid, y_grid, i, j) * tmp;
    }
  }
  return up / down;
}

double ApplyP(Matrix& g, double tau, Matrix* p) {
  double diff = 0;
  for (int i = 1; i < p->x_size - 1; ++i) {
    for (int j = 1; j < p->y_size - 1; ++j) {
      double tmp = p->at(i, j) - tau * g.at(i, j);
      diff = std::max(diff, fabs(tmp - p->at(i, j)));
      p->at(i, j) = tmp;
    }
  }
  return diff;
}

double Alpha(Matrix& r, Matrix& g, RangeMapper x_grid, RangeMapper y_grid) {
  double up = 0;
  double down = 0;
  for (int i = 1; i < r.x_size - 1; ++i) {
    for (int j = 1; j < r.y_size - 1; ++j) {
      up += Laplas(r, x_grid, y_grid, i, j) * g.at(i, j);
      down += Laplas(g, x_grid, y_grid, i, j) * g.at(i, j);
    }
  }
  return up / down;
}

void ApplyG(Matrix& r, double alpha, Matrix* g) {
  for (int i = 1; i < r.x_size - 1; ++i) {
    for (int j = 1; j < r.y_size - 1; ++j) {
      g->at(i, j) = r.at(i, j) - alpha * g->at(i, j);
    }
  } 
}

}  // namespace math


void run(int argc, char** argv) {
  CALL_MPI(MPI_Init, &argc, &argv);

  double start_time = mpi::Time();

  CALL_MPI(MPI_Comm_rank, MPI_COMM_WORLD, &state::rank);
  CALL_MPI(MPI_Comm_size, MPI_COMM_WORLD, &state::process_num);
  if (state::rank != 0) {
    if (!freopen("stdout", "w", stdout) || !freopen("stderr", "w", stderr)) {
      throw std::runtime_error("freopen error");
    }
  }
  LOG_INFO << "-- rank: " << state::rank << std::endl;
  LOG_INFO << "-- process_num: " << state::process_num << std::endl;

  state::x_points_num = atoi(argv[1]);
  state::y_points_num = atoi(argv[2]);

  DEBUG(statement::area.x.first);
  DEBUG(statement::area.x.second);
  DEBUG(statement::area.y.first);
  DEBUG(statement::area.y.second);

  RangeMapper x_grid(statement::area.x, state::x_points_num);
  RangeMapper y_grid(statement::area.y, state::y_points_num);

  Matrix p(state::x_points_num, state::y_points_num);
  Matrix r(state::x_points_num, state::y_points_num);

  LOG_INFO << "Initial iteration" << std::endl;

  using statement::Phi;

  for (int i = 0; i < state::x_points_num; ++i) {
    p.at(i, 0) = Phi(x_grid[i], y_grid[0]);
    p.at(i, state::y_points_num - 1) =
        Phi(x_grid[i], y_grid[state::y_points_num - 1]);
  }

  for (int j = 0; j < state::y_points_num; ++j) {
    p.at(0, j) = Phi(x_grid[0], y_grid[j]);
    p.at(state::x_points_num - 1, j) =
        Phi(x_grid[state::x_points_num - 1], y_grid[j]);
  }

  math::ApplyR(p, x_grid, y_grid, &r);
  double tau = math::Tau(r, r, x_grid, y_grid);
  math::ApplyP(r, tau, &p);

  Matrix g = r;
  int counter = 0;
  for (double diff = 1e100; diff > statement::eps; ) {
    math::ApplyR(p, x_grid, y_grid, &r);
    double alpha = math::Alpha(r, g, x_grid, y_grid);
    math::ApplyG(r, alpha, &g);
    tau = math::Tau(r, g, x_grid, y_grid);
    diff = sqrt(math::ApplyP(g, tau, &p) * x_grid.step * y_grid.step);
    DEBUG(counter++);
    DEBUG(diff);
  }

  double err = 0;
  for (int i = 1; i < p.x_size - 1; ++i) {
    for (int j = 1; j < p.y_size - 1; ++j) {
      err = std::max(err, fabs(statement::Phi(x_grid[i], y_grid[j]) - p.at(i, j)));
    }
  }

  DEBUG(err);

  std::ofstream out("output.dat");
  for (int i = 0; i < p.x_size; ++i) {
    for (int j = 0; j < p.y_size; ++j) {
      out << x_grid[i] << " " << y_grid[j] << " " << p.at(i, j) << std::endl;
    }
  }

  double finish_time = mpi::Time();
  LOG_INFO << "Run time: " << finish_time - start_time << std::endl;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " x_points_num y_points_num"
              << std::endl;
    return 1;
  }

  try {
    run(argc, argv);
  } catch (std::exception& e) {
    LOG_ERROR << "Exception: " << e.what() << std::endl;
    CALL_MPI(MPI_Abort, MPI_COMM_WORLD, 1);
    return 2;
  }

  return 0;
}
