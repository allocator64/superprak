#include <mpi.h>
#include <stdio.h>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

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

static int x_processors_num;
static int y_processors_num;

static int x_node;
static int y_node;

static int x_node_points_num;
static int y_node_points_num;
}

#define LOG_ERROR                                                 \
  std::cerr << "\n\033[1;31m" __FILE__ << ":" << __LINE__ << " [" \
            << state::rank << "]\033[0m "

#define LOG_INFO \
  std::cerr << "\n" __FILE__ << ":" << __LINE__ << " [" << state::rank << "] "

#define DEBUG(var) LOG_INFO << #var << ": " << var

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
  RangeMapper(Range range, int total_row_len, int row_offet = 0)
      : step((range.second - range.first) / total_row_len),
        row_offet(row_offet),
        range_(range),
        total_row_len_(total_row_len) {}

  double operator[](int idx) {
    // if (idx < 0 || idx > parts_) {
    //   std::stringstream ss;
    //   ss << "idx: " << idx << " parts: " << parts_;
    //   throw std::out_of_range(ss.str());
    // }
    return step * (idx + row_offet) + range_.first;
  }

  const double step;
  const int row_offet;

 private:
  Range range_;
  int total_row_len_;
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
      r->at(i, j) =
          Laplas(p, x_grid, y_grid, i, j) - statement::F(x_grid[i], y_grid[j]);
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

namespace utils {
void FillAllProcessors() {
  int power2 = 1;
  while ((1 << power2) < state::process_num) ++power2;
  state::x_processors_num = 1 << (power2 / 2);
  state::y_processors_num = state::process_num / state::x_processors_num;
  DEBUG(state::x_processors_num);
  DEBUG(state::y_processors_num);
}

void FillCurrentNode() {
  state::x_node = state::rank % state::x_processors_num;
  state::y_node = state::rank / state::x_processors_num;
  DEBUG(state::x_node);
  DEBUG(state::y_node);
}

RangeMapper AdjustCell(Range real_row, int real_row_points,
                       int total_row_processors, int cur_position,
                       int* cur_len) {
  int points_per_processor = real_row_points / total_row_processors;
  int additional_points = real_row_points % total_row_processors;

  int row_offet = points_per_processor * cur_position;
  if (cur_position < additional_points) {
    row_offet += cur_position;
  } else {
    row_offet += additional_points;
  }
  if (cur_position != 0) {
    row_offet--;
  }

  *cur_len = points_per_processor;
  if (cur_position < additional_points) {
    (*cur_len)++;
  }
  if (cur_position != 0) {
    (*cur_len)++;
  }
  if (cur_position != total_row_processors - 1) {
    (*cur_len)++;
  }

  return RangeMapper(real_row, real_row_points, row_offet);
}

}  // namespace utils

namespace comm {
int GetNeighborNode(int x_diff, int y_diff) {
  return (state::x_processors_num + x_diff) +
         (state::y_processors_num + y_diff) * state::x_processors_num;
}

struct Neighbor {
  struct Area {
    Area(int x_first, int x_second, int y_first, int y_second)
        : x_first(x_first),
          x_second(x_second),
          y_first(y_first),
          y_second(y_second) {
      DEBUG(x_first);
      DEBUG(x_second);
      DEBUG(y_first);
      DEBUG(y_second);
    }
    int x_first;
    int x_second;
    int y_first;
    int y_second;
  };

  Neighbor(Area area, bool need_send, int id)
      : area(area), need_send(need_send), id(id) {
    DEBUG(need_send);
    DEBUG(id);
  }

  Area area;
  bool need_send;
  int id;
};

std::vector<Neighbor> neighbors;
}  // namespace comm

namespace solution {

void Super() {
  utils::FillAllProcessors();
  utils::FillCurrentNode();

  RangeMapper x_grid = utils::AdjustCell(statement::area.x, state::x_points_num,
                                         state::x_processors_num, state::x_node,
                                         &state::x_node_points_num);
  RangeMapper y_grid = utils::AdjustCell(statement::area.y, state::y_points_num,
                                         state::y_processors_num, state::y_node,
                                         &state::y_node_points_num);
  DEBUG(state::x_node_points_num);
  DEBUG(state::y_node_points_num);
  DEBUG(x_grid.row_offet);
  DEBUG(y_grid.row_offet);

  if (state::x_node != 0) {
    int neighbor = comm::GetNeighborNode(-1, 0);
    comm::neighbors.push_back(comm::Neighbor(
        comm::Neighbor::Area(0, 1, 1, state::y_node_points_num - 1),
        /* need_send */ false, neighbor));
    comm::neighbors.push_back(comm::Neighbor(
        comm::Neighbor::Area(1, 2, 1, state::y_node_points_num - 1),
        /* need_send */ true, neighbor));
  }
  if (state::y_node != 0) {
    int neighbor = comm::GetNeighborNode(0, -1);
    comm::neighbors.push_back(comm::Neighbor(
        comm::Neighbor::Area(1, state::x_node_points_num - 1, 0, 1),
        /* need_send */ false, neighbor));
    comm::neighbors.push_back(comm::Neighbor(
        comm::Neighbor::Area(1, state::x_node_points_num - 1, 1, 2),
        /* need_send */ true, neighbor));
  }
  if (state::x_node != state::x_processors_num - 1) {
    int neighbor = comm::GetNeighborNode(1, 0);
    comm::neighbors.push_back(
        comm::Neighbor(comm::Neighbor::Area(state::x_node_points_num - 2,
                                            state::x_node_points_num - 1, 1,
                                            state::y_node_points_num - 1),
                       /* need_send */ true, neighbor));
    comm::neighbors.push_back(
        comm::Neighbor(comm::Neighbor::Area(state::x_node_points_num - 1,
                                            state::x_node_points_num, 1,
                                            state::y_node_points_num - 1),
                       /* need_send */ false, neighbor));
  }
  if (state::y_node != state::y_processors_num - 1) {
    int neighbor = comm::GetNeighborNode(0, 1);
    comm::neighbors.push_back(
        comm::Neighbor(comm::Neighbor::Area(1, state::x_node_points_num - 1,
                                            state::y_node_points_num - 2,
                                            state::y_node_points_num - 1),
                       /* need_send */ true, neighbor));
    comm::neighbors.push_back(comm::Neighbor(
        comm::Neighbor::Area(1, state::x_node_points_num - 1,
                             state::y_points_num - 1, state::y_points_num),
        /* need_send */ false, neighbor));
  }
}

void Simple() {
  RangeMapper x_grid(statement::area.x, state::x_points_num);
  RangeMapper y_grid(statement::area.y, state::y_points_num);

  Matrix p(state::x_points_num, state::y_points_num);
  Matrix r(state::x_points_num, state::y_points_num);

  LOG_INFO << "Initial iteration";

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
  for (double diff = 1e100; diff > statement::eps;) {
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
      err = std::max(err,
                     fabs(statement::Phi(x_grid[i], y_grid[j]) - p.at(i, j)));
    }
  }

  DEBUG(err);

  std::ofstream out("output.dat");
  for (int i = 0; i < p.x_size; ++i) {
    for (int j = 0; j < p.y_size; ++j) {
      out << x_grid[i] << " " << y_grid[j] << " " << p.at(i, j) << std::endl;
    }
  }
}

}  // namespace solution

void run(int argc, char** argv) {
  CALL_MPI(MPI_Init, &argc, &argv);

  double start_time = mpi::Time();

  CALL_MPI(MPI_Comm_rank, MPI_COMM_WORLD, &state::rank);
  CALL_MPI(MPI_Comm_size, MPI_COMM_WORLD, &state::process_num);
  if (state::rank != 0) {
    std::stringstream ss;
    ss << "stderr-" << std::setfill('0') << std::setw(2) << state::rank;
    if (!freopen(ss.str().c_str(), "a", stderr)) {
      throw std::runtime_error("freopen error");
    }
    std::stringstream ss2;
    ss2 << "stdout-" << std::setfill('0') << std::setw(2) << state::rank;
    if (!freopen(ss2.str().c_str(), "a", stdout)) {
      throw std::runtime_error("freopen error");
    }
  }
  LOG_INFO << "-- rank: " << state::rank;
  LOG_INFO << "-- process_num: " << state::process_num;

  state::x_points_num = atoi(argv[1]);
  state::y_points_num = atoi(argv[2]);

  DEBUG(statement::area.x.first);
  DEBUG(statement::area.x.second);
  DEBUG(statement::area.y.first);
  DEBUG(statement::area.y.second);

  if (state::process_num == 1) {
    solution::Simple();
  } else {
    solution::Super();
  }

  double finish_time = mpi::Time();
  LOG_INFO << "Run time: " << finish_time - start_time;
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
    LOG_ERROR << "Exception: " << e.what();
    CALL_MPI(MPI_Abort, MPI_COMM_WORLD, 1);
    return 2;
  }

  std::cout << std::endl;
  return 0;
}
