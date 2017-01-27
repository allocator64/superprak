#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#include <cassert>
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

#define DEBUG(var)                     \
  do {                                 \
    if (state::rank == 0) {            \
      LOG_INFO << #var << ": " << var; \
    }                                  \
  } while (false)

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

struct Request {
  Request(int buf_len) : buf(buf_len), mpi_req() {}
  std::vector<double> buf;
  MPI_Request mpi_req;

  void Done() { CALL_MPI(MPI_Wait, &mpi_req, MPI_STATUS_IGNORE); }

  void Send(int id) {
    CALL_MPI(MPI_Isend, buf.data(), buf.size(), MPI_DOUBLE, id, 0,
             MPI_COMM_WORLD, &mpi_req);
  }

  void Receive(int id) {
    CALL_MPI(MPI_Irecv, buf.data(), buf.size(), MPI_DOUBLE, id, 0,
             MPI_COMM_WORLD, &mpi_req);
  }
};

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

double Tau(Matrix& r, Matrix& g, RangeMapper x_grid, RangeMapper y_grid,
           double* dst_up = NULL, double* dst_down = NULL) {
  double up = 0;
  double down = 0;
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

double Alpha(Matrix& r, Matrix& g, RangeMapper x_grid, RangeMapper y_grid,
             double* dst_up = NULL, double* dst_down = NULL) {
  double up = 0;
  double down = 0;
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
  // DEBUG(state::x_processors_num);
  // DEBUG(state::y_processors_num);
}

void FillCurrentNode() {
  state::x_node = state::rank % state::x_processors_num;
  state::y_node = state::rank / state::x_processors_num;
  // DEBUG(state::x_node);
  // DEBUG(state::y_node);
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
  int neighbor = (state::x_node + x_diff) +
                 (state::y_node + y_diff) * state::x_processors_num;
  return neighbor;
}

struct Neighbor {
  struct Area {
    Area(int x_first, int x_second, int y_first, int y_second)
        : x_first(x_first),
          x_second(x_second),
          y_first(y_first),
          y_second(y_second) {}

    int size() { return (x_second - x_first) * (y_second - y_first); }
    int x_first;
    int x_second;
    int y_first;
    int y_second;
  };

  Neighbor(Area area, bool need_send, int id)
      : area(area), need_send(need_send), id(id) {}

  Area area;
  bool need_send;
  int id;
};

std::vector<Neighbor> neighbors;
typedef std::vector<Neighbor>::iterator neighbors_iterator;

mpi::Request* MakeRequest(Matrix& m, Neighbor& n) {
  Neighbor::Area& area = n.area;
  mpi::Request* req = new mpi::Request(area.size());
  if (n.need_send) {
    int idx = 0;
    for (int i = area.x_first; i < area.x_second; ++i) {
      for (int j = area.y_first; j < area.y_second; ++j) {
        req->buf.at(idx++) = m.at(i, j);
      }
    }
    req->Send(n.id);
  } else {
    req->Receive(n.id);
  }
  return req;
}

void FinishRequest(Neighbor& n, Matrix* m, mpi::Request* req) {
  req->Done();
  if (!n.need_send) {
    int idx = 0;
    Neighbor::Area& area = n.area;
    // assert(area.size() == req->buf.size());
    for (int i = area.x_first; i < area.x_second; ++i) {
      for (int j = area.y_first; j < area.y_second; ++j) {
        m->at(i, j) = req->buf[idx++];
      }
    }
  }
  delete req;
}

void SyncWithNeighbors(Matrix* m) {
  std::vector<mpi::Request*> requests;
  // requests.reserve(neighbors.size());
  for (neighbors_iterator n = neighbors.begin(); n != neighbors.end(); ++n) {
    requests.push_back(MakeRequest(*m, *n));
  }
  for (int idx = 0; idx < requests.size(); ++idx) {
    FinishRequest(neighbors[idx], m, requests[idx]);
  }
}

void AllReduceSum(double* up, double* down) {
  static double buf[2];
  buf[0] = *up;
  buf[1] = *down;
  CALL_MPI(MPI_Allreduce, MPI_IN_PLACE, buf, /* count */ 2, MPI_DOUBLE, MPI_SUM,
           MPI_COMM_WORLD);
  *up = buf[0];
  *down = buf[1];
}

void AllReduceMax(double* diff) {
  CALL_MPI(MPI_Allreduce, MPI_IN_PLACE, diff, /* count */ 1, MPI_DOUBLE,
           MPI_MAX, MPI_COMM_WORLD);
}

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
  // DEBUG(x_grid.row_offet);
  // DEBUG(y_grid.row_offet);

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
    comm::neighbors.push_back(
        comm::Neighbor(comm::Neighbor::Area(1, state::x_node_points_num - 1,
                                            state::y_node_points_num - 1,
                                            state::y_node_points_num),
                       /* need_send */ false, neighbor));
  }

  // for (int i = 0; i < comm::neighbors.size(); ++i) {
  //   DEBUG(i);
  //   DEBUG(comm::neighbors[i].id);
  //   DEBUG(comm::neighbors[i].need_send);
  //   DEBUG(comm::neighbors[i].area.x_first);
  //   DEBUG(comm::neighbors[i].area.x_second);
  //   DEBUG(comm::neighbors[i].area.y_first);
  //   DEBUG(comm::neighbors[i].area.y_second);
  // }

  Matrix p(state::x_node_points_num, state::y_node_points_num);
  Matrix r(state::x_node_points_num, state::y_node_points_num);

  using statement::Phi;

  for (int i = 0; i < state::x_node_points_num; ++i) {
    if (state::y_node == 0) {
      p.at(i, 0) = Phi(x_grid[i], y_grid[0]);
    }
    if (state::y_node == state::y_processors_num - 1) {
      p.at(i, state::y_node_points_num - 1) =
          Phi(x_grid[i], y_grid[state::y_node_points_num - 1]);
    }
  }

  for (int j = 0; j < state::y_node_points_num; ++j) {
    if (state::x_node == 0) {
      p.at(0, j) = Phi(x_grid[0], y_grid[j]);
    }
    if (state::x_node == state::x_processors_num - 1) {
      p.at(state::x_node_points_num - 1, j) =
          Phi(x_grid[state::x_node_points_num - 1], y_grid[j]);
    }
  }

  double up = 0;
  double down = 0;
  double global_diff = 0;

  math::ApplyR(p, x_grid, y_grid, &r);
  comm::SyncWithNeighbors(&r);

  math::Tau(r, r, x_grid, y_grid, &up, &down);
  comm::AllReduceSum(&up, &down);
  double tau = up / down;

  math::ApplyP(r, tau, &p);
  comm::SyncWithNeighbors(&p);

  Matrix g = r;
  int counter = 0;
  for (double diff = 1e100; diff > statement::eps;) {
    math::ApplyR(p, x_grid, y_grid, &r);
    comm::SyncWithNeighbors(&r);

    math::Alpha(r, g, x_grid, y_grid, &up, &down);
    comm::AllReduceSum(&up, &down);
    double alpha = up / down;

    math::ApplyG(r, alpha, &g);
    comm::SyncWithNeighbors(&g);

    math::Tau(r, g, x_grid, y_grid, &up, &down);
    comm::AllReduceSum(&up, &down);
    tau = up / down;

    global_diff = math::ApplyP(g, tau, &p);
    comm::SyncWithNeighbors(&p);
    comm::AllReduceMax(&global_diff);
    diff = global_diff;

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
  comm::AllReduceMax(&err);

  DEBUG(err);
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
    diff = math::ApplyP(g, tau, &p);
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

  // std::stringstream ss;
  // ss << "stderr-" << std::setfill('0') << std::setw(2) << state::rank;
  // if (!freopen(ss.str().c_str(), "a", stderr)) {
  //   throw std::runtime_error("freopen error");
  // }

  DEBUG(state::rank);
  DEBUG(state::process_num);

  state::x_points_num = atoi(argv[1]);
  state::y_points_num = atoi(argv[2]);

  // DEBUG(statement::area.x.first);
  // DEBUG(statement::area.x.second);
  // DEBUG(statement::area.y.first);
  // DEBUG(statement::area.y.second);

  if (state::process_num == 1) {
    solution::Simple();
  } else {
    solution::Super();
  }

  double finish_time = mpi::Time();
  DEBUG(finish_time - start_time);
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

  return 0;
}
