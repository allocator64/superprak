#include "communications.h"
#include "logger.h"
#include "math.h"
#include "matrix.h"
#include "range_mapper.h"
#include "solution.h"
#include "state.h"
#include "statement.h"
#include "utils.h"

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
  // DEBUG(state::x_node_points_num);
  // DEBUG(state::y_node_points_num);
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

    counter++;
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

}  // namespace solution
