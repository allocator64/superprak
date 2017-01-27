#pragma once

#include <vector>

#include "matrix.h"
#include "mpi.h"

namespace comm {
int GetNeighborNode(int x_diff, int y_diff);

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

extern std::vector<Neighbor> neighbors;
typedef std::vector<Neighbor>::iterator neighbors_iterator;

mpi::Request* MakeRequest(Matrix& m, Neighbor& n);

void FinishRequest(Neighbor& n, Matrix* m, mpi::Request* req);

void SyncWithNeighbors(Matrix* m);

void AllReduceSum(double* up, double* down);

void AllReduceMax(double* diff);

}  // namespace comm
