#include "communications.h"

#include <vector>

#include "matrix.h"
#include "mpi.h"
#include "state.h"

namespace comm {
int GetNeighborNode(int x_diff, int y_diff) {
  int neighbor = (state::x_node + x_diff) +
                 (state::y_node + y_diff) * state::x_processors_num;
  return neighbor;
}

std::vector<Neighbor> neighbors;

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
  for (size_t idx = 0; idx < requests.size(); ++idx) {
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
