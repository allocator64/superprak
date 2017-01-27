#pragma once

#include <mpi.h>
#include <sstream>
#include <vector>

#include "logger.h"

namespace mpi {
void Call(const char* func_name, int func_result);
}

#define CALL_MPI(func, args...) mpi::Call(#func, func(args))

namespace mpi {
inline double Time() {
  CALL_MPI(MPI_Barrier, MPI_COMM_WORLD);
  return MPI_Wtime();
}

inline void Call(const char* func_name, int func_result) {
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
