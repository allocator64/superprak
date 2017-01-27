#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "logger.h"
#include "mpi.h"
#include "solution.h"
#include "state.h"

void run(int argc, char** argv) {
  CALL_MPI(MPI_Init, &argc, &argv);

  double start_time = mpi::Time();

  CALL_MPI(MPI_Comm_rank, MPI_COMM_WORLD, &state::rank);
  CALL_MPI(MPI_Comm_size, MPI_COMM_WORLD, &state::process_num);

  if (state::rank != 0) {
    if (!freopen("/dev/null", "a", stdout)) {
      throw std::runtime_error("freopen error");
    }
    if (!freopen("/dev/null", "a", stderr)) {
      throw std::runtime_error("freopen error");
    }
  }
  // std::stringstream ss;
  // ss << "stderr-" << std::setfill('0') << std::setw(2) << state::rank;
  // if (!freopen(ss.str().c_str(), "a", stderr)) {
  //   throw std::runtime_error("freopen error");
  // }

  // DEBUG(state::rank);
  // DEBUG(state::process_num);

  state::x_points_num = atoi(argv[1]);
  state::y_points_num = atoi(argv[2]);
  state::output_path = argv[3] ? argv[3] : "";

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
