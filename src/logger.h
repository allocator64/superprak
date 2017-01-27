#pragma once

#include <iostream>
#include "mpi.h"
#include "state.h"

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
