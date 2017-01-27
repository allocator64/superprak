#pragma once

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
