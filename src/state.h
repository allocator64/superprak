#pragma once

#include <string>

namespace state {
extern int rank;
extern int process_num;
extern int x_points_num;
extern int y_points_num;

extern int x_processors_num;
extern int y_processors_num;

extern int x_node;
extern int y_node;

extern int x_node_points_num;
extern int y_node_points_num;

extern std::string output_path;
}  // namespace state
