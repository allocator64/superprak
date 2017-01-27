#pragma once

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
