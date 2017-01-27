#pragma once

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
