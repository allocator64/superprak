#pragma once

#include "geo_primitives.h"

class RangeMapper {
 public:
  RangeMapper(Range range, int total_row_len, int row_offet = 0)
      : step((range.second - range.first) / total_row_len),
        row_offet(row_offet),
        range_(range) {}

  double operator[](int idx) const {
    // if (idx < 0 || idx > parts_) {
    //   std::stringstream ss;
    //   ss << "idx: " << idx << " parts: " << parts_;
    //   throw std::out_of_range(ss.str());
    // }
    return step * (idx + row_offet) + range_.first;
  }

  const double step;
  const int row_offet;

 private:
  Range range_;
};
