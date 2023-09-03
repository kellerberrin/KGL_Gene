//
// Created by kellerberrin on 15/06/23.
//

#include "kel_interval.h"
#include "kel_exec_env.h"

namespace kel = kellerberrin;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kel::operator==(const OpenRightInterval& lhs, const OpenRightInterval &rhs) {

  return lhs.lower() == rhs.lower() and lhs.upper() == rhs.upper();

}

bool kel::operator<(const OpenRightInterval &lhs, const OpenRightInterval &rhs) {

  return lhs.lower() < rhs.lower();

}

// Modify the interval.
void kel::OpenRightInterval::resize(size_t lower, size_t upper)  {

  if (upper < lower) {

    ExecEnv::log().warn("OpenRightInterval::resize; Incorrect Initialization, Upper Offset: {} < Lower Offset: {}", upper, lower);
    std::swap(lower, upper);

  }

  lower_ = lower;
  upper_ = upper;

}

// Shift the interval without changing it's size.
void kel::OpenRightInterval::translate(int64_t shift) {

  if (shift < 0 and std::abs(shift) > static_cast<int64_t>(lower())) {

    ExecEnv::log().warn("OpenRightInterval::translate; translate shift: {} results in negative values for interval [{}, {})", shift, lower_, upper_);
    return;

  }

  lower_ += shift;
  upper_ += shift;

}


// Shift the interval so that lower() == 0. The translation (which is always non-positive) is returned.
int64_t kel::OpenRightInterval::translateZero() {

  int64_t shift =  -1 * static_cast<int64_t>(lower());
  translate(shift);

  return shift;

}


kel::OpenRightInterval kel::OpenRightInterval::intersection(const OpenRightInterval &interval) const {

  if (lower_ >= interval.upper_ or interval.lower_ >= upper_) {

    return {0, 0};

  }

  return { std::max(lower_, interval.lower_), std::min(upper_, interval.upper_ ) };

}

bool kel::OpenRightInterval::containsInterval(const OpenRightInterval &interval) const {

  return intersection(interval) == interval;

}

kel::OpenRightInterval kel::OpenRightInterval::merge(const OpenRightInterval &interval) const {

  if (intersects(interval) or adjacent(interval)) {

    return { std::min(lower_, interval.lower_), std::max(upper_, interval.upper_ ) };

  }

  return {0, 0};

}
