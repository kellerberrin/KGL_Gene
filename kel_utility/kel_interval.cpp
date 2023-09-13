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


bool kel::operator==(const OpenRightUnsigned& lhs, const OpenRightUnsigned &rhs) {

  return lhs.lower() == rhs.lower() and lhs.upper() == rhs.upper();

}

bool kel::operator<(const OpenRightUnsigned &lhs, const OpenRightUnsigned &rhs) {

  return lhs.lower() < rhs.lower();

}

// Modify the interval.
void kel::OpenRightUnsigned::resize(size_t lower, size_t upper)  {

  if (upper < lower) {

    ExecEnv::log().warn("OpenRightUnsigned::resize; Incorrect Initialization, Upper Offset: {} < Lower Offset: {}", upper, lower);
    std::swap(lower, upper);

  }

  lower_ = lower;
  upper_ = upper;

}

// Shift the interval without changing it's size.
kel::OpenRightUnsigned kel::OpenRightUnsigned::translate(int64_t shift) const {

  OpenRightUnsigned translated(*this);
  if ((static_cast<int64_t>(lower()) + shift) < 0) {

    ExecEnv::log().warn("OpenRightUnsigned::translate; translate shift: {} results in negative values for interval: {}", shift, toString());
    return {lower(), upper()};

  }

  size_t translate_lower = static_cast<int64_t>(lower()) + shift;
  size_t translate_upper = static_cast<int64_t>(upper()) + shift;

  return { translate_lower , translate_upper };

}


// Shift the interval so that lower() == 0. The translation (which is always non-positive) is returned.
kel::OpenRightUnsigned kel::OpenRightUnsigned::translateZero() const {

  int64_t shift =  -1 * static_cast<int64_t>(lower());
  return translate(shift);

}


kel::OpenRightUnsigned kel::OpenRightUnsigned::intersection(const OpenRightUnsigned &interval) const {

  if (lower_ >= interval.upper_ or interval.lower_ >= upper_) {

    return {0, 0};

  }

  return { std::max(lower_, interval.lower_), std::min(upper_, interval.upper_ ) };

}

bool kel::OpenRightUnsigned::containsInterval(const OpenRightUnsigned &interval) const {

  return intersection(interval) == interval;

}

kel::OpenRightUnsigned kel::OpenRightUnsigned::merge(const OpenRightUnsigned &interval) const {

  if (intersects(interval) or adjacent(interval)) {

    return { std::min(lower_, interval.lower_), std::max(upper_, interval.upper_ ) };

  }

  return {0, 0};

}
