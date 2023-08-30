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
    shift =  -1 * static_cast<int64_t>(lower());

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

// To insert an interval the lower() parameter of inserted interval must be within the range [lower, upper).
kel::OpenRightInterval kel::OpenRightInterval::insertInterval(const OpenRightInterval &insert_interval) const {

  if (containsOffset(insert_interval.lower())) {

    return { lower(), upper() + insert_interval.size()};

  }

  return { lower(), upper()};

}

// For a valid delete the intersection of the delete interval must be non-empty.
// Note that the delete_interval argument is modified.
kel::OpenRightInterval kel::OpenRightInterval::deleteInterval(const OpenRightInterval &delete_interval) const {

  // does not delete this interval.
  if (disjoint(delete_interval)) {

    return { lower(), upper()};

  }

  if (containsInterval(delete_interval) and delete_interval.lower() >= lower()) {

    return { lower(), upper()-delete_interval.size() };

  }

  auto intersect_interval = intersection(delete_interval);

  // If the delete interval extends beyond the modified interval.
  if (delete_interval.lower() >= lower()) {

    return { lower(), upper() - intersect_interval.size() };

  }

  // Else the delete interval is below the modified interval.
  size_t upper_adjust = (lower() - delete_interval.lower()) + intersect_interval.size();
  return { delete_interval.lower(), upper() - upper_adjust };

}

