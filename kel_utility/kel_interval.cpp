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


kel::OpenRightInterval::OpenRightInterval(size_t lower, size_t upper)  {

  if (upper <= lower) {

    ExecEnv::log().warn("OpenRightInterval::OpenRightInterval, Incorrect Initialization, Upper Offset: {} <= Lower Offset: {}", upper, lower);
    if (upper == lower) {

      ++upper;

    } else {

      std::swap(lower, upper);

    }


  }

  lower_ = lower;
  upper_ = upper;

}


// Returns true if the interval argument is contained within one of the intervals held in the set.
bool kel::IntervalSet::containsInterval(const OpenRightInterval& interval) const {

  auto iter = this->lower_bound(interval);
  if (iter != this->end()) {

    if (iter->lower() == interval.lower()) {

      return iter->containsInterval(interval);

    }

  }

  // Check if contained by the previous interval (if it exists).
  iter = std::prev(iter, 1);
  if (iter != end()) {

    return iter->containsInterval(interval);

  }

  return false;

}

