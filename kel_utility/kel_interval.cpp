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


std::pair<size_t, size_t> kel::OpenRightInterval::intersection(const OpenRightInterval &interval) const {

  if (lower_ >= interval.upper_ or interval.lower_ >= upper_) {

    return {0, 0};

  }

  return { std::max(lower_, interval.lower_), std::min(upper_, interval.upper_ ) };

}

