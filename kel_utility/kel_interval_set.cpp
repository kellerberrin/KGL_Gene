//
// Created by kellerberrin on 17/06/23.
//

#include "kel_interval_set.h"
#include "kel_exec_env.h"



namespace kel = kellerberrin;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Returns true if the interval argument is contained within one of the intervals held in the set.
bool kel::IntervalSetLower::containsInterval(const OpenRightInterval& interval) const {

  if (empty()) {

    return false;

  }

  auto const iter = lower_bound(interval);
  if (iter != end()) {

    if (iter->lower() == interval.lower()) {

      return iter->containsInterval(interval);

    }

  }

  // Look at previous intervals.
  auto prev_iter = std::ranges::prev(iter, 1, begin());
  // Check that iter was not at this->begin()
  if (prev_iter != iter) {

    while(interval.lower() < prev_iter->upper()) {

      if (prev_iter->containsInterval(interval)) {

        return true;

      }

      if (prev_iter == begin()) {

        break;

      }
      prev_iter = std::ranges::prev(prev_iter, 1, begin());

    }

  }

  return false;

}


// Returns the interval unions between this set and the argument set.
// Intervals are modified/extended as necessary for a disjoint interval minimal union.
kel::IntervalSetLower kel::IntervalSetLower::intervalSetUnion(const IntervalSetLower& interval_set) const {

  IntervalSetLower set_union;
  IntervalMultiSetLower sorted_intervals;

  if (empty() and interval_set.empty()) {

    return set_union;

  }

  // Place all intervals into a multiset to sort them.
  for (auto const& interval : *this) {

    sorted_intervals.insert(interval);

  }

  for (auto const& interval : interval_set) {

    sorted_intervals.insert(interval);

  }

  auto iter = sorted_intervals.begin();
  // Test for the empty condition.
  if (iter != sorted_intervals.end()) {

    OpenRightInterval active_interval = *iter;
    // Loop through the multiset and aggregate any overlapping intervals.
    iter = std::ranges::next(iter, 1, sorted_intervals.end());
    while (iter != sorted_intervals.end()) {

      OpenRightInterval next_interval = *iter;

      OpenRightInterval merged_interval = active_interval.merge(next_interval);

      if (not merged_interval.empty()) {

        active_interval = merged_interval;

      } else {

        auto [set_iter, result] = set_union.insert(active_interval);
        if (not result) {

          ExecEnv::log().warn("IntervalSetLower::intervalSetUnion; Cannot add duplicate interval: [{}, {})",
                              active_interval.lower(), active_interval.upper());

        }

        active_interval = next_interval;

      }

      iter = std::ranges::next(iter, 1, sorted_intervals.end());

    } // Loop through sorted intervals.

    auto [set_iter, result] = set_union.insert(active_interval);
    if (not result) {

      ExecEnv::log().warn("IntervalSetLower::intervalSetUnion; Cannot add duplicate interval: [{}, {})",
                          active_interval.lower(), active_interval.upper());

    }

  } // Empty condition.

  return set_union;

}


// Find the set intervals that intersect the argument interval.
[[nodiscard]] std::vector<kel::OpenRightInterval> kel::IntervalSetLower::findIntersectsInterval(const OpenRightInterval &interval) const {

  std::vector<OpenRightInterval> value_vector;

  if (empty()) {

    return value_vector;

  }

  auto const bound_iter = lower_bound(interval); // iter->first.lower() >= arg.lower()

  // Scroll forward to check intervals that intersect
  auto next_iter = bound_iter;
  if (next_iter != end()) {

    while (next_iter->lower() < interval.upper()) {

      auto const intersect_interval = next_iter->intersection(interval);
      if (not intersect_interval.empty()) {

        value_vector.push_back(*next_iter);

      }

      next_iter = std::ranges::next(next_iter, 1, end());
      if (next_iter == end()) {

        break;

      }

    }

  }

  // Scroll backwards to find intervals that intersect.
  auto prev_iter = std::ranges::prev(bound_iter, 1, begin());
  // Check that bound_iter was not this->begin()
  if (prev_iter != bound_iter) {

    while (interval.lower() < prev_iter->upper()) {

      // Is there an intersection between then argument interval and the previous interval.
      auto const intersect_interval = prev_iter->intersection(interval);
      // If so, then return the previous interval.
      if (not intersect_interval.empty()) {

        value_vector.push_back(*prev_iter);

      }

      if (prev_iter == begin()) {

        break;

      }
      prev_iter = std::ranges::prev(prev_iter, 1 , begin());

    } // while

  } // Not end()

  // Else just return the upper interval (which may be end())
  return value_vector;

}

