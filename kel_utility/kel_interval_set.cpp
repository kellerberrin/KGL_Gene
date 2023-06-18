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


// Find the interval that contains the argument OR the interval immediately greater (lower > arg.lower) than the argument interval.
kel::IntervalSet::const_iterator kel::IntervalSet::findUpperEqualIter(const OpenRightInterval &interval) const {

  auto iter = this->lower_bound(interval);
  // Return the previous map interval if it contains the argument interval.
  auto prev_iter = std::prev(iter, 1);
  if (prev_iter != this->end()) {

    auto const &interval_key = *prev_iter;
    if (interval_key.containsInterval(interval)) {

      return prev_iter;

    }

  }

  // Just return the upper interval
  return iter;

}

kel::IntervalSet kel::IntervalSet::intervalSetIntersection(const IntervalSet& second_set) const {

  IntervalSet set_intersection;

  auto iter = this->begin();
  auto second_iter = second_set.begin();

  while(iter != this->end() and second_iter != second_set.end()) {

    auto const [lower, upper] = iter->intersection(*second_iter);
    if (lower < upper) {

      auto [insert_iter, result] = set_intersection.insert(OpenRightInterval(lower, upper));
      if (not result) {

        ExecEnv::log().warn("IntervalSet::intervalSetIntersection; Cannot add intersection interval: [{}, {}) (duplicate)", lower, upper);

      }

    }

    if (iter->upper() < second_iter->upper()) {

      ++iter;

    } else {

      ++second_iter;

    }

  }

  return set_intersection;

}


kel::IntervalSet kel::IntervalSet::intervalSetUnion(const IntervalSet& interval_set) const {

  IntervalSet set_union;
  IntervalMultiSet sorted_intervals;

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
    ++iter;
    while (iter != sorted_intervals.end()) {

      // Have three possibilities 1. contained within active, 2. extends active, 3. is disjoint with active.
      if (not active_interval.containsInterval(*iter)) {

        if (iter->lower() <= active_interval.upper() and iter->upper() > active_interval.upper()) {

          active_interval = OpenRightInterval(active_interval.lower(), iter->upper());

        } else {
          // Disjoint
          auto [set_iter, result] = set_union.insert(active_interval);
          if (not result) {

            ExecEnv::log().warn("IntervalSet::intervalSetUnion; Cannot add duplicate interval: [{}, {})",
                                active_interval.lower(), active_interval.upper());
            auto find_iter = set_union.find(active_interval);
            if (find_iter != set_union.end()) {

              auto const& found_interval = *find_iter;
              ExecEnv::log().warn("IntervalSet::intervalSetUnion; Previously inserted interval: [{}, {})",
                                  found_interval.lower(), found_interval.upper());

              for (auto const& interval : sorted_intervals) {

                ExecEnv::log().info("IntervalSet::intervalSetUnion; Sorted interval multiset: [{}, {})", interval.lower(), interval.upper());

              }

              for (auto const& interval : set_union) {

                ExecEnv::log().info("IntervalSet::intervalSetUnion; Union interval set: [{}, {})", interval.lower(), interval.upper());

              }

            }

          }

          active_interval = *iter;

        }

      }

      ++iter;

    } // Loop through sorted intervals.

    auto [set_iter, result] = set_union.insert(active_interval);
    if (not result) {

      ExecEnv::log().warn("IntervalSet::intervalSetUnion; Cannot add duplicate interval: [{}, {})",
                          active_interval.lower(), active_interval.upper());

    }

  } // Empty condition.

  return set_union;

}

