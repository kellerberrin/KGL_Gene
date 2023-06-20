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

  if (this->empty()) {

    return false;

  }

  auto const iter = this->lower_bound(interval);
  if (iter != this->end()) {

    if (iter->lower() == interval.lower()) {

      return iter->containsInterval(interval);

    }

  }

  // Look at previous intervals.
  auto prev_iter = std::ranges::prev(iter, 1, this->begin());
  // Check that iter was not at this->begin()
  if (prev_iter != iter) {

    while(interval.lower() < prev_iter->upper()) {

      if (prev_iter->containsInterval(interval)) {

        return true;

      }

      if (prev_iter == this->begin()) {

        break;

      }
      prev_iter = std::ranges::prev(prev_iter, 1, this->begin());

    }

  }

  return false;

}


kel::IntervalSet kel::IntervalSet::intervalSetIntersection(const IntervalSet& second_set) const {

  IntervalSet set_intersection;

  if (this->empty() or second_set.empty()) {

    return set_intersection;

  }

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

      iter = std::ranges::next(iter, 1, this->end());

    } else {

      second_iter = std::ranges::next(second_iter, 1, this->end());

    }

  }

  return set_intersection;

}

// Returns the interval unions between this set and the argument set.
// Intervals are modified/extended as necessary for a disjoint interval minimal union.
kel::IntervalSet kel::IntervalSet::intervalSetUnion(const IntervalSet& interval_set) const {

  IntervalSet set_union;
  IntervalMultiSet sorted_intervals;

  if (this->empty() and interval_set.empty()) {

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

          }

          active_interval = *iter;

        }

      }

      iter = std::ranges::next(iter, 1, sorted_intervals.end());

    } // Loop through sorted intervals.

    auto [set_iter, result] = set_union.insert(active_interval);
    if (not result) {

      ExecEnv::log().warn("IntervalSet::intervalSetUnion; Cannot add duplicate interval: [{}, {})",
                          active_interval.lower(), active_interval.upper());

    }

  } // Empty condition.

  return set_union;

}


// Find the set intervals that intersect the argument interval.
[[nodiscard]] std::vector<kel::OpenRightInterval> kel::IntervalSet::findIntersectsInterval(const OpenRightInterval &interval) const {

  std::vector<OpenRightInterval> value_vector;

  if (this->empty()) {

    return value_vector;

  }

  auto const bound_iter = this->lower_bound(interval); // iter->first.lower() >= arg.lower()

  // Scroll forward to check intervals that intersect
  auto next_iter = bound_iter;
  if (next_iter != this->end()) {

    while (next_iter->lower() < interval.upper()) {

      auto const [lower, upper] = next_iter->intersection(interval);
      if (lower < upper) {

        value_vector.push_back(*next_iter);

      }

      next_iter = std::ranges::next(next_iter, 1, this->end());
      if (next_iter == this->end()) {

        break;

      }

    }

  }

  // Scroll backwards to find intervals that intersect.
  auto prev_iter = std::ranges::prev(bound_iter, 1, this->begin());
  // Check that bound_iter was not this->begin()
  if (prev_iter != bound_iter) {

    while (interval.lower() < prev_iter->upper()) {

      // Is there an intersection between then argument interval and the previous interval.
      auto const [lower, upper] = prev_iter->intersection(interval);
      // If so, then return the previous interval.
      if (lower < upper) {

        value_vector.push_back(*prev_iter);

      }

      if (prev_iter == this->begin()) {

        break;

      }
      prev_iter = std::ranges::prev(prev_iter, 1 , this->begin());

    } // while

  } // Not end()

  // Else just return the upper interval (which may be end())
  return value_vector;

}

