//
// Created by kellerberrin on 17/06/23.
//

#ifndef KEL_INTERVAL_MAP_H
#define KEL_INTERVAL_MAP_H

#include "kel_interval.h"

namespace kellerberrin {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Interval adapter for std::map.
template<typename ValueType>
using IntervalMapType = std::map<OpenRightInterval, ValueType, CompareInterval>;

template<typename ValueType>
class IntervalMap : public IntervalMapType<ValueType> {

public:

  IntervalMap() = default;
  ~IntervalMap() = default;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Interval adapter for std::multimap.
template<typename ValueType>
using IntervalMultiMapType = std::multimap<OpenRightInterval, ValueType, CompareInterval>;

template<typename ValueType>
class IntervalMultiMap : public IntervalMultiMapType<ValueType> {

public:

  IntervalMultiMap() = default;
  ~IntervalMultiMap() = default;

  // Returns vector the value for all map intervals that contain the argument
  [[nodiscard]] std::vector<ValueType> findIntervalContains(const OpenRightInterval &interval) const {

    std::vector<ValueType> value_vector;

    if (this->empty()) {

      return value_vector;

    }

    auto iter = this->lower_bound(interval);

    if (iter != this->end()) {

      auto const &[interval_key, value] = *iter;
      if (interval_key.containsInterval(interval)) {

        value_vector.push_back(value);

      }

    }

    // Look at previous intervals.
    auto prev_iter = std::ranges::prev(iter, 1, this->begin());
    // Check that iter was not at this->begin()
    if (prev_iter != iter) {

      while(interval.lower() < prev_iter->first.upper()) {

        auto const &[interval_key, value] = *prev_iter;
        if (interval_key.containsInterval(interval)) {

          value_vector.push_back(value);

        }

        if (prev_iter == this->begin()) {

          break;

        }
        prev_iter = std::ranges::prev(prev_iter, 1, this->begin());

      }

    }

    return value_vector;

  }

  // Find the interval that contains the argument OR the interval immediately greater (lower > arg.lower) than the argument interval.
  [[nodiscard]] std::vector<ValueType> findIntersectsIntervals(const OpenRightInterval &interval) const {

    std::vector<ValueType> value_vector;

    if (this->empty()) {

      return value_vector;

    }

    auto const bound_iter = this->lower_bound(interval); // iter->first.lower() >= arg.lower()

    // Scroll forward to check intervals that intersect
    auto next_iter = bound_iter;
    if (next_iter != this->end()) {

      while (next_iter->first.lower() < interval.upper()) {

        auto const [lower, upper] = next_iter->first.intersection(interval);
        if (lower < upper) {

          value_vector.push_back(next_iter->second);

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

      while (interval.lower() < prev_iter->first.upper()) {

        // Is there an intersection between then argument interval and the previous interval.
        auto const [lower, upper] = prev_iter->first.intersection(interval);
        // If so, then return the previous interval.
        if (lower < upper) {

          value_vector.push_back(prev_iter->second);

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

  [[nodiscard]] std::vector<std::tuple<OpenRightInterval, ValueType, ValueType>> intersect() const {

    std::vector<std::tuple<OpenRightInterval, ValueType, ValueType>> intersect_vec;

    if (this->empty()) {

      return intersect_vec;

    }

    auto iter = this->begin();
    while (iter != this->end()) {

      auto next_iter = std::ranges::next(iter, 1, this->end());
      while (next_iter != this->end()) {

        auto const& [interval, value] = *iter;
        auto const& [next_interval, next_value] = *next_iter;
        auto const [lower, upper] = interval.intersection(next_interval);
        if (lower < upper) {

          intersect_vec.push_back({OpenRightInterval(lower, upper), value, next_value});

        } else {

          break;

        }

        next_iter = std::ranges::next(next_iter, 1, this->end());

      }

      iter = std::ranges::next(iter, 1, this->end());

    }

    return intersect_vec;

  }

  [[nodiscard]] bool containsInterval(const OpenRightInterval &interval) const { return not findIntervalContains(interval).empty(); }

};


} // namespace



#endif //KEL_INTERVAL_MAP_H
