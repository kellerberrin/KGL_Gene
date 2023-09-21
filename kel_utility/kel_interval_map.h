//
// Created by kellerberrin on 17/06/23.
//

#ifndef KEL_INTERVAL_MAP_H
#define KEL_INTERVAL_MAP_H

#include "kel_interval_unsigned.h"

namespace kellerberrin {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Interval adapter for std::map.
template<typename ValueType>
using IntervalMapType = std::map<OpenRightUnsigned, ValueType, CompareUnsignedIntervalLower>;

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
using IntervalLowerMultiMapType = std::multimap<OpenRightUnsigned, ValueType, CompareUnsignedIntervalLower>;

template<typename ValueType>
class IntervalLowerMultiMap : public IntervalLowerMultiMapType<ValueType> {

public:

  IntervalLowerMultiMap() = default;
  ~IntervalLowerMultiMap() = default;

  // Returns vector the value for all map intervals that contain the argument
  [[nodiscard]] std::vector<ValueType> findIntervalContains(const OpenRightUnsigned &interval) const {

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
  [[nodiscard]] std::vector<ValueType> findIntersectsIntervals(const OpenRightUnsigned &interval) const {

    std::vector<ValueType> value_vector;

    if (this->empty()) {

      return value_vector;

    }

    auto const bound_iter = this->lower_bound(interval); // iter->first.lower() >= arg.lower()

    // Scroll forward to check intervals that intersect
    auto next_iter = bound_iter;
    if (next_iter != this->end()) {

      while (next_iter->first.lower() < interval.upper()) {

        auto const intersect_interval = next_iter->first.intersection(interval);
        if (not intersect_interval.empty()) {

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
        auto const intersect_interval = prev_iter->first.intersection(interval);
        // If so, then return the previous interval.
        if (not intersect_interval.empty()) {

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

  [[nodiscard]] bool containsInterval(const OpenRightUnsigned &interval) const { return not findIntervalContains(interval).empty(); }

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Interval adapter for std::map.
template<typename ValueType>
using IntervalUpperMapType = std::map<OpenRightUnsigned, ValueType, CompareUnsignedIntervalUpper>;

template<typename ValueType>
class IntervalUpperMap : public IntervalUpperMapType<ValueType> {

public:

  IntervalUpperMap() = default;
  ~IntervalUpperMap() = default;


};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Interval adapter for std::multimap.
template<typename ValueType>
using IntervalUpperMultiMapType = std::multimap<OpenRightUnsigned, ValueType, CompareUnsignedIntervalUpper>;


} // namespace



#endif //KEL_INTERVAL_MAP_H
