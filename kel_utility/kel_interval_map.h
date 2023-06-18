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

  // Returns any the map interval that contains the argument or end().
  [[nodiscard]] auto findIntervalIter(const OpenRightInterval &interval) const {

    auto iter = this->lower_bound(interval);
    if (iter != this->end()) {

      // Do the lower bound of the intervals match.
      auto const &[interval_key, value] = *iter;
      if (interval_key.lower() == interval.lower()) {

        if (interval_key.containsInterval(interval)) {

          return iter;

        } else {

          return this->end();

        }

      }

    }

    // Look at the previous interval.
    iter = std::prev(iter, 1);
    if (iter != this->end()) {

      auto const &[interval_key, value] = *iter;
      if (interval_key.containsInterval(interval)) {

        return iter;

      } else {

        return this->end();

      }

    }

    return this->end();

  }

  // Find the interval that contains the argument OR the interval immediately greater (lower > arg.lower) than the argument interval.
  [[nodiscard]] auto findUpperEqualIter(const OpenRightInterval &interval) const {

    auto iter = this->lower_bound(interval);
    // Return the previous map interval if it contains the argument interval.
    auto prev_iter = std::prev(iter, 1);
    if (prev_iter != this->end()) {

      auto const &[interval_key, value] = *prev_iter;
      if (interval_key.containsInterval(interval)) {

        return prev_iter;

      }

    }

    // Else just return the upper interval
    return iter;

  }

  [[nodiscard]] auto findUpperOffsetIter(size_t offset) const { return findUpperEqualIter(OpenRightInterval(offset, offset + 1)); }
  [[nodiscard]] bool containsInterval(const OpenRightInterval &interval) const { return findIntervalIter(interval) != this->end(); }

};


// Interval adapter for std::multimap.
template<typename ValueType>
using IntervalMultiMapType = std::multimap<OpenRightInterval, ValueType, CompareInterval>;

template<typename ValueType>
class IntervalMultiMap : public IntervalMultiMapType<ValueType> {

public:

  IntervalMultiMap() = default;
  ~IntervalMultiMap() = default;

  // Returns any the map interval that contains the argument or end().
  [[nodiscard]] auto findIntervalIter(const OpenRightInterval &interval) const {

    auto iter = this->lower_bound(interval);
    if (iter != this->end()) {

      // Do the lower bound of the intervals match.
      auto const &[interval_key, value] = *iter;
      if (interval_key.lower() == interval.lower()) {

        if (interval_key.containsInterval(interval)) {

          return iter;

        } else {

          return this->end();

        }

      }

    }

    // Look at the previous interval.
    iter = std::prev(iter, 1);
    if (iter != this->end()) {

      auto const &[interval_key, value] = *iter;
      if (interval_key.containsInterval(interval)) {

        return iter;

      } else {

        return this->end();

      }

    }

    return this->end();

  }

  // Find the interval that contains the argument OR the interval immediately greater (lower > arg.lower) than the argument interval.
  [[nodiscard]] auto findUpperEqualIter(const OpenRightInterval &interval) const {

    auto iter = this->lower_bound(interval);
    // Return the previous map interval if it contains the argument interval.
    auto prev_iter = std::prev(iter, 1);
    if (prev_iter != this->end()) {

      auto const &[interval_key, value] = *prev_iter;
      if (interval_key.containsInterval(interval)) {

        return prev_iter;

      }

    }

    // Else just return the upper interval
    return iter;

  }

  [[nodiscard]] std::vector<std::tuple<OpenRightInterval, ValueType, ValueType>> intersect() const {

    std::vector<std::tuple<OpenRightInterval, ValueType, ValueType>> intersect_vec;

    auto iter = this->begin();
    while (iter != this->end()) {

      auto next_iter = std::next(iter, 1);
      while (next_iter != this->end()) {

        auto const& [interval, value] = *iter;
        auto const& [next_interval, next_value] = *next_iter;
        auto const [lower, upper] = interval.intersection(next_interval);
        if (lower != upper) {

          intersect_vec.push_back({OpenRightInterval(lower, upper), value, next_value});

        } else {

          break;

        }

        ++next_iter;

      }

      ++iter;

    }

    return intersect_vec;

  }


  [[nodiscard]] auto findUpperOffsetIter(size_t offset) const { return findUpperEqualIter(OpenRightInterval(offset, offset + 1)); }
  [[nodiscard]] bool containsInterval(const OpenRightInterval &interval) const { return findIntervalIter(interval) != this->end(); }

};


} // namespace



#endif //KEL_INTERVAL_MAP_H
