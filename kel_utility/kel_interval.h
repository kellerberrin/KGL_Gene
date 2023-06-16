//
// Created by kellerberrin on 15/06/23.
//

#ifndef KEL_INTERVAL_H
#define KEL_INTERVAL_H


#include <set>
#include <map>
#include <string>


namespace kellerberrin {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Defines a simple right open interval [lower_, upper_)
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


class OpenRightInterval {

public:

  OpenRightInterval(size_t lower, size_t upper);
  ~OpenRightInterval() = default;

  OpenRightInterval(const OpenRightInterval &copy) = default;
  OpenRightInterval &operator=(const OpenRightInterval &copy) = default;

  [[nodiscard]] size_t lower() const { return lower_; }
  [[nodiscard]] size_t upper() const { return upper_; }
  [[nodiscard]] size_t size() const { return upper_ - lower_; }

  [[nodiscard]] bool containsOffset(size_t offset) const { return offset >= lower_ and offset < upper_; }
  [[nodiscard]] bool containsInterval(const OpenRightInterval &interval) const { return interval.lower_ >= lower_ and (interval.lower_ + interval.size()) <= upper_; }

private:

  size_t upper_{0};
  size_t lower_{0};

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Adapters for std::set and std::map to use OpenRightInterval.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Comparison operator used to order intervals within indexed containers std::set or std::map.
struct CompareInterval {

  bool operator()(const OpenRightInterval &lhs, const OpenRightInterval &rhs) const {
    return lhs.lower() < rhs.lower();
  }

};

// Interval adapter for std::set.
class IntervalSet : public std::set<OpenRightInterval, CompareInterval> {

public:

  IntervalSet() = default;

  ~IntervalSet() = default;

  // Find the interval that contains the argument OR the interval immediately greater (lower > arg.lower) than the argument interval.
  [[nodiscard]] auto findUpperEqualIter(const OpenRightInterval &interval) const {

    auto iter = this->lower_bound(interval);
    // Return the previous map interval if it contains the argument interval.
    auto prev_iter = std::prev(iter, 1);
    if (prev_iter != this->end()) {

      auto const &interval_key = *prev_iter;
      if (interval_key.containsInterval(interval)) {

        return prev_iter;

      }

    }

    // Else just return the upper interval
    return iter;

  }

  [[nodiscard]] bool containsInterval(const OpenRightInterval &interval) const;
  [[nodiscard]] bool containsOffset(size_t offset) const { return contains(OpenRightInterval(offset, offset + 1)); }

};

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




} // namespace





#endif //KEL_INTERVAL_H
