//
// Created by kellerberrin on 21/09/23.
//

#ifndef KEL_INTERVAL_TYPE_H
#define KEL_INTERVAL_TYPE_H


#include <set>
#include <map>
#include <vector>
#include <string>

#include "kel_exec_env.h"


namespace kellerberrin {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Define the concepts for the underlying numeric types used in the OpenRightInterval template class.
// This includes floats. However, floats may require work on the comparison operators.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
concept IntervalNumericType = std::integral<T> || std::floating_point<T>;

template < class T >
concept SignedIntervalType = (std::integral<T> && std::is_signed_v<T>) || std::floating_point<T>;

template < class T >
concept UnsignedIntervalType = std::integral<T> && !std::signed_integral<T>;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Defines a simple right open interval [lower_, upper_)
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <IntervalNumericType IntervalValue, SignedIntervalType SignedInterval>
class OpenRightInterval {

public:

  OpenRightInterval(IntervalValue lower, IntervalValue upper) { resize(lower, upper); }
  ~OpenRightInterval() = default;

  OpenRightInterval(const OpenRightInterval &copy) = default;

  OpenRightInterval &operator=(const OpenRightInterval &copy) = default;

  void resize(IntervalValue lower, IntervalValue upper)  {

    if (upper < lower) {

      ExecEnv::log().warn("OpenRightInterval::resize; Incorrect Initialization, Upper Offset: {} < Lower Offset: {}", upper, lower);
      std::swap(lower, upper);

    }

    lower_ = lower;
    upper_ = upper;

  }


  // Shift the interval without changing it's size.
  [[nodiscard]] OpenRightInterval translate(SignedInterval shift) const {

    OpenRightInterval translated(*this);

    if constexpr (UnsignedIntervalType<IntervalValue>) {

      if ((static_cast<SignedInterval>(lower()) + shift) < 0) {

        ExecEnv::log().warn("OpenRightInterval::translate; translate shift: {} results in negative values for interval: {}", shift, toString());
        return {lower(), upper()};

      }

    }

    IntervalValue translate_lower = static_cast<SignedInterval>(lower()) + shift;
    IntervalValue translate_upper = static_cast<SignedInterval>(upper()) + shift;

    return { translate_lower , translate_upper };

  }

  // Return the zero-translated interval so that lower() == 0.
  [[nodiscard]] OpenRightInterval translateZero() const {

    SignedInterval shift =  -1 * static_cast<SignedInterval>(lower());
    return translate(shift);

  }

  bool operator==(const OpenRightInterval &rhs) const { return this->lower() == rhs.lower() and this->upper() == rhs.upper(); }

  bool operator<(const OpenRightInterval &rhs) const { return this->lower() < rhs.lower(); }

  [[nodiscard]] IntervalValue lower() const { return lower_; }

  [[nodiscard]] IntervalValue upper() const { return upper_; }

  [[nodiscard]] size_t size() const { return upper_ - lower_; }

  // Returns the intersection interval or the empty [0, 0) interval indicating no intersection.
  [[nodiscard]] OpenRightInterval intersection(const OpenRightInterval &interval) const {

    if (lower_ >= interval.upper_ or interval.lower_ >= upper_) {

      return {0, 0};

    }

    return { std::max<IntervalValue>(lower_, interval.lower_), std::min<IntervalValue>(upper_, interval.upper_ ) };

  }

  // Merge intersecting or adjacent intervals. If the argument intervals are disjoint and not adjacent
  // then the empty [0, 0) interval is returned.
  // Note that the merging of the empty intervals will also produce an empty interval.
  [[nodiscard]] OpenRightInterval merge(const OpenRightInterval &interval) const {

    if (intersects(interval) or adjacent(interval)) {

      return { std::min<IntervalValue>(lower_, interval.lower_), std::max<IntervalValue>(upper_, interval.upper_ ) };

    }

    return {0, 0};

  }

  [[nodiscard]] bool empty() const { return size() == 0; }

  [[nodiscard]] bool containsOffset(size_t offset) const { return offset >= lower_ and offset < upper_; }

  [[nodiscard]] bool containsInterval(const OpenRightInterval &interval) const { return intersection(interval) == interval; }

  // Note that empty [0, 0) intervals adjoin each other, other empty intervals such as [k, k) and [l, l) do not adjoin if k != l.
  [[nodiscard]] bool adjacent(const OpenRightInterval &interval) const { return lower_ == interval.upper_ or interval.lower_ == upper_; }

  [[nodiscard]] bool intersects(const OpenRightInterval &interval) const { return not disjoint(interval); }

  [[nodiscard]] bool disjoint(const OpenRightInterval &interval) const { return intersection(interval).empty(); }

  // Convenience routine to convert an interval to a string.
  [[nodiscard]] std::string toString() const { return "[ " + std::to_string(lower_) + ", " + std::to_string(upper_) + ")"; }

private:

  IntervalValue upper_{0};
  IntervalValue lower_{0};

};


} // namespace


#endif //KEL_INTERVAL_TYPE_H
